import ast
import h5py
import numpy
import scipy.linalg
from pyscf import fci

def write_wfn_mol(scf_data, ortho_ao, filename, wfn=None):
    """Generate QMCPACK trial wavefunction.

    Parameters
    ----------
    scf_data : dict
        Dictionary containing scf data extracted from pyscf checkpoint file.
    ortho_ao : bool
        Whether we are working in orthogonalised AO basis or not.
    filename : string
        HDF5 file path to store wavefunction to.
    wfn : tuple
        User defined wavefunction. Not fully supported. Default None.

    Returns
    -------
    wfn : :class:`numpy.ndarray`
        Wavefunction as numpy array. Format depends on wavefunction.
    """
    ghf = False
    mol = scf_data['mol']
    nalpha, nbeta = mol.nelec
    C = scf_data['mo_coeff']
    X = scf_data['X']
    uhf = scf_data['isUHF']
    # For RHF only nalpha entries will be filled.
    if uhf:
        norb = C[0].shape[0]
    else:
        norb = C.shape[0]
    if wfn is None:
        wfn = numpy.zeros((norb,nalpha+nbeta), dtype=numpy.float64)
        wfn_type = 'NOMSD'
        if ortho_ao:
            Xinv = scipy.linalg.inv(X)
            if uhf:
                # We are assuming C matrix is energy ordered.
                wfn[:,:nalpha] = numpy.dot(Xinv, C[0])[:,:nalpha]
                wfn[:,nalpha:] = numpy.dot(Xinv, C[1])[:,:nbeta]
            else:
                wfn[:,:nalpha] = numpy.dot(Xinv, C)[:,:nalpha]
        else:
            # Assuming we are working in MO basis, only works for RHF, ROHF trials.
            I = numpy.identity(C.shape[-1], dtype=numpy.float64)
            wfn[:,:nalpha] = I[:,:nalpha]
            if uhf:
                print(" # Warning: UHF trial wavefunction can only be used of "
                      "working in ortho AO basis.")
    else:
        # Multi determinant style trial wavefunction.
        coeffs, cis = wfn
        wfn_type = 'PHMSD'
    with h5py.File(filename, 'r+') as fh5:
        # TODO: FIX for GHF eventually.
        if ghf:
            walker_type = 'NONCOLLINEAR'
        elif uhf:
            walker_type = 'COLLINEAR'
        else:
            walker_type = 'CLOSED'
        if wfn_type == 'PHMSD':
            walker_type = 'COLLINEAR'
        fh5['Wavefunction/type'] = wfn_type
        fh5['Wavefunction/walker_type'] = walker_type
        # TODO: Fix for multideterminant case.
        if wfn_type == 'PHMSD':
            fh5['Wavefunction/occ_a'] = cis[0]
            fh5['Wavefunction/occ_b'] = cis[1]
            fh5['Wavefunction/coeffs'] = coeffs
            nci = len(cis)
            dims = [nci, 0, len(cis[0])]
        else:
            fh5['Wavefunction/orbs'] = wfn
            write_nomsd_wfn('wfn.dat', wfn, nalpha, uhf)
            dims = [1, wfn.shape[0], wfn.shape[1]]
        fh5['Wavefunction/dims'] = numpy.array(dims, dtype=numpy.int32)

    return wfn

#
# Graveyard. Old QMCPACK wavefunction plain text format.
# Keep around for backwards compatability.
#
def write_nomsd_wfn(filename, wfn, nalpha, uhf, coeffs=[1.0]):
    if len(wfn.shape) == 2:
        wfn = wfn.reshape((1,wfn.shape[0],wfn.shape[1]))
    namelist = qmcpack_wfn_namelist(wfn.shape[0], uhf)
    with open(filename, 'w') as f:
        f.write(namelist)
        f.write('Coefficients: ' + ' '.join(str(c) for c in coeffs) +'\n')
        for (i,d) in enumerate(wfn):
            f.write('Determinant: {}\n'.format(i+1))
            if uhf:
                write_single(f, d[:,:nalpha])
                write_single(f, d[:,nalpha:])
            else:
                write_single(f, d[:,:nalpha])


def qmcpack_wfn_namelist(nci, uhf):
    return ("&FCI\n UHF = {}\n CMajor\n "
            "NCI = {}\n TYPE = matrix\n/\n".format(int(uhf),nci))

def write_single(out, mos):
    for j in range(0, mos.shape[1]):
        for i in range(0, mos.shape[0]):
            val = mos[i,j]
            out.write('(%.10e,%.10e) '%(val.real, val.imag))
        out.write('\n')

def gen_multi_det_wavefunction(mc, weight_cutoff=0.95, verbose=False,
                               max_ndets=1e5, norb=None,
                               filename=None):
    """Generate multi determinant particle-hole trial wavefunction.

    Format adopted to be compatable with QMCPACK PHMSD type wavefunction.

    Parameters
    ----------
    mc : pyscf CI solver type object
        Input object containing multi determinant coefficients.
    weight_cutoff : float, optional
        Print determinants until accumulated weight equals weight_cutoff.
        Default 0.95.
    verbose : bool
        Print information about process. Default False.
    max_ndets : int
        Max number of determinants to print out. Default 1e5.
    norb : int or None, optional
        Total number of orbitals in simulation. Used if we want to run CI within
        active space but QMC in full space. Deault None.
    filename : string
        Output filename. Default "multi_det.dat"
    """
    occlists = fci.cistring._gen_occslst(range(mc.ncas), mc.nelecas[0])

    ci_coeffs = mc.ci.ravel()
    # Sort coefficients in terms of increasing absolute weight.
    ix_sort = numpy.argsort(numpy.abs(ci_coeffs))[::-1]
    cweight = numpy.cumsum(ci_coeffs[ix_sort]**2)
    max_det = numpy.searchsorted(cweight, weight_cutoff)
    ci_coeffs = ci_coeffs[ix_sort]
    ndets = min(max_det,max_ndets)
    if verbose:
        print(" # Number of dets in CI expansion: {:d}".format(ndets))

    output = open(filename, 'w')
    namelist = "&FCI\n UHF = 0\n NCI = %d\n TYPE = occ\n&END" % ndets
    output.write(namelist+'\n')
    output.write("Configurations:"+'\n')
    if norb is None:
        norb = mc.ncas

    occups = []
    occdns = []
    coeffs = []
    for idet in range(min(max_det,max_ndets)):
        if mc.ncore > 0:
            ocore_up = ' '.join('{:d}'.format(x+1) for x in range(mc.ncore))
            ocore_dn = ' '.join('{:d}'.format(x+1+norb) for x in range(mc.ncore))
        else:
            ocore_up = ' '
            ocore_dn = ' '
        coeff = '%.13f'%ci_coeffs[idet]
        coeffs.append(ci_coeffs[idet])
        ix_alpha = ix_sort[idet] // len(occlists)
        ix_beta = ix_sort[idet] % len(occlists)
        ia = occlists[ix_alpha]
        ib = occlists[ix_beta]
        oup = ' '.join('{:d}'.format(x+1+mc.ncore) for x in ia)
        odown = ' '.join('{:d}'.format(x+norb+1+mc.ncore) for x in ib)
        occups.append([int(o) for o in oup.split()])
        occdns.append([int(o) for o in odown.split()])
        output.write(coeff+' '+ocore_up+' '+oup+' '+ocore_dn+' '+odown+'\n')
    return coeffs, [occups,occdns]


def read_qmcpack_wfn(filename, nskip=9):
    with open(filename) as f:
        content = f.readlines()[nskip:]
    useable = numpy.array([c.split() for c in content]).flatten()
    tuples = [ast.literal_eval(u) for u in useable]
    orbs = [complex(t[0], t[1]) for t in tuples]
    return numpy.array(orbs)

def write_phmsd_wfn(filename, occs, nmo, ncore=0):
    output = open(filename, 'w')
    namelist = "&FCI\n UHF = 0\n NCI = %d\n TYPE = occ\n&END" % len(occs)
    output.write(namelist+'\n')
    output.write("Configurations:"+'\n')
    corea = [i + 1 for i in range(ncore)]
    coreb = [i + nmo + 1 for i in range(ncore)]
    for c, da, db in occs:
        occup = corea + [ncore + oa + 1 for oa in da.tolist()]
        occdn = coreb + [ncore + nmo + ob + 1 for ob in db.tolist()]
        # print(occup, occdn)
        occstra = ' '.join('{:d} '.format(x) for x in occup)
        occstrb = ' '.join('{:d}'.format(x) for x in occdn)
        output.write('%13.8e '%c + occstra + occstrb + '\n')
