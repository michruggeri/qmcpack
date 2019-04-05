#include "Interfaces/ESHDF5/ESHDF5Interface.h"
#include <fftw3.h>
#include <Utilities/ProgressReportEngine.h>
//#include <QMCWaveFunctions/einspline_helper.hpp>
#include "ParticleIO/HDFParticleAttrib.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "simd/vmath.hpp"
#include "qmc_common.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
namespace qmcplusplus
{

void ESHDF5Interface::initialize()
{
   std::cerr << "Beginning the initialization... \n";
   H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
   app_log()<<"Opened HDF5 File for reading... " << H5FileName.c_str() << "\n";
   if (H5FileID < 0)
   {
     app_error() << "Could not open HDF5 file \"" << H5FileName
                << "\" in ESHDF5Interface::initialize().  Aborting.\n";
     APP_ABORT("ESHDF5Interface::initialize()");
   }
}
bool ESHDF5Interface::put(xmlNodePtr cur)
{
  OhmmsAttributeSet a;
  a.add (H5FileName, "href");
  a.put (cur);
  return true;
}

void ESHDF5Interface::getVersion()
{
  HDFAttribIO<TinyVector<int,3> > h_version(Version);
  h_version.read (H5FileID, "/version");
  app_log() << "  USING INTERFACE:  ESHDF orbital file version "
            << Version[0] << "." << Version[1] <<"." << Version[2] << std::endl;
}

void ESHDF5Interface::getPrimVecs(Tensor<double,OHMMS_DIM> & Lattice)
{
  HDFAttribIO<Tensor<double,OHMMS_DIM> > h_Lattice(Lattice);
  h_Lattice.read      (H5FileID, "/supercell/primitive_vectors");
  return;
}
int ESHDF5Interface::getNumBands()
{
  int numBands(0);
  HDFAttribIO<int> h_NumBands(numBands);
  h_NumBands.read(H5FileID, "/electrons/kpoint_0/spin_0/number_of_states");
  return numBands;
}
void ESHDF5Interface::getNumElectrons(std::vector<int>& num_spins)
{
  HDFAttribIO<std::vector<int> > h_numspins(num_spins);
  h_numspins.read(H5FileID,"/electrons/number_of_electrons");
  return;
}
int ESHDF5Interface::getNumSpins()
{
  int NumSpins(0);
  HDFAttribIO<int> h_NumSpins(NumSpins); 
  h_NumSpins.read(H5FileID, "/electrons/number_of_spins");
  return NumSpins;
}
int ESHDF5Interface::getNumTwists()
{
  int NumTwists(0);
  HDFAttribIO<int> h_NumTwists(NumTwists);
  h_NumTwists.read(H5FileID, "/electrons/number_of_kpoints");
  return NumTwists;
}
int ESHDF5Interface::getNumCoreStates()
{
  int NumCore(0);
  HDFAttribIO<int> h_NumCore(NumCore); 
  h_NumCore.read(H5FileID, "/electrons/kpoint_0/spin_0/number_of_core_states");
  return NumCore;
}
int ESHDF5Interface::getNumMuffinTins()
{
  int NumMuffinTins(0);
  HDFAttribIO<int> h_NumMuffinTins(NumMuffinTins);
  h_NumMuffinTins.read(H5FileID, "/muffin_tins/number_of_tins");
  return NumMuffinTins;
}
int ESHDF5Interface::getHaveDPsi()
{
  int have_dpsi(0);
  HDFAttribIO<int> h_have_dpsi(have_dpsi); 
  h_have_dpsi.read(H5FileID, "/electrons/have_dpsi");
  return have_dpsi;
}
int ESHDF5Interface::getNumAtomicOrbitals()  
{
  int NumAtomicOrbitals(0);
  HDFAttribIO<int>  h_NumAtomicOrbitals(NumAtomicOrbitals);
  h_NumAtomicOrbitals.read(H5FileID, "/electrons/number_of_atomic_orbitals");
  return NumAtomicOrbitals;
}

int ESHDF5Interface::getNumSpecies()
{
  int num_species(0);
  HDFAttribIO<int> h_num_species(num_species);
  h_num_species.read (H5FileID, "/atoms/number_of_species");
  return num_species;
}

int ESHDF5Interface::getNumAtoms()
{
  int num_atoms(0);
  HDFAttribIO<int> h_num_atoms(num_atoms);
  h_num_atoms.read (H5FileID, "/atoms/number_of_atoms");
  return num_atoms;
}

//Warning.  This will not resize 
void ESHDF5Interface::getSpeciesIDs(ParticleIndex_t& species_ids)
{
  int natom = getNumAtoms();
  species_ids.resize(natom);
  HDFAttribIO<ParticleIndex_t> h_species_ids(species_ids);
  h_species_ids.read (H5FileID, "/atoms/species_ids");
  return;
}

void ESHDF5Interface::getSpeciesData(Vector<int>& am, Vector<int>& charge, Vector<double>& mass, Vector<std::string>& names)
{
  int num_species=getNumSpecies();
  am.resize(num_species);
  charge.resize(num_species);
  mass.resize(num_species);
  names.resize(num_species);

  for (int isp=0; isp<num_species; isp++)
  {
    std::ostringstream anumname;
    std::ostringstream chargename;
    std::ostringstream massname;
    std::ostringstream speciesname;

    anumname    << "/atoms/species_" << isp << "/atomic_number";
    chargename  << "/atoms/species_" << isp << "/valence_charge";
    massname    << "/atoms/species_" << isp << "/mass";
    speciesname << "/atoms/species_" << isp << "/name";

    HDFAttribIO<int> h_atomic_number (am[isp]);
    HDFAttribIO<int> h_charge_number (charge[isp]);
    HDFAttribIO<double> h_mass_number (mass[isp]);
    HDFAttribIO<std::string> h_species_name (names[isp]);
    
    h_atomic_number.read(H5FileID, anumname.str().c_str());
    h_charge_number.read(H5FileID, chargename.str().c_str());
    h_mass_number.read(H5FileID, massname.str().c_str());
    h_species_name.read(H5FileID, speciesname.str().c_str());
    mass[isp]=1.0;
  }

}

void ESHDF5Interface::getAtomicNumbers(Vector<int> & am)
{
  int num_species=getNumSpecies();
  am.resize(num_species);
  for (int isp=0; isp<num_species; isp++)
  {
    std::ostringstream name;
    name << "/atoms/species_" << isp << "/atomic_number";
    HDFAttribIO<int> h_atomic_number (am[isp]);
    h_atomic_number.read(H5FileID, name.str().c_str());
  }

  return;
}

void ESHDF5Interface::getIonPositions( ParticleSet::ParticlePos_t & R)
{
  int natom=getNumAtoms();
  R.resize(natom);
  HDFAttribIO<ParticleSet::ParticlePos_t> h_IonPos(R);
  h_IonPos.read   (H5FileID, "/atoms/positions");

  return;
}

void ESHDF5Interface::getAtomicOrbitals(std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals)
{
  int NumAtomicOrbitals=getNumAtomicOrbitals();

  AtomicOrbitals.resize(NumAtomicOrbitals);
  for (int iat=0; iat<NumAtomicOrbitals; iat++)
  {
    AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
    int lmax, polynomial_order, spline_points;
    double cutoff_radius, polynomial_radius, spline_radius;
    PosType position;
    HDFAttribIO<int> h_lmax(lmax), h_polynomial_order(polynomial_order),
                h_spline_points(spline_points);
    HDFAttribIO<double> h_cutoff_radius(cutoff_radius),
                h_polynomial_radius(polynomial_radius),
                h_spline_radius(spline_radius);
    HDFAttribIO<PosType> h_position(position);
    std::ostringstream groupstream;
    groupstream << "/electrons/atomic_orbital_" << iat << "/";
    std::string groupname = groupstream.str();
    h_lmax.read              (H5FileID, (groupname + "lmax"             ).c_str());
    h_polynomial_order.read  (H5FileID, (groupname + "polynomial_order" ).c_str());
    h_spline_points.read     (H5FileID, (groupname + "spline_points"    ).c_str());
    h_cutoff_radius.read     (H5FileID, (groupname + "cutoff_radius"    ).c_str());
    h_polynomial_radius.read (H5FileID, (groupname + "polynomial_radius").c_str());
    h_spline_radius.read     (H5FileID, (groupname + "spline_radius"    ).c_str());
    h_position.read          (H5FileID, (groupname + "position"         ).c_str());
    orb.set_pos (position);
    orb.set_lmax (lmax);
    orb.set_cutoff (cutoff_radius);
    orb.set_spline (spline_radius, spline_points);
    orb.set_polynomial (polynomial_radius, polynomial_order);
  }

  return;
}

void ESHDF5Interface::getOrbEigenvals(const int spin, const int kid, std::vector<double>& eigenvals)
{
    std::ostringstream ePath;
    ePath << "/electrons/kpoint_" << kid << "/spin_"
          << spin << "/eigenvalues";
    Vector<double> eigvals;
    HDFAttribIO<std::vector<double> > h_eigvals(eigenvals);
    std::cerr << "adadasdasdasda\n";
    app_log() << ePath.str() << "\n";
    h_eigvals.read(H5FileID, ePath.str().c_str());
    std::cerr << "qweqweqweqeqwe\n";

    return;
}

void ESHDF5Interface::getTwistData(std::vector<PosType>& TwistAngles,
		      std::vector<double>& TwistWeight, 
 		      std::vector<int>& TwistSymmetry)
{
  int NumTwists = getNumTwists();
  
  TwistAngles.resize(NumTwists);
  TwistSymmetry.resize(NumTwists);
  TwistWeight.resize(NumTwists);
  for (int ti=0; ti<NumTwists; ti++)
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/reduced_k";
    HDFAttribIO<PosType> h_Twist(TwistAngles[ti]);
    h_Twist.read (H5FileID, path.str().c_str());
    if ((Version[0] >= 2) and (Version[1] >= 1))
    {
      std::ostringstream sym_path;
      sym_path << "/electrons/kpoint_" << ti << "/symgroup";
      HDFAttribIO<int> h_Sym(TwistSymmetry[ti]);
      h_Sym.read (H5FileID, sym_path.str().c_str());
      std::ostringstream nsym_path;
      nsym_path << "/electrons/kpoint_" << ti << "/numsym";
      HDFAttribIO<double> h_Nsym(TwistWeight[ti]);
      h_Nsym.read (H5FileID, nsym_path.str().c_str());
    }
    // Early versions from wfconv had wrong sign convention for
    // k-points.  EinsplineSet uses the opposite sign convention
    // from most DFT codes.
    if (Version[0] >= 2)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        TwistAngles[ti][dim] *= -1.0;
    //       snprintf (buff, 1000, "  Found twist angle (%6.3f, %6.3f, %6.3f)\n",
    //           TwistAngles[ti][0], TwistAngles[ti][1], TwistAngles[ti][2]);
    //       app_log() << buff;
  }
  return;
}
int ESHDF5Interface::getHasPsiG()
{
    int hasPsig(1);
#if defined(__bgp__)||(__bgq__)
    hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_g"); //Ray:  HWI
    if(gid<0)   
      hasPsig=0;
    H5Dclose(gid);
    return hasPsig;
#else
    TinyVector<int,3> MeshSize(0);
    HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize); //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");  //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/mesh"); //Ray:  HWI
    hasPsig = (MeshSize[0] == 0);
    return hasPsig;
#endif

}
    
void ESHDF5Interface::getMeshSize(TinyVector<int,3> & mesh)
{
#if defined(__bgp__)||(__bgq__)
#else
    HDFAttribIO<TinyVector<int,3> > h_mesh(mesh); //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");  //Ray:  HWI
    h_mesh.read (H5FileID, "/electrons/mesh"); //Ray:  HWI
    return;
#endif
}
     
void ESHDF5Interface::getReducedGVecs(std::vector<std::vector<TinyVector<int,3> > > & gvecs, int index)
{
    std::ostringstream Gpath;  //Ray:  HWI
    Gpath    << "/electrons/kpoint_"<<index<<"/gvectors";  //Ray:  HWI 
    HDFAttribIO<std::vector<TinyVector<int,3> > > h_Gvecs(gvecs[index]); //Ray:  HWI
    h_Gvecs.read (H5FileID, Gpath.str().c_str()); //Ray:  HWI

    return;
}

//bool ESHDF5Interface::getPsi_kspace(Vector<std::complex<double> > & Cg, int orbid, int twistid)
//{
//        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex); //Ray:  HWI (Handle with interface)
//        return h5f.read(cG,s); //RAY:  HWI  (Handle with interface).
//}

bool ESHDF5Interface::getPsi_kspace(Vector<std::complex<double> > & cG,int spin, int orbid, int twistid)
{
        std::string s=psi_g_path(twistid,spin,orbid); //Ray:  HWI (Handle with interface)
        HDFAttribIO<Vector<std::complex<double> > > h_cg(cG);
	h_cg.read(H5FileID, s.c_str());
	return true;	
//        return h5f.read(cG,s); //RAY:  HWI  (Handle with interface).
}

}
