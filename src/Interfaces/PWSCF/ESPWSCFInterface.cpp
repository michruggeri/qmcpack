#include "Interfaces/PWSCF/ESPWSCFInterface.h"
#include <fftw3.h>
#include <Utilities/ProgressReportEngine.h>
//#include <QMCWaveFunctions/einspline_helper.hpp>
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "simd/vmath.hpp"
#include "qmc_common.h"
#include "Configuration.h"

#include "Interfaces/PWSCF/pwinterface.h"
#include "mpi.h"
#include <iostream>
namespace qmcplusplus
{

void ESPWSCFInterface::initialize()
{
  app_log()<<"ESPWSCFInterface::initialize().  initialized="<<initialized<<std::endl;
  if (initialized==false)
  {
    //int commID = myComm->getMPI();
    int commID = myComm->rank()+1;
    const char * pwname = PWFileName.c_str();
    app_log()<<" ESPWSCFInterface::initialize() commID = "<<commID<<std::endl;
    app_log()<<" MPI_WORLD_COMM = "<<MPI_COMM_WORLD<<std::endl;
    pwlib_init_(&commID,pwname);
    app_log()<<" pwlib is initialized "<<std::endl;
//  pwlib_setatompos();
//    ParticlePos_t aaa;
//    setIonPositions(aaa);
    pwlib_scf_();
    initialized=true;
  }
 
 return;
   
}
bool ESPWSCFInterface::put(xmlNodePtr cur)
{
  OhmmsAttributeSet a;
  a.add (PWFileName, "href");
  a.put (cur);
  return true;
}

void ESPWSCFInterface::getVersion()
{
  app_log() << "  USING PWINTERFACE:  PW version 5.1.3 ";
}

void ESPWSCFInterface::getPrimVecs(Tensor<double,OHMMS_DIM> & Lattice)
{
 // HDFAttribIO<Tensor<double,OHMMS_DIM> > h_Lattice(Lattice);
 // h_Lattice.read      (H5FileID, "/supercell/primitive_vectors");
  
  pwlib_getbox_data_(Lattice.data());
  return;
}
int ESPWSCFInterface::getNumBands()
{
  
//  std::cout<<"QMCPACKRANK="<<OHMMS::Controller->rank()<<std::endl; 
//  OHMMS::Controller->barrier();
  if (nbands<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
  return nbands;
}

int ESPWSCFInterface::getNumSpins()
{
  int ns=0;
  pwlib_getelectron_info_(&nelec, &nup, &ndown);
  
  if(nup!=ndown and nup>0 and ndown>0)
    return 2;
  else
    return 1; 
}

int ESPWSCFInterface::getNumTwists()
{
  if (nktot<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
  return nktot;
}
int ESPWSCFInterface::getNumAtoms()
{

  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);

  return nat;
}
int ESPWSCFInterface::getNumCoreStates()
{
  return 0;
}
int ESPWSCFInterface::getNumMuffinTins()
{
  return 0;
}
int ESPWSCFInterface::getHaveDPsi()
{
  return 0;
}
int ESPWSCFInterface::getNumAtomicOrbitals()  
{ 
  return 0;
}
void ESPWSCFInterface::getNumElectrons(std::vector<int>& numspins)
{
  numspins.resize(2);
  int ntot_tmp=0;
  pwlib_getelectron_info_(&ntot_tmp, &numspins[0], &numspins[1]);
  return;
}
int ESPWSCFInterface::getNumSpecies()
{
  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);
  return nsp;
}

void ESPWSCFInterface::getSpeciesIDs(ParticleIndex_t& species_ids)
{
  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);
  //double Rtmp[nat][3];
  double **Rtmp;
  Rtmp = new double*[nat];
  for (int i=0;i<nat;i++)
    Rtmp[i] = new double[3]; 
  std::cout << "In ESPWSCFInterface::getSpeciesIDs; number of atoms:\t" << nat << ", number of species:\t" << nsp << std::endl;
  species_ids.resize(nat);
 
  pwlib_getatom_data_(Rtmp[0], &species_ids[0]);
  std::cout << "After pwlib_getatom_data_; number of atoms:\t" << nat << ", number of species:\t" << nsp << std::endl;

  for (int i=0; i<nat; i++){
    species_ids[i]-=1; //this is because pwscf is a fotran code, and so indexing starts at 1 instead of 0.
    std::cout << "atom " << i << "  species " << species_ids[i] << std::endl;
  delete [] Rtmp;
  };
  return;
}

//returns all species data pwscf has on file.  this is am: atomic numbers, q: charges, mass: masses, and species names.  
void ESPWSCFInterface::getSpeciesData(Vector<int> & am, Vector<int> & q, Vector<double>& mass, Vector<std::string>& speciesNames)
{
  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);
  char pwnames[nsp][PW_NAME_LEN]; //pwscf allocates three characters (two name chars terminated by a space)
  
  am.resize(nsp);
  q.resize(nsp);
  mass.resize(nsp);
  speciesNames.resize(nsp);
  pwlib_getspecies_data_(&am[0], &q[0], &mass[0], pwnames[0] );

  speciesNames=convertCharToString(pwnames,nsp);
  return;
}

void ESPWSCFInterface::getAtomicNumbers(Vector<int> & am)
{
  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);
  std::cout << "In getAtomicNumbers, checking nat and nsp (as above)\t" << nat << "   " << nsp << std::endl;
  am.resize(nsp);
  nsp+=7; // This guy and its partner below are here because using nsp for the sizing leads to stack smashing/segfaulting; this will 
           // do for now but it should really be fixed in a proper way at some point -- For my test this is the least increment that avoids crashes
  int vcharg[nsp];
  double mass[nsp];
  char names[nsp*3];
  nsp-=7; // See comment above
  for(int i=0;i<nsp;i++){
    vcharg[i] = -1;
    mass[i] = -0.5;
    am[i]   =3;
  };
  std::cout << "Initialized arrays before pwlib_getspecies_data_\t";
  std::cout << nsp << " " << vcharg[0] << "  " << mass[0] << " " << am[0] << " " <<  names[0] << std::endl;
  pwlib_getspecies_data_(&am[0], vcharg, mass, names);
  std::cout << "After pwlib_getspecies_data_:\t";
  std::cout << nsp << " " << vcharg[0] << "  " << mass[0] << " " << am[0] << " " << names[0] << std::endl;
  return;
}

void ESPWSCFInterface::getIonPositions(ParticlePos_t & R)
{
  
  if (nat==-1 || nsp==-1) pwlib_getatom_info_(&nat, &nsp);
  double Rtmp[nat][OHMMS_DIM];
  int ityp[nat]; 
  
  pwlib_getatom_data_(Rtmp[0], ityp);

  R.resize(nat);

  for (int i=0; i<nat; i++)
  {
    for(int j=0; j<OHMMS_DIM; j++)
    {
       R[i][j]=Rtmp[i][j];
       
    }
  }
  return;
}
//Quantum Espresso uses Bohr internal units.  
void ESPWSCFInterface::setIonPositions(const ParticlePos_t & ionpos)
//void ESPWSCFInterface::setIonPositions(Vector<TinyVector<double,OHMMS_DIM> > & R)
{
  int nat=ionpos.size();
  double Rtmp[nat][OHMMS_DIM];

  for (int i=0; i<nat; i++)
  {
    for (int j=0; j<OHMMS_DIM; j++)
    {
       Rtmp[i][j]=ionpos[i][j];
    }
  } 
  pwlib_setatom_pos_(Rtmp[0]);
  return;
}

void ESPWSCFInterface::getAtomicOrbitals(std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals)
{
  app_log()<<"WARNING!!!  ESPWSCFInterface::getAtomicOrbitals not implemented"<<std::endl;
  return;
}

void ESPWSCFInterface::getOrbEigenvals(const int spin, const int kid, std::vector<double>& eigenvals)
{
  int pwkid=kid+1; //c++ to fortran index conversion.  pwscf index starts at 1.   
  if (nbands<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );

  eigenvals.resize(nbands);
  if(getNumSpins()==2)
    pwkid+=nktot*(spin);
  pwlib_getwfn_eigenvals_(&eigenvals[0], &pwkid); //&eigenvals[0] is the address of the first element in the vector array.
	                                       //spin not used.  Should be put in k-id index calculation.     
//    std::cout << "Energies:" << std::endl;
//  for (int i=0;i<8;i++)
//    std::cout << i << "  " << eigenvals[i] << std::endl;
  return;
}

void ESPWSCFInterface::getTwistData(std::vector<PosType>& TwistAngles,
		      std::vector<double>& TwistWeight, 
 		      std::vector<int>& TwistSymmetry)
{
   if(nktot<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
   
   double klist_tmp[nktot][3];
   TwistAngles.resize(nktot);
   TwistWeight.resize(nktot);
   TwistSymmetry.resize(nktot);
   std::cout << "calling getwfn_kpoints... I have " << nktot << " kpoints\n";
   pwlib_getwfn_kpoints_(&klist_tmp[0][0], &TwistWeight[0]);
   std::cout << "Done!\n";
   for (int i=0; i<nktot; i++)
   {
     for(int j=0; j<OHMMS_DIM; j++) TwistAngles[i][j]=-klist_tmp[i][j]; 
     TwistSymmetry[i]=0;
   }


  
   return;
}
int ESPWSCFInterface::getHasPsiG()
{
  return 1; //its pwscf.  of course it has psi_g.
}
    
void ESPWSCFInterface::getMeshSize(TinyVector<int,3> & mesh_out)
{
  pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
  mesh_out[0]=mesh[0];
  mesh_out[1]=mesh[1];
  mesh_out[2]=mesh[2];
  return;  
}
     
void ESPWSCFInterface::getReducedGVecs(std::vector<std::vector<TinyVector<int,3> > > & gvecs, int index)
{
  if (ngtot<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
  
  int gvec_tmp[ngtot][3];
  gvecs.resize(1);
  gvecs[0].resize(ngtot);

  pwlib_getwfn_gvecs_(gvec_tmp[0]);

  for (int i=0; i<ngtot; i++)
  {
    for(int j=0; j<OHMMS_DIM; j++) gvecs[0][i][j]=gvec_tmp[i][j];
  }
  /*std::cerr << "writing gvecs... \n";
  for (int i=0; i<ngtot; i++)
  {
    for(int j=0; j<OHMMS_DIM; j++) std::cerr << gvecs[0][i][j] << "\t" ;
    std::cerr << std::endl;
  }
  */
  return;
}

//bool ESPWSCFInterface::getPsi_kspace(Vector<complex<double> > & Cg, int orbid, int twistid)
//{
//        std::string s=psi_g_path(ti,spin,cur_bands[iorb].BandIndex); //Ray:  HWI (Handle with interface)
//        return h5f.read(cG,s); //RAY:  HWI  (Handle with interface).
//}

bool ESPWSCFInterface::getPsi_kspace(Vector<std::complex<double> > & cG,int spin, int orbid, int twistid)
{
  std::cout<<"Made it into ESPWSSCF Interface. "<<cG.size()<<" "<<spin<<" "<<orbid<<" "<<twistid<<std::endl;
  std::cout<<"ngtot="<<ngtot<<std::endl;
  int pwtwistid=twistid+1;
  int pworbid=orbid+1;
/////
  if(getNumSpins()==2)
    pwtwistid+=nktot*(spin);
/////
  if (ngtot<0) pwlib_getwfn_info_(&nbands, &nktot, &nkloc, mesh, &ngtot, &npw, &npwx );
  
  fftw_complex cgtmp[ngtot];
  cG.resize(ngtot);

  pwlib_getwfn_band_((double*) cgtmp[0], &pworbid, &pwtwistid);

  for(int i=0; i<ngtot; i++)  cG[i]=std::complex<double>(cgtmp[i][0], cgtmp[i][1]); //fftw_complex is a typedef of double x[2].  
  std::cout<<"The complex coefficient cG[0] is \t" << cG[0] << std::endl; 
  return true;   
}

void getParticleSet(ParticleSet& P){
  return;
};
}
