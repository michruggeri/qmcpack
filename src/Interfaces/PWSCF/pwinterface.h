
#ifndef PWINTERFACE_H
#define PWINTERFACE_H
//extern "C" void pwscf_();
//extern "C" void pwlib_init_(int * comm, char* inputFile,  char* outputFile, int * npool_lib,  int *p_id,  int * p_rk);

struct
{
  float r;
  float i;
} f90complex;

extern "C" void pwlib_init_(int * mpicomm);
extern "C" void pwlib_scf_();
//extern "C" void pwlib_getinfo(int * nat, int * charge, int * nup, int * ndown , int * nktot, 
//				int * nkloc_, double * box, double * kpts , double * wkp, 
//				int * ngrid , int * nbds, int * ngridx);
extern "C" void pwlib_getinfo_(int * nat, int * charge, int * nup, int * ndown, 
				int * nkloc, double * box, int * nktot,double * kpts) ;

extern "C" void pwlib_getbox_data_( double * box );
extern "C" void pwlib_getatom_info_( int * nat, int * nsp );
extern "C" void pwlib_getatom_data_( double * R, int * ityp);
extern "C" void pwlib_getspecies_data_( int * atomic_nums, int * valence_charge, double * masses, char * names );
extern "C" void pwlib_getelectron_info_( int * nelec, int * nup, int * ndown);
extern "C" void pwlib_getwfn_info_(int * nbands, int * nktot, int * nkloc, int * mesh, int * ngtot, int * npw, int * npwx );
extern "C" void pwlib_getwfn_kpoints_( double * klist, double * weights);
extern "C" void pwlib_getwfn_gvecs_(int * gvecs);
extern "C" void pwlib_getwfn_band_(double * cg, int * ibnd, int * ik);
extern "C" void pwlib_getwfn_eigenvals_(double * eig_, int * ik);
extern "C" void pwlib_getmpifft_(int * nnr, int * mesh);
#endif


