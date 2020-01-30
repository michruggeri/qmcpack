
#ifndef QMCPLUSPLUS_ESINTERFACE_BASE_H
#define QMCPLUSPLUS_ESINTERFACE_BASE_H

#include "Interfaces/InterfaceBase.h"
#include <vector>
#include <complex>
#include <string>
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "Particle/ParticleSet.h"
namespace qmcplusplus
{

/**
* @brief This class derives from InterfaceBase, and is specialized
*        for the IO required for the EinsplineSet builder.   
*     
* @author Raymond Clay
*
* @ 
*/
  class ESInterfaceBase: public InterfaceBase
  {
  public:
    ESInterfaceBase(Communicate* mycomm) : InterfaceBase(mycomm),myComm(mycomm),atomic_number_tag("atomic_number"),mass_tag("mass") {ifacename="ES";};
    ESInterfaceBase() :atomic_number_tag("atomic_number"),mass_tag("mass") {ifacename="ES";};
    ~ESInterfaceBase(){};

    typedef typename ParticleSet::ParticleIndex_t ParticleIndex_t;
    typedef typename ParticleSet::ParticlePos_t   ParticlePos_t;   
    ///Interface startup and shutdown methods///
    ///
      
    /**
    * @brief Initialize runs the necessary initialization steps before
    *        the interface can be used.  This can be checking and opening
    *        an HDF5 file, or it could be starting up an external code.  
    *
    * @param mpicomm In case the interface has a nontrivial/distributed
    *        structure (as an external code might), hand it that appropraite
    *        mpi communicator.   
    */
    virtual void initialize()=0;
    virtual void finalize()=0;
    virtual void update()=0;
    /**
    * @brief Put() is called before initialization.  It is responsible for
    *        reading all relevant interface state info prior to initializaiton.
    *        This could be things like the hdf5 file path, or an external code's
    *        input and output files.   
    *
    * @param xml Pointer to the relevant xml node.  
    *
    * @return True/False if the initialization is successful.  
    */
    virtual bool put(xmlNodePtr cur)=0;

    ///Interface IO methods
      
    virtual void getParticleSet(ParticleSet& P); ///<Overwrites
    virtual void setParticleSet(ParticleSet& P){};
    
    virtual void getElectronParticleSet(ParticleSet& P);
   
    virtual void bcastParticleSet(ParticleSet& P);  
//    void getSimulationCell(SimulationCell);
//    void setSimulationCell(SimulationCell);
     
    ///All following get/set methods are specific to EinsplineSet objects.

    ///CELL GEOMETRY  
    virtual void getPrimVecs(Tensor<double,OHMMS_DIM>& primvec)=0;
    virtual void setPrimVecs(const Tensor<double,OHMMS_DIM>& primvec)=0;

    virtual int getNumTwists()=0;
    virtual void setNumTwists(){};

    virtual void getReducedTwistAngles(){};
    virtual void setReducedTwistAngles(){};

    virtual void getTwistSymmetry(){};
    virtual void setTwistSymmetry(){};

    virtual void getTwistWeight(){};
    virtual void setTwistWeight(){};

    virtual void getTwistData(std::vector<PosType>& TwistAngles,
		      std::vector<double>& TwistWeight, 
 		      std::vector<int>& TwistSymmetry)=0;
    
    ///ION INFO
    virtual void getSpeciesIDs(ParticleIndex_t& )=0;
    virtual void setSpeciesIDs(const ParticleIndex_t& ){};
    
    virtual void getSpeciesData(Vector<int>& am, Vector<int>& charges, Vector<double>& masses, Vector<std::string>& names)=0;
   
    virtual int getNumSpecies()=0;
    virtual void setNumSpecies(const int num){};
   
    virtual void getAtomicNumbers(Vector<int> & am)=0;
    virtual void setAtomicNumbers(Vector<int> & am){}; 
   
    virtual  int getNumAtoms()=0;
    virtual void setNumAtoms(const int num){};

    virtual void getIonPositions( ParticlePos_t & ionpos)=0;
    virtual void setIonPositions(const ParticlePos_t & ionpos){};
    
    ///Electron Info
    virtual void getNumElectrons(std::vector<int>& num_spins)=0;
    virtual void setNumElectrons(const std::vector<int>& num_spins){};

     ///ORBITAL INFO
    virtual int getNumBands()=0;
    virtual void setNumBands(const int num){};

    virtual int getNumCoreStates()=0;
    virtual void setNumCoreStates(const int){};

    virtual int getNumSpins()=0;
    virtual void setNumSpins(const int){};

    virtual int getNumMuffinTins()=0;
    virtual void setNumMuffinTins(const int){};

    virtual int getNumAtomicOrbitals()=0;
    virtual void setNumAtomicOrbitals(const int){};
    
    virtual void getOrbEigenvals(const int spin, const int kid, std::vector<double> & eigenvals)=0;
    virtual void setOrbEigenvals(const int spin, const int kid, const std::vector<double> & eigenvals){};
    ///ORBITAL DATA
    virtual int getHasPsiG()=0;
    
    virtual void getMeshSize(TinyVector<int,3> & mesh)=0;
     
    virtual void getReducedGVecs(std::vector<std::vector<TinyVector<int,3> > > & gvecs, int index=0)=0;
    
         
    virtual bool getPsi_kspace(Vector<std::complex<double> > & Cg, int spin, int orbid, int twistid)=0;
//    virtual bool getPsi_rspace(data, int orbid, int twistid){return true;};
    
     
    //DENSITY INFO
 //   void getDensityReducedGvecs(){};
 //   void setDensityReducedGvecs(const int){};

    //DENSITY 
 
    //ORBITAL DATA
    virtual void getAtomicOrbitals(std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals)=0;
    virtual void setAtomicOrbitals(const std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals){};
   
    virtual int getHaveDPsi(){return 0;};
    virtual void setHaveDPsi(int have_dpsi){}; 
    virtual void getVersion(){};  
      
  private:
///From InterfaceBase:
//    bool box_changed;  
//    bool pos_changed; 
//    bool orb_changed;   
    Communicate* myComm;  
    bool orb_in_kspace; ///< If orbital coefficients are in k-space, true
    bool orb_in_rspace; ///< If orbital coefficients are in real-space, true 

    bool density_in_kspace; ///< If density info is in k-space, true
    bool density_in_rspace; ///< if density info is in r-space, true

    std::string atomic_number_tag;
    std::string mass_tag;

  };
}

#endif
