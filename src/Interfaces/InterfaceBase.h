
#ifndef QMCPLUSPLUS_INTERFACE_BASE_H
#define QMCPLUSPLUS_INTERFACE_BASE_H

#include <Particle/ParticleSet.h>
#include "Configuration.h"
#include "Message/MPIObjectBase.h"
namespace qmcplusplus
{

  //class MPIObjectBase;
  //class ParticleSet;
/**
* @brief This class is a generalization/abstraction of the 
*        file IO that occurs when setting up the particlesets, 
*        lattices, and orbitals. Provides a uniform class for
*        hdf5 and external codes. 
*     
* @author Raymond Clay
*
* @ 
*/
  class InterfaceBase: public QMCTraits
  //class InterfaceBase: public MPIObjectBase, public QMCTraits
  {
    ///Derives from MPIObjectBase so that we can get access to the Communicate object
    ///and all the good mpi functions.  Just in case.  
    ///
    ///Derives from QMCTraits for things like double, int, PosType, etc.
    ///  See src/Configuration.h 
    ///
    /// Henceforth, don't use std::vector or vector.  Use the PETE Vector, found in 
    ///   src/OhmmsPETE/OhmmsVector.h
   ///
  public:
    //typedef ParticleSet::ParticleLayout_t Lattice_t ///< typedef to ParticleLayout_t found in src/Particle/ParticleBase.h

    InterfaceBase(Communicate* mpicomm):initialized("false"),myComm(mpicomm),ifacename("base") {};
    InterfaceBase():initialized("false"), ifacename("base") {};
    ~InterfaceBase(){};

    ///Interface startup and shutdown methods///
    ///
      
    /**
    * @brief Initialize runs the necessary initialization steps before
    *        the interface can be used.  This can be checking and opening
    *        an HDF5 file, or it could be starting up an external code.  
    *
    */
    virtual void initialize()=0;
    virtual void finalize()=0;
    
    /** 
    *@brief Updates the current state of the interface.  Presumably, changes could be made to
    *       The particle positions, the box, etc., which the interface would then respond to
    *       to give proper wavefunctions or other info.  
    */
    virtual void update()=0;
    /**
    * @brief put() is called before initialization.  It is responsible for
    *        reading all relevant interface state info prior to initializaiton.
    *        This could be things like the hdf5 file path, or an external code's
    *        input and output files.   
    *
    * @param xml Pointer to the relevant xml node.  
    *
    * @return True/False if the initialization is successful.  
    */
    virtual bool put(xmlNodePtr cur)=0;

    std::string getInterfaceName(){return ifacename;};
    ///Interface IO methods
    
    //Ion Particle Set utilities  
    virtual void getParticleSet(ParticleSet& P)=0;
    virtual void setParticleSet(ParticleSet& P){};
    
    //Electron Particle Set Utilities
    virtual void getElectronParticleSet(ParticleSet& P)=0;
    virtual void setElectronParticleSet(ParticleSet& P){};

    virtual void bcastParticleSet(ParticleSet& P)=0;  
//    virtual void getSimulationCell(SimulationCell)=0;
//    virtual void setSimulationCell(SimulationCell);

   protected:
    bool initialized;
    Communicate* myComm;
//    bool box_changed; ///< Flag for if box is modified. 
    bool pos_changed; ///< If the ion positions have changed.
    bool orb_changed; ///< if the SPO's change.  
    
    std::string ifacename; 
     
  };
}

#endif
