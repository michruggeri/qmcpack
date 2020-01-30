
#ifndef QMCPLUSPLUS_ESPWSCFINTERFACE_BASE_H
#define QMCPLUSPLUS_ESPWSCFINTERFACE_H


#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/BandInfo.h"
#include "QMCWaveFunctions/AtomicOrbital.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Interfaces/ESInterfaceBase.h"
#include <map>

#define PW_NAME_LEN 3

class Communicate;

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
  class ESPWSCFInterface: public ESInterfaceBase
  {
  public:
    ESPWSCFInterface(Communicate* mycomm):ESInterfaceBase(mycomm), myComm(mycomm),PWFileName("pwscf.in"),
   					  nat(-1), nsp(-1), nbands(-1), nktot(-1), nkloc(-1),
					  nelec(-1), nup(-1), ndown(-1),npw(-1), npwx(-1), ngtot(-1),  box_changed(false),pos_changed(false) 
					 {mesh[0]=mesh[1]=mesh[2]=-1; ifacename="ESPWSCF";initialized=false;}; 
   ~ESPWSCFInterface(){};

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
    void initialize();
    void finalize(){};
    void update(){};
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
    bool put(xmlNodePtr cur);

    ///Interface IO methods
      
//    void getParticleSet(ParticleSet& P){}; ///<Overwrites
//    void setParticleSet(ParticleSet& P){};
      
//    void getSimulationCell(SimulationCell){};
//    void setSimulationCell(SimulationCell){};
     

    void getPrimVecs(Tensor<double,OHMMS_DIM>& primvec);
    void setPrimVecs(const Tensor<double,OHMMS_DIM>& primvec){};
   
    int getNumBands();
    int getNumSpins();
    int getNumTwists();
    int getNumCoreStates();
    int getNumMuffinTins();
    int getHaveDPsi();
    int getNumAtomicOrbitals();

    int getNumAtoms();
    void getNumElectrons(std::vector<int>& numspins);
 
    void setHaveDPsi(int have_dpsi){};
    void getVersion();

    int getNumSpecies();
    void setNumSpecies(int ns){};

    void getSpeciesIDs(ParticleIndex_t& sid);
    void setSpeciesIDs(const ParticleIndex_t& sid){};
    
    void getSpeciesData(Vector<int>& am, Vector<int>& charges, Vector<double>& masses, Vector<std::string>& names);
     
    void getAtomicNumbers(Vector<int> & am);
    void setAtomicNumbers(Vector<int> & am){}; 

    void getIonPositions(ParticlePos_t & ionpos);
    void setIonPositions(const ParticlePos_t & ionpos);

    //void printVersion();  
    void getTwistData(std::vector<PosType>& TwistAngles,
		      std::vector<double>& TwistWeight, 
 		      std::vector<int>& TwistSymmetry);
    
    void getAtomicOrbitals(std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals);
    void setAtomicOrbitals(const std::vector<AtomicOrbital<std::complex<double> > > & AtomicOrbitals){};

    void getOrbEigenvals(const int spin, const int kid, std::vector<double> & eigenvals);
    void setOrbEigenvals(const int spin, const int kid, const std::vector<double> & eigenvals){}; 
    
    int getHasPsiG();
    
    void getMeshSize(TinyVector<int,3> & mesh);
     
    void getReducedGVecs(std::vector<std::vector<TinyVector<int,3> > > & gvecs, int index=0);
    
    bool getPsi_kspace(Vector<std::complex<double> > & Cg, int spin, int orbid, int twistid);
    
  private:
    Communicate* myComm;
    std::string PWFileName;

    TinyVector<int,3> Version;  
    bool initialized;
    //Following variables are for persistence purposes;
    int nat;  //Number of atoms
    int nsp;  //number of species
    int nbands;  //number of bands
    int nktot; //total number of twists
    int nkloc; //number of twists on this node.  

    int nelec;
    int nup;
    int ndown;

    int npw; //number of plane wave coefficients 
    int npwx; //max allowed.  
    int ngtot; //number of nonzero pw coefficients < ecut

    int mesh[3];

    bool box_changed; ///< Flag for if box is modified. 
    bool pos_changed; ///< If the ion positions have changed.
    bool orb_changed; ///< if the SPO's change.  
    
    std::vector<int> species_ids;
    std::vector<double> species_masses; 
    std::vector<double> species_charges;
    std::vector<std::string> species_names;   
  /** return the path name in hdf5
   */
    inline std::string psi_g_path(int ti, int spin, int ib)
    {
      std::ostringstream path;
      path << "/electrons/kpoint_" << ti
           << "/spin_" << spin << "/state_" << ib << "/psi_g";
      return path.str();
    }

  /** return the path name in hdf5
   */
    inline std::string psi_r_path(int ti, int spin, int ib)
    {
      std::ostringstream path;
      path << "/electrons/kpoint_" << ti
           << "/spin_" << spin << "/state_" << ib << "/psi_r";
      return path.str();
    }
   
   //PW_NAME_LEN macro definition at top of this file.  
   inline Vector<std::string> convertCharToString(char namelist[][PW_NAME_LEN], int num_names)
   {

      Vector<std::string> namevec(num_names);
      for(int i=0; i<num_names; i++)
      {
        std::string tmpstr(namelist[i]);

        for (int j=tmpstr.size(); j>0; j--)
          if (tmpstr[j]==' ' || tmpstr[j]=='\n' || tmpstr[j]=='\t') tmpstr.erase(j);

        namevec[i]=tmpstr;
      }

      return namevec;
    }
  };
}

#endif
