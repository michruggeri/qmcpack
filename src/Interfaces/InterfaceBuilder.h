
#ifndef QMCPLUSPLUS_INTERFACE_BUILDER_H
#define QMCPLUSPLUS_INTERFACE_BUILDER_H
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

class Communicate;

namespace qmcplusplus
{

/**
* @brief This class is a builder class for the generalized
*        interface to HDF5, external code, or other.   
*     
* @author Raymond Clay
*
* @ 
*/
  class InterfaceBase;

  class InterfaceBuilder
  {
  public:
    InterfaceBuilder(Communicate* mycomm): spointerface(0), ptclinterface(0), 
			                   myComm(mycomm), ptclNode(NULL), wfnNode(NULL) {};
    ~InterfaceBuilder(){};

    void put(xmlNodePtr cur);
    bool initialize();
   
    xmlNodePtr getPtclNode(){return ptclNode;};
    xmlNodePtr getWfnNode(){return wfnNode;};

    InterfaceBase* getSPOInterface(){return spointerface;};
    InterfaceBase* getPtclInterface(){return ptclinterface;};
      
  private:
     InterfaceBase* spointerface;
     InterfaceBase* ptclinterface;  

     Communicate* myComm;

     void putParticleSet(xmlNodePtr cur);
     void putSPOSet(xmlNodePtr cur);

     xmlNodePtr ptclNode;
     xmlNodePtr wfnNode;
  };
}

#endif
