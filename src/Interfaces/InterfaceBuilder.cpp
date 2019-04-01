#include "Interfaces/InterfaceBuilder.h"
#include "Interfaces/InterfaceBase.h"
#include "Message/CommOperators.h"
#include "qmc_common.h"
#include "Configuration.h"
#include "OhmmsData/AttributeSet.h"

#include "ESHDF5/ESHDF5Interface.h"
#include "PWSCF/ESPWSCFInterface.h"

namespace qmcplusplus
{


//Here's an example interface tag.  
//  <interface>
//     <sposet name=" " type=" " href=" "/>
//     <particleset  name=" " type=" " href=" "/>
//  </interface>
void InterfaceBuilder::put(xmlNodePtr cur)
{


  xmlNodePtr tcur=cur->children;
  while(tcur != NULL)
  { //check <determinantset/> or <sposet_builder/> to extract the ionic and electronic structure
    std::string cname((const char*)tcur->name);
    if(cname == "sposet")
    {
        app_log()<<"  found sposet tag in interface.  put\n";
	putSPOSet(tcur);
    }
    else if (cname == "particleset")
    {
        app_log()<<"  found particleset tag in interface.  put\n";
   	putParticleSet(tcur);
    }
    tcur=tcur->next;
  }
}

void InterfaceBuilder::putParticleSet(xmlNodePtr cur)
{
  ptclNode=cur;
  OhmmsAttributeSet a;
  std::string ifacetype;
  std::string sourcepath;
  std::string ptclsetname;

  a.add(ifacetype,"type");
  a.add(sourcepath,"href");
  a.add(ptclsetname,"name");
  a.put(cur);
   
  app_log()<<" putParticleSet(..).  Found the following parameters...\n";
  app_log()<<"  ifacetype="<<ifacetype<<" sourcepath="<<sourcepath<<std::endl;
 
  if(ifacetype=="ESHDF")
  {
    if (ptclinterface !=0  && ptclinterface->getInterfaceName() != "ESHDF")
    {
      app_log()<<"ESHDF. ptclinterface!=0 and interfacename !=ESHDF\n";  
      delete ptclinterface;
      ptclinterface=new ESHDF5Interface(myComm);
    }
    else
    {
      ptclinterface=new ESHDF5Interface(myComm);
    }
  }
/*  else if(ifacetype=="ESPWSCF")
  {
    app_log() << "In ESPWSCF...\n";
    app_log()<< "ptclinterface="<<ptclinterface<<std::endl;
    if (ptclinterface !=0)
    {
       if (ptclinterface->getInterfaceName() != "ESPWSCF")
       {
          app_log()<<"ESPWSCF. ptclinterface!=0 and interfacename !=ESPWSCF\n";  
          delete ptclinterface;
          ptclinterface=new ESPWSCFInterface(myComm);
       }
    }
    else
    {
      app_log()<<"Create ESPWSCF from scratch...\n";
      ptclinterface=new ESPWSCFInterface(myComm);
    }
  }*/
  else
  {
	APP_ABORT("ERROR:  Interface type for ParticleSet not valid");
  }

  ptclinterface->put(cur);

}

void InterfaceBuilder::putSPOSet(xmlNodePtr cur)
{
  wfnNode=cur;
 
  OhmmsAttributeSet a;
  std::string ifacetype("");
  std::string sourcepath("");
  std::string sposetname("myspo");
  std::string fromscratch("no");

  a.add(ifacetype,"type");
  a.add(sourcepath,"href");
  a.add(sposetname,"name");
  a.add(fromscratch,"fromscratch");
  a.put(cur);

  if(sourcepath=="") APP_ABORT("InterfaceBuilder ERROR: No path for inputfile specified for interface");
 
  //First, check to see if the SPOset interface is already instantiated in the ptcl set builder.  
  //
  if ( ptclinterface !=0 )
  {
    if (ptclinterface->getInterfaceName() == ifacetype)
    {
      spointerface=ptclinterface;
      return;
    }
  }

  //Now, if sposet isn't already initialized in ptclinterface, then go through
  //and initialize it if it doesn't already exist.   
  if(ifacetype=="ESHDF")
  {
    if (spointerface !=0  && spointerface->getInterfaceName() != "ESHDF")
    {
      delete spointerface;
      spointerface=new ESHDF5Interface(myComm);
    }
    else
    {
      spointerface=new ESHDF5Interface(myComm);
    }
  }
/*  else if(ifacetype=="ESPWSCF")
  {
    if (spointerface !=0  && spointerface->getInterfaceName() != "ESPWSCF")
    {
      delete spointerface;
      spointerface=new ESPWSCFInterface(myComm);
    }
    else
    {
      spointerface=new ESPWSCFInterface(myComm);
    }
  }
*/  else
  {
	APP_ABORT("ERROR:  Interface type for SPOset not valid");
  }
  spointerface->put(cur);
}

bool InterfaceBuilder::initialize()
{
   std::cerr << "Beginning the initialization... \n";
  app_log()<<"InterfaceBuilder::initialize() started\n";
  if(ptclinterface!=0)
  {
    ptclinterface->initialize();    
  }
 
  if(spointerface!=0)
  {
    spointerface->initialize();
  }
  app_log()<<"InterfaceBuilder::initialize() completed\n";
}
  

}
