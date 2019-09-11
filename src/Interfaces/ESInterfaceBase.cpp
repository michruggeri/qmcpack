#include "Interfaces/ESInterfaceBase.h"
#include "ParticleIO/ParticleIOUtility.h"
#include "ParticleBase/ParticleUtility.h"
//#include "ParticleIO/HDFParticleAttrib.h"
//#include "Numerics/HDFNumericAttrib.h"
//#include "OhmmsData/HDFStringAttrib.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/IteratorUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
namespace qmcplusplus
{
/*
void ESInterfaceBase::getParticleSet(ParticleSet& P)
{
  app_log()<<"ESInterfaceBase::getParticleSet()\n";
  //add basic attributes of the speciesset
  SpeciesSet& tspecies(P.getSpeciesSet());
  app_log()<<"Adding attributes...\n";
  int icharge= tspecies.addAttribute("charge");//charge_tag);
  int iatnumber= tspecies.addAttribute(atomic_number_tag);
  int membersize= tspecies.addAttribute("membersize");
  int massind= tspecies.addAttribute(mass_tag);

  //and now we initialize everything...
  Vector<int> atomicNumbers(0);
  Vector<string> speciesNames(0);
  Vector<double> speciesMass(0);
  Vector<int> speciesCharge(0);
  app_log()<<"getSpeciesData()\n";
  getSpeciesData(atomicNumbers,speciesCharge,speciesMass,speciesNames);
  app_log()<<"getNumSpecies()\n";
  int num_species=atomicNumbers.size();

  //set the species information
  for (int i=0; i< num_species; i++)
  {
    int ii=tspecies.addSpecies(speciesNames[i]);
    int q=-1;
    tspecies(icharge,ii)=speciesCharge[i];
    tspecies(iatnumber,ii)=atomicNumbers[i];
    tspecies(massind,ii)=speciesMass[i];
  }
  app_log()<<"grab Prim vecs\n";
  //get the cell dimensions
  Tensor<double,3> alat;
  getPrimVecs(alat);
  
  P.Lattice.set(alat);
  
  app_log()<<"get Species ID's\n";
  //get particle information
  getSpeciesIDs(P.GroupID);
  app_log()<<"ESInterfaceBase get natom\n"; 
  int natom=P.GroupID.size();
  P.create(natom);
  P.R.InUnit=PosUnit::CartesianUnit;
  getIonPositions(P.R);
  P.applyBC(P.R);
 
  P.PrimitiveLattice=P.Lattice;
  for(int i=0; i<P.getTotalNum(); ++i) P.ID[i]=i;
  P.PCID=P.ID;
  P.resetGroups();
  app_log()<<"done ....\n";
 
}

void ESInterfaceBase::bcastParticleSet(ParticleSet& P)
{
  if(myComm->size()==1) return;

  SpeciesSet& tspecies(P.getSpeciesSet());
  int nspecies=tspecies.getTotalNum();
  int natoms=P.getTotalNum();
  ostringstream o;
  if(myComm->rank()==0)
  {
    int i=0;
    for(; i<nspecies-1; ++i)
      o<<tspecies.speciesName[i]<<",";
    o<<tspecies.speciesName[i];
  }
  TinyVector<int,3> bsizes(nspecies,natoms,o.str().size()+1);
  myComm->bcast(bsizes);
  //send the names: UGLY!!!!
  nspecies=bsizes[0];
  char *species_names=new char[bsizes[2]];
  if(myComm->rank()==0)
    snprintf(species_names, bsizes[2], "%s",o.str().c_str());
  myComm->bcast(species_names,bsizes[2]);
  if(myComm->rank())
  {
    std::vector<string> vlist;
    parsewords(species_names,vlist);
    for(int i=0; i<vlist.size(); ++i)
      tspecies.addSpecies(vlist[i]);
    //create natoms particles
    P.create(bsizes[1]);
  }
  delete [] species_names;
  ParticleSet::Tensor_t lat(P.Lattice.R);
  ParticleSet::Buffer_t pbuffer;
  for(int i=0; i<tspecies.numAttributes(); ++i)
    pbuffer.add(tspecies.d_attrib[i]->begin(),tspecies.d_attrib[i]->end());
  pbuffer.add(lat.begin(),lat.end());
  pbuffer.add(get_first_address(P.R),get_last_address(P.R));
  pbuffer.add(P.GroupID.begin(),P.GroupID.end());
  myComm->bcast(pbuffer);
  P.R.InUnit=PosUnit::CartesianUnit;
  if(myComm->rank())
  {
    pbuffer.rewind();
    for(int i=0; i<tspecies.numAttributes(); ++i)
      pbuffer.get(tspecies.d_attrib[i]->begin(),tspecies.d_attrib[i]->end());
    pbuffer.get(lat.begin(),lat.end());
    pbuffer.get(get_first_address(P.R),get_last_address(P.R));
    pbuffer.get(P.GroupID.begin(),P.GroupID.end());
    P.Lattice.set(lat);
  }
}*/

void ESInterfaceBase::getParticleSet(ParticleSet & P)
{
  SpeciesSet& tspecies(P.getSpeciesSet());
  int icharge= tspecies.addAttribute("charge");//charge_tag);
  int iatnumber= tspecies.addAttribute(atomic_number_tag);
  int membersize= tspecies.addAttribute("membersize");
  int massind= tspecies.addAttribute(mass_tag);

  Vector<int> atomicNumbers(0);
  Vector<std::string> speciesNames(0);
  Vector<double> speciesMass(0);
  Vector<int> speciesCharge(0);
  getSpeciesData(atomicNumbers,speciesCharge,speciesMass,speciesNames);
  int nspecies=atomicNumbers.size();
  //add charge
  //atomic_number is optional
  for(int i=0; i<nspecies; ++i)
  {
    int ii=tspecies.addSpecies(speciesNames[i]);
    tspecies(icharge,ii)=speciesCharge[i];
    tspecies(iatnumber,ii)=atomicNumbers[i];
    //close the group
  }
  for(int ig=0; ig<nspecies; ++ig)
    tspecies(massind,ig)=1.0;
  //just for checking
  // tspecies(icharge,0)=15.0;
  // tspecies(icharge,1)=6.0;
  {
    //get the unit cell
    Tensor<double,3> alat;
    getPrimVecs(alat);
    P.Lattice.set(alat);
  }
  {
    //get the unit cell
    int natoms=getNumAtoms();
    
    P.create(natoms);
//    P.R.InUnit=PosUnit::CartesianUnit;
    P.R.InUnit=PosUnit::Cartesian;
    getIonPositions(P.R);
    P.applyBC(P.R); //force them [0,1)

    getSpeciesIDs(P.GroupID);
  }

  app_log()<<"atomicNumbers.size()="<<atomicNumbers.size()<<std::endl;
  app_log()<<"speciesNames.size()="<<speciesNames.size()<<std::endl;
  app_log()<<"speciesMass.size()="<<speciesMass.size()<<std::endl;
  app_log()<<"speciesCharge.size()="<<speciesCharge.size()<<std::endl;
  app_log()<<"P.groupID.size()="<<P.GroupID.size()<<std::endl;
  app_log()<<"P.R.size()="<<P.R.size()<<std::endl;
  app_log()<<"P.ID.size()="<<P.ID.size()<<std::endl;
  app_log()<<"P.getTotalNum()="<<P.getTotalNum()<<std::endl;
  P.PrimitiveLattice=P.Lattice;
  
  
  for(int i=0; i<P.getTotalNum(); ++i) P.ID[i]=i;
  P.PCID=P.ID;
  
  P.resetGroups();

}

void ESInterfaceBase::getElectronParticleSet(ParticleSet& P)
{

    std::vector<int> num_spin(0);
    getNumElectrons(num_spin);

    {
      //create species
      SpeciesSet& species=P.getSpeciesSet();
      //add up and down
      species.addSpecies("u");
      if(num_spin.size()>1)
        species.addSpecies("d");
      int chid=species.addAttribute("charge");
      for(int i=0; i<num_spin.size(); ++i)
        species(chid,i)=-1.0;
      int mid=species.addAttribute("membersize");
      for(int i=0; i<num_spin.size(); ++i)
        species(mid,i)=num_spin[i];
      mid=species.addAttribute("mass");
      for(int i=0; i<num_spin.size(); ++i)
        species(mid,i)=1.0;
      P.create(num_spin);
    }
    if(P.Lattice.SuperCellEnum == SUPERCELL_BULK)
    {
      makeUniformRandom(P.R);
//      P.R.setUnit(PosUnit::LatticeUnit);
      P.R.setUnit(PosUnit::Lattice);
      P.convert2Cart(P.R);
    }
    else
    {
      APP_ABORT("Error:  ESInterfaceBase::getElectronParticleSet only works for Bulk calcluations right now")
      //assign non-trivial positions for the quanmtum particles
    }

}
void ESInterfaceBase::bcastParticleSet(ParticleSet& P)
{
  myComm->barrier();
  std::cout<<"ESInterfaceBase::DEBUGGGG mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  SpeciesSet& tspecies(P.getSpeciesSet());
  if(myComm->size()==1)
    return; 
  int nspecies=tspecies.getTotalNum();
  int natoms=P.getTotalNum();
  int bnamesize=0;
  std::ostringstream o;
  if(myComm->rank()==0)
  {
    int i=0;
    for(; i<nspecies-1; ++i)
      o<<tspecies.speciesName[i]<<",";
    o<<tspecies.speciesName[i];
  }

  TinyVector<int,3> bsizes(nspecies,natoms,o.str().size()+1);
  myComm->bcast(bsizes);
  myComm->barrier();
  std::cout<<"ESInterfaceBase:: made it through barrier o= "<<o.str()<<" "<<bsizes[2]<<" "<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  
  //send the names: UGLY!!!!
  nspecies=bsizes[0];
  char *species_names=new char[bsizes[2]];
  if(myComm->rank()==0)
    snprintf(species_names, bsizes[2], "%s",o.str().c_str());

  std::cout<<"ESInterfaceBase:: before bcast species_names "<< species_names<<" mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  
  myComm->bcast(species_names,bsizes[2]);
  std::cout<<"ESInterfaceBase:: after bcast species_names  mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  myComm->barrier();
  if(myComm->rank()==0)
  {
    std::vector<std::string> vlist;
    parsewords(species_names,vlist);
    for(int i=0; i<vlist.size(); ++i)
      tspecies.addSpecies(vlist[i]);
    //create natoms particles
 //   P.create(bsizes[1]);
  }
  delete [] species_names;
  ParticleSet::Tensor_t lat(P.Lattice.R);
  ParticleSet::Buffer_t pbuffer;
  
  for(int i=0; i<tspecies.numAttributes(); ++i)
    pbuffer.add(tspecies.d_attrib[i]->begin(),tspecies.d_attrib[i]->end());
  pbuffer.add(lat.begin(),lat.end());
  std::cout<<"ESInterfaceBase:: lat len = "<<lat.end()-lat.begin()<<" mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  pbuffer.add(get_first_address(P.R),get_last_address(P.R));
  std::cout<<"ESInterfaceBase:: p.r len = "<<get_last_address(P.R)-get_first_address(P.R)<<" mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  pbuffer.add(P.GroupID.begin(),P.GroupID.end());
  std::cout<<"ESInterfaceBase:: pgroupid len = "<<P.GroupID.end()-P.GroupID.begin()<<" mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;
  std::cout<<"ESInterfaceBase::Before bcast pbuffer. pbuffer.size="<<pbuffer.size()<<"  mycomm="<<myComm<<" size="<<myComm->size()<<" rank="<<myComm->rank()<<std::endl;

  myComm->bcast(pbuffer);
  myComm->barrier();
//  P.R.InUnit=PosUnit::CartesianUnit;
  P.R.InUnit=PosUnit::Cartesian;
  if(myComm->rank())
  {
    pbuffer.rewind();
    for(int i=0; i<tspecies.numAttributes(); ++i)
      pbuffer.get(tspecies.d_attrib[i]->begin(),tspecies.d_attrib[i]->end());
    pbuffer.get(lat.begin(),lat.end());
    pbuffer.get(get_first_address(P.R),get_last_address(P.R));
    pbuffer.get(P.GroupID.begin(),P.GroupID.end());
    P.Lattice.set(lat);
  }
  std::cout<<"end barrier\n";
  myComm->barrier();
  return ;

}

}
