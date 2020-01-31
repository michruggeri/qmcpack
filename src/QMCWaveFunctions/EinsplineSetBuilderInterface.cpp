//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
#include "qmc_common.h"

//#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <vector>
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/EinsplineSetBuilderInterface.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderInterface.h"
//#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorBase.h"

//#include "qmc_common.h"

namespace qmcplusplus
{

// Constructor (from EinsplineSetBuilderCommon.cpp)

EinsplineSetBuilderInterface::EinsplineSetBuilderInterface(ParticleSet& p, PtclPoolType& psets, Communicate *comm, xmlNodePtr cur,std::string label)
: SPOSetBuilder(comm), TargetPtcl(p),ParticleSets(psets), MixedSplineReader(0),esinterface(0), XMLRoot(cur),
intlabel(label),H5FileID(-1),Format(QMCPACK),
NumBands(0), NumElectrons(0), NumSpins(0), NumTwists(0),
NumCoreStates(0),
 MeshFactor(1.0), MeshSize(0,0,0),TileFactor(1,1,1),TwistNum(0),LastSpinSet(-1),NumOrbitalsRead(-1),NumMuffinTins(0),makeRotations(false)
{
  //assume one, not safe!! 
  myTableIndex=1;

  MatchingTol=10*std::numeric_limits<float>::epsilon();
  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
  TileMatrix(i,j) = 0;
  MyToken=0;

  //invalidate states by the basis class
  delete_iter(states.begin(),states.end());
  states.clear();
  states.resize(p.groups(),0);

  //create vectors with nullptr
  FullBands.resize(p.groups(),0);

  //HACK HACK HACK
  //For interface debugging.
  //Proper construction of interface should be done in some put(xml) command,
  //and input file should be parsed.
//  ESHDF5Interface* myint;  // To be replaced with the proper type of interface (e.g. pwscf)
//  myint = new ESHDF5Interface(myComm);
  if(intlabel=="PWSCF"){
#if !defined(QE_INTERFACE)
      app_error() << "Cannot use the Espresso interface \n"
                  << "Please recompile with QE_INTERFACE=1.\n";
      app_error().flush();
      APP_ABORT("EinsplineSetBuilderInterface::EinsplineSetBuilderInterface");
#else
    ESPWSCFInterface* myint;  // To be replaced with the proper type of interface (e.g. pwscf)
    myint = new ESPWSCFInterface(myComm);
    myint->put(cur);
    esinterface=static_cast<ESInterfaceBase*>(myint);
#endif
  }
  if(intlabel=="ESHDF"){
    ESHDF5Interface* myint;  // To be replaced with the proper type of interface (e.g. pwscf)
    myint = new ESHDF5Interface(myComm);
    myint->put(cur);
    esinterface=static_cast<ESInterfaceBase*>(myint);
  }
//  esinterface->initialize();
}

template<typename T>
inline TinyVector<T,3>
IntPart (const TinyVector<T,3>& twist)
{
  return TinyVector<T,3> (round(twist[0]-1.0e-6),
		  round(twist[1]-1.0e-6),
		  round(twist[2]-1.0e-6));
}

template<typename T>
inline TinyVector<T,3>
FracPart (const TinyVector<T,3>& twist)
{
  return twist - IntPart (twist);
}
// Destructor (from EinsplineSetBuilderCommon.cpp)

EinsplineSetBuilderInterface::~EinsplineSetBuilderInterface()
{
  DEBUG_MEMORY("EinsplineSetBuilderInterface::~EinsplineSetBuilderInterface");
  if(MixedSplineReader) delete MixedSplineReader;
  if(H5FileID>=0)
  H5Fclose(H5FileID);
}

void EinsplineSetBuilderInterface::BroadcastOrbitalInfo()
{
  if(myComm->size() == 1)
    return;
  int numIons = IonTypes.size();
  int numAtomicOrbitals = AtomicOrbitals.size();
  int numDensityGvecs = TargetPtcl.DensityReducedGvecs.size();
  PooledData<double> abuffer;
  PooledData<int>       aibuffer;
  aibuffer.add(Version.begin(),Version.end()); //myComm->bcast(Version);
  aibuffer.add(Format);
  abuffer.add(Lattice.begin(),Lattice.end());//myComm->bcast(Lattice);
  abuffer.add(RecipLattice.begin(),RecipLattice.end()); //myComm->bcast(RecipLattice);
  abuffer.add(SuperLattice.begin(),SuperLattice.end()); //myComm->bcast(SuperLattice);
  abuffer.add(LatticeInv.begin(),LatticeInv.end()); //myComm->bcast(LatticeInv);
  aibuffer.add(NumBands); //myComm->bcast(NumBands);
  aibuffer.add(NumElectrons); //myComm->bcast(NumElectrons);
  aibuffer.add(NumSpins); //myComm->bcast(NumSpins);
  aibuffer.add(NumTwists); //myComm->bcast(NumTwists);
  aibuffer.add(numIons); //myComm->bcast(numIons);
  aibuffer.add(NumMuffinTins);
  aibuffer.add(numAtomicOrbitals);
  aibuffer.add(numDensityGvecs);
  aibuffer.add(HaveOrbDerivs); 
  myComm->bcast(abuffer);
  myComm->bcast(aibuffer);
  if(myComm->rank())
  {
    abuffer.rewind();
    aibuffer.rewind();
    aibuffer.get(Version.begin(),Version.end());
    aibuffer.get(Format);
    abuffer.get(Lattice.begin(),Lattice.end());
    abuffer.get(RecipLattice.begin(),RecipLattice.end());
    abuffer.get(SuperLattice.begin(),SuperLattice.end());
    abuffer.get(LatticeInv.begin(),LatticeInv.end());
    aibuffer.get(NumBands);
    aibuffer.get(NumElectrons);
    aibuffer.get(NumSpins);
    aibuffer.get(NumTwists);
    aibuffer.get(numIons);
    aibuffer.get(NumMuffinTins);
    aibuffer.get(numAtomicOrbitals);
    aibuffer.get(numDensityGvecs);
    aibuffer.get(HaveOrbDerivs);
    MT_APW_radii.resize(NumMuffinTins);
    MT_APW_lmax.resize(NumMuffinTins);
    MT_APW_rgrids.resize(NumMuffinTins);
    MT_APW_num_radial_points.resize(NumMuffinTins);
    MT_centers.resize(NumMuffinTins);
    TargetPtcl.DensityReducedGvecs.resize(numDensityGvecs);
    TargetPtcl.Density_G.resize(numDensityGvecs);
    AtomicOrbitals.resize(numAtomicOrbitals);
  }
  std::vector<int> rgrids_sizes(NumMuffinTins);
  for (int tin=0; tin<NumMuffinTins; tin++)
    rgrids_sizes[tin] = MT_APW_rgrids[tin].size();
  myComm->bcast(rgrids_sizes);
  if (myComm->rank())
    for (int tin=0; tin<NumMuffinTins; tin++)
      MT_APW_rgrids[tin].resize(rgrids_sizes[tin]);
  if (IonTypes.size() != numIons)
  {
    IonTypes.resize(numIons);
    IonPos.resize(numIons);
  }
  //new buffer
  PooledData<double> bbuffer;
  PooledData<int> bibuffer;
  for(int i=0; i<numIons; ++i)
    bibuffer.add(IonTypes[i]);
  //myComm->bcast(IonTypes);
  bbuffer.add(&IonPos[0][0],&IonPos[0][0]+OHMMS_DIM*numIons);
  //myComm->bcast(IonPos);
  if (TwistAngles.size() != NumTwists)
    TwistAngles.resize(NumTwists);
  bbuffer.add(&TwistAngles[0][0],&TwistAngles[0][0]+OHMMS_DIM*NumTwists);
  //myComm->bcast(TwistAngles);
  if (TwistSymmetry.size() != NumTwists)
    TwistSymmetry.resize(NumTwists);
  bibuffer.add(&TwistSymmetry[0],&TwistSymmetry[0]+NumTwists);
  if (TwistWeight.size() != NumTwists)
    TwistWeight.resize(NumTwists);
  bibuffer.add(&TwistWeight[0],&TwistWeight[0]+NumTwists);
  bbuffer.add(MT_APW_radii.begin(), MT_APW_radii.end());
  bibuffer.add(MT_APW_lmax.begin(),  MT_APW_lmax.end());
  bibuffer.add(MT_APW_num_radial_points.begin(),
       MT_APW_num_radial_points.end());
  bbuffer.add(&(MT_centers[0][0]), &(MT_centers[0][0])+OHMMS_DIM*NumMuffinTins);
  for (int i=0; i<NumMuffinTins; i++)
    bbuffer.add(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
  bibuffer.add(&(TargetPtcl.DensityReducedGvecs[0][0]),
       &(TargetPtcl.DensityReducedGvecs[0][0])+numDensityGvecs*OHMMS_DIM);
  bbuffer.add(&(TargetPtcl.Density_G[0]),
      &(TargetPtcl.Density_G[0]) + numDensityGvecs);
  for (int iat=0; iat<numAtomicOrbitals; iat++)
  {
    AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
    bibuffer.add (orb.SplinePoints);
    bibuffer.add (orb.PolyOrder);
    bibuffer.add (orb.lMax);
    bibuffer.add (orb.Numlm);
    bbuffer.add  (&orb.Pos[0], &orb.Pos[0]+OHMMS_DIM);
    bbuffer.add  (orb.CutoffRadius);
    bbuffer.add  (orb.SplineRadius);
    bbuffer.add  (orb.PolyRadius);
  }
  myComm->bcast(bbuffer);
  myComm->bcast(bibuffer);
  if(myComm->rank())
  {
    bbuffer.rewind();
    bibuffer.rewind();
    for(int i=0; i<numIons; ++i)
      bibuffer.get(IonTypes[i]);
    bbuffer.get(&IonPos[0][0],&IonPos[0][0]+OHMMS_DIM*numIons);
    bbuffer.get(&TwistAngles[0][0],&TwistAngles[0][0]+OHMMS_DIM*NumTwists);
    bibuffer.get(&TwistSymmetry[0],&TwistSymmetry[0]+NumTwists);
    bibuffer.get(&TwistWeight[0],&TwistWeight[0]+NumTwists);
    bbuffer.get(MT_APW_radii.begin(), MT_APW_radii.end());
    bibuffer.get(MT_APW_lmax.begin(),  MT_APW_lmax.end());
    bibuffer.get(MT_APW_num_radial_points.begin(),
	 MT_APW_num_radial_points.end());
    bbuffer.get(&(MT_centers[0][0]),
	&(MT_centers[0][0])+OHMMS_DIM*NumMuffinTins);
    for (int i=0; i<NumMuffinTins; i++)
      bbuffer.get(MT_APW_rgrids[i].begin(), MT_APW_rgrids[i].end());
    bibuffer.get(&(TargetPtcl.DensityReducedGvecs[0][0]),
	 &(TargetPtcl.DensityReducedGvecs[0][0])+
	 numDensityGvecs*OHMMS_DIM);
    bbuffer.get(&(TargetPtcl.Density_G[0]),
	&(TargetPtcl.Density_G[0]) + numDensityGvecs);
    for (int iat=0; iat<numAtomicOrbitals; iat++)
    {
      AtomicOrbital<std::complex<double> > &orb = AtomicOrbitals[iat];
      bibuffer.get (orb.SplinePoints);
      bibuffer.get (orb.PolyOrder);
      bibuffer.get (orb.lMax);
      bibuffer.get (orb.Numlm);
      bbuffer.get  (&orb.Pos[0], &orb.Pos[0]+OHMMS_DIM);
      bbuffer.get  (orb.CutoffRadius);
      bbuffer.get  (orb.SplineRadius);
      bbuffer.get  (orb.PolyRadius);
    }
  }
//buffer to bcast hybrid representation atomic orbital info
  PooledData<double> cbuffer;
  PooledData<int> cibuffer;
  myComm->bcast(cbuffer);
  myComm->bcast(cibuffer);
  AtomicCentersInfo.resize(numIons);
  Super2Prim.resize(SourcePtcl->R.size());
  cbuffer.add(AtomicCentersInfo.inner_cutoff.begin(), AtomicCentersInfo.inner_cutoff.end());
  cbuffer.add(AtomicCentersInfo.non_overlapping_radius.begin(), AtomicCentersInfo.non_overlapping_radius.end());
  cbuffer.add(AtomicCentersInfo.cutoff.begin(), AtomicCentersInfo.cutoff.end());
  cbuffer.add(AtomicCentersInfo.spline_radius.begin(), AtomicCentersInfo.spline_radius.end());
  cibuffer.add(Super2Prim.begin(),Super2Prim.end());
  cibuffer.add(AtomicCentersInfo.lmax.begin(), AtomicCentersInfo.lmax.end());
  cibuffer.add(AtomicCentersInfo.GroupID.begin(), AtomicCentersInfo.GroupID.end());
  cibuffer.add(AtomicCentersInfo.spline_npoints.begin(), AtomicCentersInfo.spline_npoints.end());
  myComm->bcast(cbuffer);
  myComm->bcast(cibuffer);
  if(myComm->rank())
  {
    cbuffer.rewind();
    cibuffer.rewind();
    cbuffer.get(AtomicCentersInfo.inner_cutoff.begin(), AtomicCentersInfo.inner_cutoff.end());
    cbuffer.get(AtomicCentersInfo.non_overlapping_radius.begin(), AtomicCentersInfo.non_overlapping_radius.end());
    cbuffer.get(AtomicCentersInfo.cutoff.begin(), AtomicCentersInfo.cutoff.end());
    cbuffer.get(AtomicCentersInfo.spline_radius.begin(), AtomicCentersInfo.spline_radius.end());
    cibuffer.get(Super2Prim.begin(),Super2Prim.end());
    cibuffer.get(AtomicCentersInfo.lmax.begin(), AtomicCentersInfo.lmax.end());
    cibuffer.get(AtomicCentersInfo.GroupID.begin(), AtomicCentersInfo.GroupID.end());
    cibuffer.get(AtomicCentersInfo.spline_npoints.begin(), AtomicCentersInfo.spline_npoints.end());
    for (int i=0; i<numIons; i++)
      AtomicCentersInfo.ion_pos[i]=IonPos[i];
  }
  return;
}

bool EinsplineSetBuilderInterface::CheckLattice()
{
  update_token(__FILE__,__LINE__,"CheckLattice");

  double diff=0.0;
  for (int i=0; i<OHMMS_DIM; i++)
  for (int j=0; j<OHMMS_DIM; j++)
  {
    double max_abs=std::max(std::abs(SuperLattice(i,j)),static_cast<double>(std::abs(TargetPtcl.Lattice.R(i,j))));
    if(max_abs>MatchingTol)
      diff=std::max(diff,std::abs(SuperLattice(i,j) - TargetPtcl.Lattice.R(i,j))/max_abs);
  }
  if(diff>MatchingTol)
  {
    std::ostringstream o;
    o.setf(std::ios::scientific, std::ios::floatfield);
    o.precision(6);
    o << "EinsplineSetBuilder::ReadOrbitalInfo_ESHDF \n"
      << "Mismatched supercell lattices.\n";
    o << " Lattice in ESHDF5 " << std::endl;
    o << SuperLattice << std::endl;
    o << " Lattice in xml" << std::endl;
    o << TargetPtcl.Lattice.R << std::endl;
    o << " Difference " << std::endl;
    o << SuperLattice-TargetPtcl.Lattice.R << std::endl;
    o << " Max relative error = "<< diff << std::endl;
    o << " Tolerance      = "<< MatchingTol << std::endl;
    app_error() << o.str();
    return false;
  }
  return true;
}

// Should I get the other stuff here too?

void EinsplineSetBuilderInterface::set_metadata(int numOrbs, int TwistNum_inp)
{
// 1. set a lot of internal parameters in the EinsplineSetBuilder class
//  e.g. TileMatrix, UseRealOrbitals, DistinctTwists, MakeTwoCopies.
// 2. this is also where metadata for the orbitals are read from the wavefunction hdf5 file
//  and broadcast to MPI groups. Variables broadcasted are listed in 
//  EinsplineSetBuilderCommon.cpp EinsplineSetBuilder::BroadcastOrbitalInfo()
//   
//  if(myComm->rank()==0)
{

  Timer orb_info_timer;
// The tiling can be set by a simple vector, (e.g. 2x2x2), or by a
// full 3x3 matrix of integers.  If the tilematrix was not set in
// the input file...
  bool matrixNotSet = true;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      matrixNotSet = matrixNotSet && (TileMatrix(i,j) == 0);
// then set the matrix to what may have been specified in the
// tiling vector
  if (matrixNotSet)
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        TileMatrix(i,j) = (i==j) ? TileFactor[i] : 0;
  char buff[1000];
  if (myComm->rank() == 0)
  {
    snprintf (buff, 1000, "  TileMatrix = \n [ %2d %2d %2d\n   %2d %2d %2d\n   %2d %2d %2d ]\n",
       TileMatrix(0,0), TileMatrix(0,1), TileMatrix(0,2),
       TileMatrix(1,0), TileMatrix(1,1), TileMatrix(1,2),
       TileMatrix(2,0), TileMatrix(2,1), TileMatrix(2,2));
    app_log() << buff;
  }  
  if (numOrbs == 0)
  {
    app_error() << "You must specify the number of orbitals in the input file.\n";
    APP_ABORT("EinsplineSetBuilder::createSPOSet");
  }
  else
    app_log() << "  Generating " << numOrbs << " orbitals from Quantum Espresso file.\n";
  orb_info_timer.restart();
/////////////////////////////////////////////////////////////////
// Read the basic orbital information, without reading all the //
// orbitals themselves.                                        //
/////////////////////////////////////////////////////////////////
//  myComm->barrier();
//  if (myComm->rank() == 1 || myComm->size()==1) // This is because Fortran wants to work with rank 1 and not 0 (or is it?)
  if (myComm->rank() == 0 || myComm->size()==1) 
    if (!ReadOrbitalInfo())
    {
      app_error() << "Error reading orbital info from Espresso interface.  Aborting.\n";
      APP_ABORT("EinsplineSetBuilder::createSPOSet");
    }
  app_log() <<  "TIMER  EinsplineSetBuilder::ReadOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  myComm->barrier();
  orb_info_timer.restart();

//  if (myComm->rank() == 0 || myComm->size()==1) // This is because Fortran wants to work with rank 1 and not 0 (or is it?)
  BroadcastOrbitalInfo();

//  app_log() <<  "TIMER  EinsplineSetBuilder::BroadcastOrbitalInfo " << orb_info_timer.elapsed() << std::endl;
  app_log().flush();

// setup primitive cell and supercell
  PrimCell.set(Lattice);
  SuperCell.set(SuperLattice);
  GGt=dot(transpose(PrimCell.G), PrimCell.G);
//  for (int iat=0; iat<AtomicOrbitals.size(); iat++)
//    AtomicOrbitals[iat].Lattice = Lattice;

// Now, analyze the k-point mesh to figure out the what k-points  are needed
  TwistNum = TwistNum_inp;
  AnalyzeTwists2();
  }/////////////////////////////////////////
}

SPOSet* EinsplineSetBuilderInterface::createSPOSetFromXML(xmlNodePtr cur)
{
  update_token(__FILE__,__LINE__,"createSPOSetFromXML");
//use 2 bohr as the default when truncated orbitals are used based on the extend of the ions
  SPOSet *OrbitalSet;
  int numOrbs = 0;
  int sortBands(1);
  int spinSet = 0;
  int TwistNum_inp=0;

  std::string sourceName;
  std::string spo_prec("double");
  std::string truncate("no");
  std::string hybrid_rep("no");
  std::string use_einspline_set_extended("no"); // use old spline library for high-order derivatives, e.g. needed for backflow optimization
#if defined(QMC_CUDA) || defined(ENABLE_OFFLOAD)
  std::string useGPU="yes";
#else
  std::string useGPU="no";
#endif
  std::string GPUsharing="no";
//  NewTimer* spo_timer = new NewTimer("einspline::CreateSPOSetFromXML", timer_level_medium);
//  TimerManager.addTimer(spo_timer);
  NewTimer* spo_timer = TimerManager.createTimer("einspline::CreateSPOSetFromXML", timer_level_medium);
  spo_timer->start();
// Why the brace here?
  {
    OhmmsAttributeSet a;
    a.add (H5FileName, "href");
    a.add (TileFactor, "tile");
    a.add (sortBands,  "sort");
    a.add (TileMatrix, "tilematrix");
    a.add (TwistNum_inp,   "twistnum");
    a.add (givenTwist,   "twist");
    a.add (sourceName, "source");
    a.add (MeshFactor, "meshfactor");
    a.add (hybrid_rep, "hybridrep");
    a.add (useGPU,     "gpu");
    a.add (GPUsharing, "gpusharing"); // split spline across GPUs visible per rank
    a.add (spo_prec,   "precision");
    a.add (truncate,   "truncate");
    a.add (use_einspline_set_extended,"use_old_spline");
    a.add (myName, "tag");
#if defined(QMC_CUDA)
    a.add (gpu::MaxGPUSpineSizeMB, "Spline_Size_Limit_MB");
#endif

    a.put (XMLRoot);
    a.add (numOrbs,    "size");
    a.add (numOrbs,    "norbs");
    a.add(spinSet,"spindataset"); a.add(spinSet,"group");
    a.put (cur);

    if(myName.empty()) myName="einspline";

  }

  SourcePtcl=ParticleSets[sourceName];
  if(SourcePtcl==0)
  {
    APP_ABORT("Einspline needs the source particleset");
  }
  else
  { //keep the one-body distance table index 
#if defined(ENABLE_SOA)
    myTableIndex=TargetPtcl.addTable(*SourcePtcl,DT_SOA_PREFERRED);
#else
    myTableIndex=TargetPtcl.addTable(*SourcePtcl,DT_AOS);
#endif
    SourcePtcl->addTable(*SourcePtcl,DT_SOA);
  }

///////////////////////////////////////////////
// Read occupation information from XML file //
///////////////////////////////////////////////
  std::vector<int> Occ_Old(0,0);
  Occ.resize(0,0);
  bool NewOcc(false);

  {
    OhmmsAttributeSet oAttrib;
    oAttrib.add(spinSet,"spindataset");
    oAttrib.put(cur);
  }

  xmlNodePtr spo_cur=cur;
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      std::string occ_mode("ground");
      occ_format="energy";
      particle_hole_pairs=0;
      OhmmsAttributeSet oAttrib;
      oAttrib.add(occ_mode,"mode");
      oAttrib.add(spinSet,"spindataset");
      oAttrib.add(occ_format,"format");
      oAttrib.add(particle_hole_pairs,"pairs");
      oAttrib.put(cur);
      if(occ_mode == "excited")
      {
        putContent(Occ,cur);
      }
      else if(occ_mode != "ground")
      {
        app_error() << "Only ground state occupation currently supported "
                    << "in EinsplineSetBuilder.\n";
        APP_ABORT("EinsplineSetBuilder::createSPOSet");
      }
    }
    cur = cur->next;
  }
  if (Occ != Occ_Old)
  {
    NewOcc=true;
    Occ_Old = Occ;
  }
  else
    NewOcc=false;
#if defined(QMC_CUDA)
  if(hybrid_rep=="yes") APP_ABORT("The 'hybridrep' feature of spline SPO has not been enabled on GPU. Stay tuned.");
  app_log() << "\t  QMC_CUDA=1 Overwriting the einspline storage on the host to double precision.\n";
  spo_prec="double"; //overwrite
  truncate="no"; //overwrite
#endif
#if defined(MIXED_PRECISION)
  app_log() << "\t  MIXED_PRECISION=1 Overwriting the einspline storage to single precision.\n";
  spo_prec="single"; //overwrite
#endif
  H5OrbSet aset(H5FileName, spinSet, numOrbs);
  std::map<H5OrbSet,SPOSet*,H5OrbSet>::iterator iter;
  iter = SPOSetMap.find (aset);
  if ((iter != SPOSetMap.end() ) && (!NewOcc))
  {
    app_log() << "SPOSet parameters match in EinsplineSetBuilder:  "
        << "cloning EinsplineSet object.\n";
    return iter->second->makeClone();
  }

  if(FullBands[spinSet]==0) FullBands[spinSet]=new std::vector<BandInfo>;

// Ensure the first SPO set must be spinSet==0
// to correctly initialize key data of EinsplineSetBuilder
  if ( SPOSetMap.size()==0 && spinSet!=0 )
  {
    app_error() << "The first SPO set must have spindataset=\"0\"" << std::endl;
    abort();
  }

// set the internal parameters
  if (spinSet == 0) set_metadata(numOrbs,TwistNum_inp);
//if (use_complex_orb == "yes") UseRealOrbitals = false; // override given user input

// look for <backflow>, would be a lot easier with xpath, but I cannot get it to work
  bool has_backflow = false;

  xmlNodePtr wf  = XMLRoot->parent; // <wavefuntion>
  xmlNodePtr kid = wf->children;
  while (kid != NULL)
  {
    std::string tag((const char*)(kid->name));
    if (tag=="determinantset" || tag=="sposet_builder")
    {
      xmlNodePtr kid1 = kid->children;
      while (kid1 != NULL)
      {
        std::string tag1((const char*)(kid1->name));
        if (tag1=="backflow")
        {
          has_backflow = true;
        }
        kid1 = kid1->next;
      } 
    }
    kid = kid->next; 
  }

  if (has_backflow && use_einspline_set_extended=="yes" && UseRealOrbitals) APP_ABORT("backflow optimization is broken with UseRealOrbitals");

//////////////////////////////////
// Create the OrbitalSet object
//////////////////////////////////
//  app_log() << "So far so good\n";
  Timer mytimer;
  mytimer.restart();
  OccupyBands(spinSet, sortBands, numOrbs);
//  app_log() << "So far so good\n";
  if(spinSet==0) TileIons();
//  app_log() << "So far so good\n";

  bool use_single= (spo_prec == "single" || spo_prec == "float");

// safeguard for a removed feature
  if(truncate=="yes") APP_ABORT("The 'truncate' feature of spline SPO has been removed. Please use hybrid orbital representation.");

//  app_log() << "So far so good\n";
#if !defined(QMC_COMPLEX)
  if (UseRealOrbitals)
  {
//if(TargetPtcl.Lattice.SuperCellEnum != SUPERCELL_BULK && truncate=="yes")
    if(MixedSplineReader==0)
    {
      if(use_single)
        MixedSplineReader= createBsplineRealSingle(this, hybrid_rep=="yes", useGPU);
      else
        MixedSplineReader= createBsplineRealDouble(this, hybrid_rep=="yes", useGPU);
    }
  }
  else
#endif
  {
    if(MixedSplineReader==0)
    {
      if(use_single)
        MixedSplineReader= createBsplineComplexSingle(this, hybrid_rep=="yes", useGPU);
      else
        MixedSplineReader= createBsplineComplexDouble(this, hybrid_rep=="yes", useGPU);
    }
  }

  MixedSplineReader->setCommon(XMLRoot);
// temporary disable the following function call, Ye Luo
// RotateBands_ESHDF(spinSet, dynamic_cast<EinsplineSetExtended<std::complex<double> >*>(OrbitalSet));
  HasCoreOrbs=bcastSortBands(spinSet,NumDistinctOrbitals,myComm->rank()==0);
  SPOSet* bspline_zd=MixedSplineReader->create_spline_set(spinSet,spo_cur);
  if(!bspline_zd)
    APP_ABORT_TRACE(__FILE__,__LINE__,"Failed to create SPOSet*");
  OrbitalSet = bspline_zd;
#if defined(MIXED_PRECISION)
  if(use_einspline_set_extended=="yes")
  {
    app_error() << "Option use_old_spline is not supported by the mixed precision build!" << std::endl;
    abort();
  }
/*
#else
#ifndef QMC_CUDA
if(use_einspline_set_extended=="yes")
#endif
{
EinsplineSet *new_OrbitalSet;
if (UseRealOrbitals)
{
EinsplineSetExtended<double> *temp_OrbitalSet;
#if defined(QMC_CUDA)
if (AtomicOrbitals.size() > 0)
temp_OrbitalSet = new EinsplineSetHybrid<double>;
else
#endif
temp_OrbitalSet = new EinsplineSetExtended<double>;
MixedSplineReader->export_MultiSpline(&(temp_OrbitalSet->MultiSpline));
temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
//set the flags for anti periodic boundary conditions
temp_OrbitalSet->HalfG = dynamic_cast<SplineAdoptorBase<double,3> *>(OrbitalSet)->HalfG;
new_OrbitalSet = temp_OrbitalSet;
}
else
{
EinsplineSetExtended<std::complex<double> > *temp_OrbitalSet;
#if defined(QMC_CUDA)
if (AtomicOrbitals.size() > 0)
temp_OrbitalSet = new EinsplineSetHybrid<std::complex<double> >;
else
#endif
temp_OrbitalSet = new EinsplineSetExtended<std::complex<double> >;
MixedSplineReader->export_MultiSpline(&(temp_OrbitalSet->MultiSpline));
temp_OrbitalSet->MultiSpline->num_splines = NumDistinctOrbitals;
temp_OrbitalSet->resizeStorage(NumDistinctOrbitals, NumValenceOrbs);
for (int iorb=0, num=0; iorb<NumDistinctOrbitals; iorb++)
{
int ti = (*FullBands[spinSet])[iorb].TwistIndex;
temp_OrbitalSet->kPoints[iorb] = PrimCell.k_cart(TwistAngles[ti]);
temp_OrbitalSet->MakeTwoCopies[iorb] = (num < (numOrbs-1)) && (*FullBands[spinSet])[iorb].MakeTwoCopies;
num += temp_OrbitalSet->MakeTwoCopies[iorb] ? 2 : 1;
}
new_OrbitalSet = temp_OrbitalSet;
}
//set the internal parameters
setTiling(new_OrbitalSet,numOrbs);
OrbitalSet = new_OrbitalSet;
}*/
#endif 
  app_log() <<  "Time spent in creating B-spline SPOs " << mytimer.elapsed() << "sec" << std::endl;
#ifdef Ye_debug
#ifndef QMC_COMPLEX
  if (myComm->rank()==0 && OrbitalSet->MuffinTins.size() > 0)
  {
    FILE *fout  = fopen ("TestMuffins.dat", "w");
    Vector<double> phi(numOrbs), lapl(numOrbs);
    Vector<PosType> grad(numOrbs);
    ParticleSet P;
    P.R.resize(6);
    for (int i=0; i<P.R.size(); i++)
      P.R[i] = PosType (0.0, 0.0, 0.0);
    PosType N = 0.25*PrimCell.a(0) + 0.25*PrimCell.a(1) + 0.25*PrimCell.a(2);
    for (double x=-1.0; x<=1.0; x+=0.0000500113412)
    {
// for (double x=-0.003; x<=0.003; x+=0.0000011329343481381) {
      P.R[0] = x * (PrimCell.a(0) + 0.914*PrimCell.a(1) +
	    0.781413*PrimCell.a(2));
      double r = std::sqrt(dot(P.R[0], P.R[0]));
      double rN = std::sqrt(dot(P.R[0]-N, P.R[0]-N));
      OrbitalSet->evaluate(P, 0, phi, grad, lapl);
// OrbitalSet->evaluate(P, 0, phi);
      fprintf (fout, "%1.12e ", r*x/std::abs(x));
      for (int j=0; j<numOrbs; j++)
      {
        double gmag = std::sqrt(dot(grad[j],grad[j]));
        fprintf (fout, "%16.12e ",
	 /*phi[j]*phi[j]**/(-5.0/r  -0.5*lapl[j]/phi[j]));
// double E = -5.0/r -0.5*lapl[j]/phi[j];
        fprintf (fout, "%16.12e ", phi[j]);
        fprintf (fout, "%16.12e ", gmag);
      }
      fprintf (fout, "\n");
    }
    fclose(fout);
  }
#endif
#endif
//if (sourceName.size() && (ParticleSets.find(sourceName) == ParticleSets.end()))
//{
//  app_log() << "  EinsplineSetBuilder creates a ParticleSet " << sourceName << std::endl;
//  ParticleSet* ions=new ParticleSet;
//  ions->Lattice=TargetPtcl.Lattice;
//  ESHDFIonsParser ap(*ions,H5FileID,myComm);
//  ap.put(XMLRoot);
//  ap.expand(TileMatrix);
//  ions->setName(sourceName);
//  ParticleSets[sourceName]=ions;
//  //overwrite the lattice and assign random
  //  if(TargetPtcl.Lattice.SuperCellEnum)
  //  {
  //    TargetPtcl.Lattice=ions->Lattice;
  //    makeUniformRandom(TargetPtcl.R);
  //    TargetPtcl.R.setUnit(PosUnit::LatticeUnit);
  //    TargetPtcl.convert2Cart(TargetPtcl.R);
  //    TargetPtcl.createSK();
  //  }
  //}
#ifdef QMC_CUDA
  if (useGPU == "yes" || useGPU == "1")
  {
    if ((GPUsharing == "yes" || GPUsharing == "1"))
    {
      if (!gpu::cudamps)
      {
        app_log() << "Warning: GPU spline sharing cannot be enabled due to missing Cuda MPS service.\n";
        gpu::device_group_size=1;
      }
      if (gpu::device_group_size>1)
        app_log() << "1/" << gpu::device_group_size << " of GPU spline data stored per rank.\n";
      else
        app_log() << "Full GPU spline data stored per rank.\n";
    } else
    {
      if (gpu::device_group_size>1)
        app_log() << "Full GPU spline data stored per rank.\n";
      gpu::device_group_size=1;
    }
  }
#endif
  OrbitalSet->finalizeConstruction();
  SPOSetMap[aset] = OrbitalSet;
  spo_timer->stop();
  return OrbitalSet;
}
extern bool sortByIndex(BandInfo leftB, BandInfo rightB);

 //!!!!!!!!!!!!!!!!
 /*
 Here there is all the experimental stuff for the hdf5 interface for the current version of qmcpack;
 most of it is taken from the old ion_mover, especially from the ReadOrbitalInfo in EinsplineSetBuilder_Interface.cpp;
 for the moment everything is dumbly locally defined; when (if?) everthing works all the inclusions & references 
 are to be fixed!
 */
bool
 EinsplineSetBuilderInterface::ReadOrbitalInfo()
 {
   //int NumCoreStates , NumMuffinTins , NumTwists , NumSpins , NumBands , NumAtomicOrbitals, NumElectrons;
   //Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv, GGt;
   //Tensor<int,OHMMS_DIM> TileMatrix;
   //std::cerr << "Declare the pointer to ESHDF5interface...";
//  if(myComm->size() > 1)
//    APP_ABORT(" EinsplineSetBuilderInterface::ReadOrbitalInfo(); The QMCQEPack interface at the moment works only for serial runs!\n");
//  OHMMS::Controller->barrier();
///  ESPWSCFInterface* myint;  // To be replaced with the proper type of interface (e.g. pwscf)
///  myint = new ESPWSCFInterface(myComm);
///  esinterface=static_cast<ESInterfaceBase*>(myint);

   esinterface->initialize();


     esinterface->getVersion();
     esinterface->getPrimVecs(Lattice);
     esinterface->getPrimVecs(SuperLattice);
 
//   OHMMS::Controller->barrier();
//   std::cerr << "Barrier: seems to be working\n";

   for(int i=0;i<3;i++)
     for(int j=0;j<3;j++)
       TileMatrix(i,j) = (i==j) ? 1 : 0;
   RecipLattice = 2.0*M_PI*inverse(Lattice);
   SuperLattice = dot(TileMatrix, Lattice);
   char buff[1000];
   snprintf (buff, 1000,
             "  Lattice = \n    [ %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f ]\n",
             Lattice(0,0), Lattice(0,1), Lattice(0,2),
             Lattice(1,0), Lattice(1,1), Lattice(1,2),
             Lattice(2,0), Lattice(2,1), Lattice(2,2));
//   OHMMS::Controller->barrier();
   app_log() << buff;
   snprintf (buff, 1000,
             "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f ]\n",
             SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2),
             SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2),
             SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
//   OHMMS::Controller->barrier();
//   if(myComm->rank()==0)
     CheckLattice();
//   OHMMS::Controller->barrier();
   app_log() << buff;
   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
       LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);
   int have_dpsi = false;
   int NumAtomicOrbitals = 0;
   NumCoreStates = NumMuffinTins = NumTwists = NumSpins = NumBands = NumAtomicOrbitals = 0;
   NumElectrons=TargetPtcl.getTotalNum();
//   OHMMS::Controller->barrier();

   NumBands          = esinterface->getNumBands();
   NumSpins          = esinterface->getNumSpins();
   NumTwists         = esinterface->getNumTwists();
   NumCoreStates     = esinterface->getNumCoreStates();
   NumMuffinTins     = esinterface->getNumMuffinTins();
   have_dpsi         = esinterface->getHaveDPsi();
   NumAtomicOrbitals = esinterface->getNumAtomicOrbitals();
   int num_species   = esinterface->getNumSpecies();
   int NumAtoms      = esinterface->getNumAtoms();
   HaveOrbDerivs = have_dpsi;
   app_log() << "bands=" << NumBands << ", elecs=" << NumElectrons
             << ", spins=" << NumSpins << ", twists=" << NumTwists
             << ", muffin tins=" << NumMuffinTins
             << ", core states=" << NumCoreStates << std::endl;
   app_log() << "atomic orbital=" << NumAtomicOrbitals << std::endl;
   if (TileFactor[0]!=1 || TileFactor[1]!=1 || TileFactor[2]!=1)
     app_log() << "  Using a " << TileFactor[0] << "x" << TileFactor[1]
               << "x" << TileFactor[2] << " tiling factor.\n";

  //////////////////////////////////
  // Read ion types and locations //
  //////////////////////////////////
  app_log()<<"Species ID's"<<std::endl;
  ParticleSet::ParticleIndex_t species_ids;
  esinterface->getSpeciesIDs(species_ids);  
  //int num_species = species_ids.size();  // this is actually the number of atoms
  app_log() << "#Species: " << num_species << std::endl;
  for(int i=0;i<num_species;i++)
    app_log()<<"speciesids: "<<species_ids[i] << std::endl;
  Vector<int> atomic_numbers(num_species);
//  for (int isp=0; isp<num_species; isp++)
//  {
//    std::ostringstream name;
//    name << "/atoms/species_" << isp << "/atomic_number";
//    HDFAttribIO<int> h_atomic_number (atomic_numbers[isp]);
//    h_atomic_number.read(H5FileID, name.str().c_str());
//  }
//  app_log()<<"Get Atomic Numbers...  ";

  esinterface->getAtomicNumbers(atomic_numbers);
//  app_log()<<"Done!\n";
  for (int isp=0; isp<num_species; isp++)
    app_log()<<"speciesids: "<<species_ids[isp] << "\t"<< atomic_numbers[isp] <<std::endl;
  IonTypes.resize(species_ids.size());
  for (int i=0; i<species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
  //HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
  //h_IonPos.read   (H5FileID, "/atoms/positions");

//   for (int i=0; i<IonTypes.size(); i++)
//     app_log() << "Atom type(" << i << ") = " << IonTypes[i] << std::endl;
//   app_log()<<"get Ion Positions"<<std::endl;
   esinterface->getIonPositions(IonPos);
//   app_log() <<"got teh Positions! I think"<<std::endl;
//   for(int i=0;i<NumAtoms;i++)
//     app_log() << i << "\t" <<  IonPos[i][0] << "\t" << IonPos[i][1] << "\t" << IonPos[i][2] << "\n"; 

   //esinterface->getAtomicOrbitals(AtomicOrbitals);

   std::vector<double> dummy(1);
   esinterface->getTwistData(TwistAngles, dummy, TwistSymmetry);

  if(qmc_common.use_density)
  {
    APP_ABORT("So...  Density not implemented yet in interface.  Because I'm lazy")
    //////////////////////////////////////////////////////////
    // Only if it is bulk: If the density has not been set in TargetPtcl, and   //
    // the density is available, read it in and save it     //
    // in TargetPtcl.                                       //
    //////////////////////////////////////////////////////////
    if(TargetPtcl.Lattice.SuperCellEnum == SUPERCELL_BULK)
    {
      // FIXME:  add support for more than one spin density
      if (!TargetPtcl.Density_G.size())
      {
        HDFAttribIO<std::vector<TinyVector<int,OHMMS_DIM> > >
        h_reduced_gvecs(TargetPtcl.DensityReducedGvecs);
        HDFAttribIO<Array<RealType,OHMMS_DIM> >
        h_density_r (TargetPtcl.Density_r);
        TinyVector<int,3> mesh;
        h_reduced_gvecs.read (H5FileID, "/electrons/density/gvectors");
        int numG = TargetPtcl.DensityReducedGvecs.size();
        // Convert primitive G-vectors to supercell G-vectors
        // Also, flip sign since ESHDF format uses opposite sign convention
        #pragma omp parallel for
        for (int iG=0; iG < numG; iG++)
          TargetPtcl.DensityReducedGvecs[iG] =
            -1 * dot(TileMatrix, TargetPtcl.DensityReducedGvecs[iG]);
        app_log() << "  Read " << numG << " density G-vectors.\n";
        for (int ispin=0; ispin<NumSpins; ispin++)
        {
          std::ostringstream density_r_path, density_g_path;
          density_r_path << "/electrons/density/spin_" << ispin << "/density_r";
          density_g_path << "/electrons/density/spin_" << ispin << "/density_g";
          h_density_r.read (H5FileID, density_r_path.str().c_str());
          if (TargetPtcl.DensityReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found density in the HDF5 file.\n";
            std::vector<ComplexType> density_G;
            HDFAttribIO<std::vector<ComplexType > > h_density_G (density_G);
            h_density_G.read (H5FileID, density_g_path.str().c_str());
            if (!density_G.size())
            {
              app_error() << "  Density reduced G-vectors defined, but not the"
                          << " density.\n";
              abort();
            }
            else
            {
              if (ispin == 0)
                TargetPtcl.Density_G = density_G;
              else
                for (int iG=0; iG<density_G.size(); iG++)
                  TargetPtcl.Density_G[iG] += density_G[iG];
            }
          }
        }
      }
      //////////////////////////////////////////////////////////
      // If the density has not been set in TargetPtcl, and   //
      // the density is available, read it in and save it     //
      // in TargetPtcl.                                       //
      //////////////////////////////////////////////////////////
      // FIXME:  add support for more than one spin potential
      if (!TargetPtcl.VHXC_r[0].size())
      {
        HDFAttribIO<std::vector<TinyVector<int,OHMMS_DIM> > >
        h_reduced_gvecs(TargetPtcl.VHXCReducedGvecs);
        TinyVector<int,3> mesh;
        h_reduced_gvecs.read (H5FileID, "/electrons/VHXC/gvectors");
        int numG = TargetPtcl.VHXCReducedGvecs.size();
        // Convert primitive G-vectors to supercell G-vectors
        // Also, flip sign since ESHDF format uses opposite sign convention
        #pragma omp parallel for
        for (int iG=0; iG < numG; iG++)
          TargetPtcl.VHXCReducedGvecs[iG] =
            -1 * dot(TileMatrix, TargetPtcl.VHXCReducedGvecs[iG]);
        app_log() << "  Read " << numG << " VHXC G-vectors.\n";
        for (int ispin=0; ispin<NumSpins; ispin++)
        {
          HDFAttribIO<Array<RealType,OHMMS_DIM> >
          h_VHXC_r (TargetPtcl.VHXC_r[ispin]);
          std::ostringstream VHXC_r_path, VHXC_g_path;
          VHXC_r_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_r";
          VHXC_g_path << "/electrons/VHXC/spin_" << ispin << "/VHXC_g";
          h_VHXC_r.read (H5FileID, VHXC_r_path.str().c_str());
          if (TargetPtcl.VHXCReducedGvecs.size())
          {
            app_log() << "  EinsplineSetBuilder found VHXC in the HDF5 file.\n";
            std::vector<ComplexType> VHXC_G;
            HDFAttribIO<std::vector<ComplexType > > h_VHXC_G (VHXC_G);
            h_VHXC_G.read (H5FileID, VHXC_g_path.str().c_str());
            if (!VHXC_G.size())
            {
              app_error() << "  VHXC reduced G-vectors defined, but not the"
                          << " VHXC.\n";
              abort();
            }
            else
              TargetPtcl.VHXC_G[ispin] = VHXC_G;
          }
        }
      }
    }
  }
  else
  {
    app_log() << "   Skip initialization of the density" << std::endl;
  }

   // AAAAAAAAAAAAAAAAA
//   esinterface->getReducedGVecs(Gvecs,0);
//   exit(0);
//   } // THIS REFERS TO THE RANK=0 CONDITION
   return true;
 }

void EinsplineSetBuilderInterface::OccupyBands(int spin, int sortBands, int numOrbs)
{
  update_token(__FILE__,__LINE__,"OccupyBandsInterface");
  if (myComm->rank() != 0)
    return;

  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  SortBands.clear(); //??? can exit if SortBands is already made?
  int maxOrbs(0);
  app_log()<<"DistincTwists.size()="<<DistinctTwists.size()<<" in OccupyBandsInterface\n";
  app_log()<<" spin="<<spin<<" sortBands="<<sortBands<<" numOrbs="<<numOrbs<<std::endl;
  for (int ti=0; ti<DistinctTwists.size(); ti++)
  {
    app_log()<<"Doop\n";
    app_log()<<"DistinctTwist[ti]="<<ti<<" "<<DistinctTwists[ti]<<std::endl;

    int tindex = DistinctTwists[ti];
    std::vector<double> eigvals;
    esinterface->getOrbEigenvals(spin,ti,eigvals);
    for (int bi=0; bi<NumBands; bi++)
    {
      BandInfo band;
      band.IsCoreState = false;
      band.TwistIndex = tindex;
      band.BandIndex  = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      band.Energy = eigvals[bi]; 
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs+=2;
      else
        maxOrbs++;
    }
  }
   // Now sort the bands by energy
  if (sortBands==2)
  {
    app_log() << "Sorting the bands by index now:\n";
    sort (SortBands.begin(), SortBands.end(), sortByIndex);
  }
  else if (sortBands==1)
  {
    app_log() << "Sorting the bands now:\n";
    sort (SortBands.begin(), SortBands.end());
  }

  std::vector<int> gsOcc(maxOrbs);
  int N_gs_orbs=numOrbs;
  int nocced(0), ntoshift(0);
  for (int ti=0; ti<SortBands.size(); ti++)
  {
    if (nocced<N_gs_orbs)
    {
      if (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced>1))
      {
        nocced+=2;
        ntoshift++;
        gsOcc[ti]=2;
      }
      else if ( (SortBands[ti].MakeTwoCopies && (N_gs_orbs-nocced==1)) || !SortBands[ti].MakeTwoCopies )
      {
        nocced+=1;
        ntoshift++;
        gsOcc[ti]=1;
      }
    }
  }
  if (occ_format=="energy")
  {
    // To get the occupations right.
    std::vector<int> Removed(0,0);
    std::vector<int> Added(0,0);
    for(int ien=0; ien<Occ.size(); ien++)
    {
      if (Occ[ien]<0)
        Removed.push_back(-Occ[ien]);
      else if (Occ[ien]>0)
        Added.push_back(Occ[ien]);
    }
    if(Added.size()-Removed.size() != 0)
    {
      app_log()<<"need to add and remove same number of orbitals. "<< Added.size()<<" "<<Removed.size()<<std::endl;
      APP_ABORT("ChangedOccupations");
    }
    std::vector<int> DiffOcc(maxOrbs,0);
    //Probably a cleaner way to do this.
    for(int i=0; i<Removed.size(); i++)
      DiffOcc[Removed[i]-1]-=1;
    for(int i=0; i<Added.size(); i++)
      DiffOcc[Added[i]-1]+=1;
    std::vector<int> SumOrb(SortBands.size(),0);
    int doi(0);
    for(int i=0; i<SumOrb.size(); i++)
    {
      if(SortBands[i].MakeTwoCopies)
      {
        SumOrb[i]=  gsOcc[i]+DiffOcc[doi++];
        SumOrb[i]+= DiffOcc[doi++];
      }
      else
        SumOrb[i]=gsOcc[i]+DiffOcc[doi++];
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for(int i=0; i<SumOrb.size(); i++)
    {
      if(SumOrb[i]==2)
      {
        SortBands[i].MakeTwoCopies=true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i]==1)
      {
        SortBands[i].MakeTwoCopies=false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (SumOrb[i]==0)
      {
        SortBands[i].MakeTwoCopies=false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<<std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
    SortBands=ReOrderedBands;
  }
  else if (occ_format=="band")
  {
    app_log()<<"  Occupying bands based on (bi,ti) data."<<std::endl;
    if(Occ.size() != particle_hole_pairs*4)
    {
      app_log()<<" Need Occ = pairs*4. Occ is (ti,bi) of removed, then added."<<std::endl;
      app_log()<<Occ.size()<<" "<<particle_hole_pairs<<std::endl;
      APP_ABORT("ChangedOccupations");
    }
    int cnt(0);
    for(int ien=0; ien<SortBands.size(); ien++)
    {
      if((Occ[cnt] == SortBands[ien].TwistIndex)&&(Occ[cnt+1] == SortBands[ien].BandIndex))
        if(cnt<particle_hole_pairs*2)
        {
          gsOcc[ien]-=1;
          cnt+=2;
          app_log()<<"removing orbital "<<ien<<std::endl;
        }
        else
        {
          gsOcc[ien]+=1;
          app_log()<<"adding orbital "<<ien<<std::endl;
          cnt+=2;
        }
    }
    std::vector<BandInfo> ReOrderedBands;
    std::vector<BandInfo> RejectedBands;
    for(int i=0; i<SortBands.size(); i++)
    {
      if(gsOcc[i]==2)
      {
        SortBands[i].MakeTwoCopies=true;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i]==1)
      {
        SortBands[i].MakeTwoCopies=false;
        ReOrderedBands.push_back(SortBands[i]);
      }
      else if (gsOcc[i]==0)
      {
        SortBands[i].MakeTwoCopies=false;
        RejectedBands.push_back(SortBands[i]);
      }
      else
      {
        app_log()<<" Trying to add the same orbital ("<<i<<") less than zero or more than 2 times."<<std::endl;
        APP_ABORT("Sorting Excitation");
      }
    }
    ReOrderedBands.insert(ReOrderedBands.end(),RejectedBands.begin(),RejectedBands.end());
    SortBands=ReOrderedBands;
  }
  //for(int sw=0;sw<Removed.size();sw++){
  //  app_log()<<" Swapping two orbitals "<<Removed[sw]<<" and "<<Added[sw]<<std::endl;
  //  BandInfo tempband(SortBands[Removed[sw]-1]);
  //  SortBands[Removed[sw]-1] = SortBands[Added[sw]-1];
  //  SortBands[Added[sw]-1] = tempband;
  //}
  int orbIndex = 0;
  int numOrbs_counter = 0;
  NumValenceOrbs=0;
  NumCoreOrbs=0;
  while (numOrbs_counter < numOrbs)
  {
    if (SortBands[orbIndex].MakeTwoCopies)
      numOrbs_counter += 2;
    else
      numOrbs_counter++;
    if (SortBands[orbIndex].IsCoreState)
      NumCoreOrbs++;
    else
      NumValenceOrbs++;
    orbIndex++;
  }
  NumDistinctOrbitals = orbIndex;
  app_log() << "We will read " << NumDistinctOrbitals << " distinct orbitals.\n";
  app_log() << "There are " << NumCoreOrbs << " core states and "
            << NumValenceOrbs << " valence states.\n";
}

bool EinsplineSetBuilderInterface::TwistPair (PosType a, PosType b)
{
  bool pair = true;
  for (int n=0; n<OHMMS_DIM; n++)
  {
    double d = a[n] + b[n];
    if (std::abs(d - round(d)) > MatchingTol)
      pair = false;
  }
  return pair;
}

void
EinsplineSetBuilderInterface::AnalyzeTwists2()
{
  Tensor<double,3> S;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      S(i,j) = (double)TileMatrix(i,j);
  std::vector<PosType> superFracs;
  // This holds to which supercell kpoint each primitive k-point belongs
  std::vector<int> superIndex;
  int numPrimTwists = TwistAngles.size();
  for (int ki=0; ki<numPrimTwists; ki++)
  {
    PosType primTwist = TwistAngles[ki];
    PosType superTwist = dot (S, primTwist);
    PosType kp = PrimCell.k_cart(primTwist);
    PosType ks = SuperCell.k_cart(superTwist);
    if (dot(ks-kp, ks-kp) > 1.0e-6)
    {
      app_error() << "Primitive and super k-points do not agree.  Error in coding.\n";
      app_error().flush();
      APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
    }
    PosType frac = FracPart (superTwist);
    bool found = false;
    for (int j=0; j<superFracs.size(); j++)
    {
      PosType diff = frac - superFracs[j];
      if (dot(diff,diff)<1.0e-6)
      {
        found = true;
        superIndex.push_back(j);
      }
    }
    if (!found)
    {
      superIndex.push_back(superFracs.size());
      superFracs.push_back(frac);
    }
  }
  int numSuperTwists = superFracs.size();
  app_log() << "Found " << numSuperTwists << " distinct supercell twists.\n";
  // For each supercell twist, create a list of primitive twists which
  // belong to it.
  std::vector<std::vector<int> > superSets;
  superSets.resize(numSuperTwists);
  for (int ki=0; ki<numPrimTwists; ki++)
    superSets[superIndex[ki]].push_back(ki);
  app_log()<<"number of things"<< std::endl;
  app_log()<<TwistSymmetry.size()<< std::endl;
  app_log()<<TwistWeight.size()<< std::endl;
//     for (int ki=0; ki<TwistSymmetry.size(); ki++)
//       fprintf (stderr, "%d %d %d\n",ki,TwistSymmetry[ki],TwistWeight[ki]);
  if (myComm->rank() == 0)
  {
    int n_tot_irred(0);
    for (int si=0; si<numSuperTwists; si++)
    {
//      bool irreducible(false);
      int irrep_wgt(0);
// 	 for (int i=0; i<superSets[si].size(); i++)
      if(TwistSymmetry[superSets[si][0]]==1)
      {
        irrep_wgt=TwistWeight[superSets[si][0]];
//        irreducible=true;
        n_tot_irred++;
      }
//	if((irreducible) and ((Version[0] >= 2) and (Version[1] >= 0)))
//	  fprintf (stderr, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]  IRREDUCIBLE-K %d  %d \n",
// 		 si, superFracs[si][0], superFracs[si][1], superFracs[si][2], si, irrep_wgt);
//	else
      char buf[1000];
      snprintf (buf, 1000, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n",
               si, superFracs[si][0], superFracs[si][1], superFracs[si][2]);
      app_log() << buf;
      app_log().flush();
//  	fprintf (stderr, "  Using k-points: ");
//  	for (int i=0; i<superSets[si].size(); i++)
//  	  fprintf (stderr, " %d", superSets[si][i]);
//  	fprintf (stderr, "\n");
    }
//        fprintf (stderr, "Number in irredicible twist grid: %d \n", n_tot_irred);
  }
  if(TwistNum<0)
  {
    PosType gtFrac = FracPart(givenTwist);
    float eps=1e-5;
    for (int si=0; si<numSuperTwists; si++) {
      PosType locDiff = gtFrac - superFracs[si];
      if (dot(locDiff,locDiff) < eps)
	TwistNum=si;
    }
  }
  // Check supertwist for this node
  if (!myComm->rank()) {
    char buf[1000];
    snprintf (buf, 1000, "  Using supercell twist %d:  [ %9.5f %9.5f %9.5f]\n",
             TwistNum, superFracs[TwistNum][0], superFracs[TwistNum][1],
             superFracs[TwistNum][2]);
    app_log() << buf;
    app_log().flush();
  }
  TargetPtcl.setTwist(superFracs[TwistNum]);
#ifndef QMC_COMPLEX
  // Check to see if supercell twist is okay to use with real wave
  // functions
  for (int dim=0; dim<OHMMS_DIM; dim++)
  {
    double t = 2.0*superFracs[TwistNum][dim];
    if (std::abs(t - round(t)) > MatchingTol*100)
    {
      app_error() << "Cannot use this super twist with real wavefunctions.\n"
                  << "Please recompile with QMC_COMPLEX=1.\n";
      app_error().flush();
      APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");
    }
  }
#endif
  // Now check to see that each supercell twist has the right twists
  // to tile the primitive cell orbitals.
  int numTwistsNeeded = std::abs(det(TileMatrix));
  for (int si=0; si<numSuperTwists; si++)
  {
    // First make sure we have enough points
    if (superSets[si].size() != numTwistsNeeded)
    {
      char buf[1000];
      snprintf (buf, 1000, "Super twist %d should own %d k-points, but owns %d.\n",
               si, numTwistsNeeded, static_cast<int>(superSets[si].size()));
      app_error() << buf;
      if(si==TwistNum)
        {APP_ABORT("EinsplineSetBuilder::AnalyzeTwists2");}
      else
        continue;
    }
    // Now, make sure they are all distinct
    int N = superSets[si].size();
    for (int i=0; i<N; i++)
    {
      PosType twistPrim_i  = TwistAngles[superSets[si][i]];
      PosType twistSuper_i = dot (S, twistPrim_i);
      PosType superInt_i   = IntPart (twistSuper_i);
      for (int j=i+1; j<N; j++)
      {
        PosType twistPrim_j  = TwistAngles[superSets[si][j]];
        PosType twistSuper_j = dot (S, twistPrim_j);
        PosType superInt_j   = IntPart (twistSuper_j);
        if (dot(superInt_i-superInt_j, superInt_i-superInt_j) < 1.0e-6)
        {
          app_error() << "Identical k-points detected in super twist set "
                      << si << std::endl;
          APP_ABORT_TRACE(__FILE__,__LINE__, "AnalyzeTwists2");
        }
      }
    }
  }
  app_log().flush();
  if (TwistNum >= superSets.size())
  {
    app_error() << "Trying to use supercell twist " << TwistNum
                << " when only " << superSets.size() << " sets exist.\n"
                << "Please select a twist number between 0 and "
                << superSets.size()-1 << ".\n";
    APP_ABORT_TRACE(__FILE__,__LINE__, "Invalid TwistNum");
  }
  // Finally, record which k-points to include on this group of
  // processors, which have been assigned supercell twist TwistNum
  IncludeTwists.clear();
  for (int i=0; i<superSets[TwistNum].size(); i++)
    IncludeTwists.push_back(superSets[TwistNum][i]);
  // Now, find out which twists are distinct
  DistinctTwists.clear();
#ifndef QMC_COMPLEX
  std::vector<int> copyTwists;
  for (int i=0; i<IncludeTwists.size(); i++)
  {
    int ti        = IncludeTwists[i];
    PosType twist_i = TwistAngles[ti];
    bool distinct=true;
    for (int j=i+1; j<IncludeTwists.size(); j++)
    {
      int tj = IncludeTwists[j];
      PosType twist_j = TwistAngles[tj];
      PosType sum  = twist_i + twist_j;
      PosType diff = twist_i - twist_j;
      if (TwistPair (twist_i, twist_j))
        distinct = false;
    }
    if (distinct)
      DistinctTwists.push_back(ti);
    else
      copyTwists.push_back(ti);
  }
  // Now determine which distinct twists require two copies
  MakeTwoCopies.resize(DistinctTwists.size());
  for (int i=0; i<DistinctTwists.size(); i++)
  {
    MakeTwoCopies[i] = false;
    int ti = DistinctTwists[i];
    PosType twist_i = TwistAngles[ti];
    for (int j=0; j<copyTwists.size(); j++)
    {
      int tj = copyTwists[j];
      PosType twist_j = TwistAngles[tj];
      if (TwistPair(twist_i, twist_j))
        MakeTwoCopies[i] = true;
    }
    if (myComm->rank() == 0) {
      char buf[1000];
      
      snprintf (buf, 1000, "Using %d copies of twist angle [%6.3f, %6.3f, %6.3f]\n",
		MakeTwoCopies[i] ? 2 : 1, twist_i[0], twist_i[1], twist_i[2]);
      app_log() << buf;
      app_log().flush();
    }
  }
  // Find out if we can make real orbitals
  UseRealOrbitals = true;
  for (int i=0; i < DistinctTwists.size(); i++)
  {
    int ti = DistinctTwists[i];
    PosType twist = TwistAngles[ti];
    for (int j=0; j<OHMMS_DIM; j++)
      if (std::abs(twist[j]-0.0) > MatchingTol &&
          std::abs(twist[j]-0.5) > MatchingTol &&
          std::abs(twist[j]+0.5) > MatchingTol)
        UseRealOrbitals = false;
  }
  if (UseRealOrbitals && (DistinctTwists.size() > 1))
  {
    app_log() << "***** Use of real orbitals is possible, but not currently implemented\n"
              << "      with more than one twist angle.\n";
    UseRealOrbitals = false;
  }
  if (UseRealOrbitals)
    app_log() << "Using real orbitals.\n";
  else
    app_log() << "Using complex orbitals.\n";
#else
  DistinctTwists.resize(IncludeTwists.size());
  MakeTwoCopies.resize(IncludeTwists.size());
  for (int i=0; i<IncludeTwists.size(); i++)
  {
    DistinctTwists[i] = IncludeTwists[i];
    MakeTwoCopies[i] = false;
  }
  UseRealOrbitals = false;
#endif
}

void
EinsplineSetBuilderInterface::TileIons()
{
  update_token(__FILE__,__LINE__,"TileIons");

  //set the primitive lattice
  SourcePtcl->PrimitiveLattice.set(Lattice);

  for(int j=0; j<IonPos.size(); ++j)
    IonPos[j]=FracPart(SourcePtcl->PrimitiveLattice.toUnit(IonPos[j]));

  for(int i=0; i<SourcePtcl->R.size(); ++i)
  {
    PosType u=FracPart(SourcePtcl->PrimitiveLattice.toUnit(SourcePtcl->R[i]));
    int j=0;
    bool foundit=false;
    while(!foundit &&j<IonPos.size())
    {
      PosType d=u-IonPos[j++];
      foundit=(dot(d,d)<MatchingTol);
    }
    SourcePtcl->PCID[i]=j-1;
  }

  IonPos.resize(SourcePtcl->getTotalNum());
  IonTypes.resize(SourcePtcl->getTotalNum());
  std::copy(SourcePtcl->R.begin(),SourcePtcl->R.end(),IonPos.begin());
  std::copy(SourcePtcl->GroupID.begin(),SourcePtcl->GroupID.end(),IonTypes.begin());

  //app_log() << "  Primitive Cell\n";
  //SourcePtcl->PrimitiveLattice.print(app_log());
  //app_log() << "  Super Cell\n";
  //SourcePtcl->Lattice.print(app_log());

//Don't need to do this, already one by ParticleSetPool.cpp
//  Vector<TinyVector<double, OHMMS_DIM> > primPos   = IonPos;
//  Vector<int>                            primTypes = IonTypes;
//  int numCopies = std::abs(det(TileMatrix));
//  IonTypes.resize(primPos.size()*numCopies);
//  IonPos.resize  (primPos.size()*numCopies);
//  int maxCopies = 10;
//  typedef TinyVector<double,3> Vec3;
//  int index=0;
//  for (int i0=-maxCopies; i0<=maxCopies; i0++)
//    for (int i1=-maxCopies; i1<=maxCopies; i1++)
//      for (int i2=-maxCopies; i2<=maxCopies; i2++)
//        for (int iat=0; iat < primPos.size(); iat++)
//        {
//          Vec3 r     = primPos[iat];
//          Vec3 uPrim = PrimCell.toUnit(r);
//          for (int i=0; i<3; i++)
//            uPrim[i] -= std::floor(uPrim[i]);
//          r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0) +
//              (double)i1*PrimCell.a(1) + (double)i2*PrimCell.a(2);
//          Vec3 uSuper = SuperCell.toUnit(r);
//          if ((uSuper[0] >= -1.0e-4) && (uSuper[0] < 0.9999) &&
//              (uSuper[1] >= -1.0e-4) && (uSuper[1] < 0.9999) &&
//              (uSuper[2] >= -1.0e-4) && (uSuper[2] < 0.9999))
//          {
//            IonPos[index]= r;
//            IonTypes[index]= primTypes[iat];
//            index++;
//          }
//        }
//  if (index != primPos.size()*numCopies)
//  {
//    app_error() << "The number of tiled ions, " << IonPos.size()
//                << ", does not match the expected number of "
//                << primPos.size()*numCopies << " or the index "<< index <<".  Aborting.\n";
//    APP_ABORT("EinsplineSetBuilder::TileIons()");
//  }
//  if (myComm->rank() == 0)
//  {
//    char buf[1000];
//    snprintf (buf, 1000, "Supercell reduced ion positions = \n");
//    app_log() << buf;
//    app_log().flush();
//    for (int i=0; i<IonPos.size(); i++)
//    {
//      PosType u = SuperCell.toUnit(IonPos[i]);
//      char buf2[1000];
//      snprintf (buf2, 1000, "   %14.10f %14.10f %14.10f\n",
//               u[0], u[1], u[2]);
//      app_log() << buf2;
//      app_log().flush();
//      //		 IonPos[i][0], IonPos[i][1], IonPos[i][2]);
//    }
//  }
}
  inline bool EinsplineSetBuilderInterface::bcastSortBands(int spin, int n, bool root)
  {
    update_token(__FILE__,__LINE__,"bcastSortBands");

    std::vector<BandInfo>& SortBands(*FullBands[spin]);

    TinyVector<int,4> nbands(int(SortBands.size()),n,NumValenceOrbs,NumCoreOrbs);
    mpi::bcast(*myComm,nbands);

    //buffer to serialize BandInfo
    PooledData<OHMMS_PRECISION_FULL> misc(nbands[0]*5);
    bool isCore=false;
    n=NumDistinctOrbitals=nbands[1];
    NumValenceOrbs=nbands[2];
    NumCoreOrbs=nbands[3];

    if(root)
    {
      misc.rewind();
      //misc.put(NumValenceOrbs);
      //misc.put(NumCoreOrbs);
      for(int i=0; i<n; ++i)
      {
        misc.put(SortBands[i].TwistIndex);
	misc.put(SortBands[i].BandIndex);
	misc.put(SortBands[i].Energy);
        misc.put(SortBands[i].MakeTwoCopies);
        misc.put(SortBands[i].IsCoreState);

        isCore |= SortBands[i].IsCoreState;
      }

      for(int i=n; i<SortBands.size(); ++i)
      {
        misc.put(SortBands[i].TwistIndex);
	misc.put(SortBands[i].BandIndex);
	misc.put(SortBands[i].Energy);
        misc.put(SortBands[i].MakeTwoCopies);
        misc.put(SortBands[i].IsCoreState);
      }
    }
    myComm->bcast(misc);

    if(!root)
    {
      SortBands.resize(nbands[0]);
      misc.rewind();
      //misc.get(NumValenceOrbs);
      //misc.get(NumCoreOrbs);
      for(int i=0; i<n; ++i)
      {
        misc.get(SortBands[i].TwistIndex);
	misc.get(SortBands[i].BandIndex);
	misc.get(SortBands[i].Energy);
        misc.get(SortBands[i].MakeTwoCopies);
        misc.get(SortBands[i].IsCoreState);

        isCore |= SortBands[i].IsCoreState;
      }
      for(int i=n; i<SortBands.size(); ++i)
      {
        misc.get(SortBands[i].TwistIndex);
	misc.get(SortBands[i].BandIndex);
	misc.get(SortBands[i].Energy);
        misc.get(SortBands[i].MakeTwoCopies);
        misc.get(SortBands[i].IsCoreState);
      }
    }

    //char fname[64];
    //sprintf(fname,"debug.%d",myComm->rank());
    //ofstream fout(fname);
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //for(int i=0; i<misc.size();++i)
    //  fout << misc[i] << std::endl;
    return isCore;
  }
bool EinsplineSetBuilderInterface::ReadGvectors_ESHDF()
{
  std::cerr << "YOU SHOULD NOT BE HERE \n";
  update_token(__FILE__,__LINE__,"ReadGvectors_ESHDF");
  bool root=myComm->rank() ==0;
  //this is always ugly
  MeshSize = 0;
  int hasPsig=1;
#if defined(__bgp__)||(__bgq__)
  if(root)
  {
    hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_g");
    if(gid<0)
      hasPsig=0;
    H5Dclose(gid);
  }
  myComm->bcast(hasPsig);
#else
  if(root)
  {
    HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize);
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
    h_mesh.read (H5FileID, "/electrons/mesh");
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
#endif
  if(hasPsig)
  {
    int nallowed=257;
    int allowed[] =
    {
      72,75,80,81,90,96,100,108,120,125,
      128,135,144,150,160,162,180,192,200,216,
      225,240,243,250,256,270,288,300,320,324,
      360,375,384,400,405,432,450,480,486,500,
      512,540,576,600,625,640,648,675,720,729,
      750,768,800,810,864,900,960,972,1000,1024,
      1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,
      1458,1500,1536,1600,1620,1728,1800,1875,1920,1944,
      2000,2025,2048,2160,2187,2250,2304,2400,2430,2500,
      2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,
      3375,3456,3600,3645,3750,3840,3888,4000,4050,4096,
      4320,4374,4500,4608,4800,4860,5000,5120,5184,5400,
      5625,5760,5832,6000,6075,6144,6250,6400,6480,6561,
      6750,6912,7200,7290,7500,7680,7776,8000,8100,8192,
      8640,8748,9000,9216,9375,9600,9720,10000,10125,10240,
      10368,10800,10935,11250,11520,11664,12000,12150,12288,12500,
      12800,12960,13122,13500,13824,14400,14580,15000,15360,15552,
      15625,16000,16200,16384,16875,17280,17496,18000,18225,18432,
      18750,19200,19440,19683,20000,20250,20480,20736,21600,21870,
      22500,23040,23328,24000,24300,24576,25000,25600,25920,26244,
      27000,27648,28125,28800,29160,30000,30375,30720,31104,31250,
      32000,32400,32768,32805,33750,34560,34992,36000,36450,36864,
      37500,38400,38880,39366,40000,40500,40960,41472,43200,43740,
      45000,46080,46656,46875,48000,48600,49152,50000,50625,51200,
      51840,52488,54000,54675,55296,56250,57600,58320,59049,60000,
      60750,61440,62208,62500,64000,64800,65536
    };
    MaxNumGvecs=0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int,3> maxIndex(0,0,0);
    Gvecs.resize(NumTwists);
    {
      int numg=0;
      if(root)
      {
        std::ostringstream Gpath;
        Gpath    << "/electrons/kpoint_0/gvectors";
        HDFAttribIO<std::vector<TinyVector<int,3> > > h_Gvecs(Gvecs[0]);
        h_Gvecs.read (H5FileID, Gpath.str().c_str());
        numg=Gvecs[0].size();
      }
      myComm->bcast(numg);
      if(!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs=Gvecs[0].size();
      for (int ig=0; ig<Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0*MeshFactor*maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0*MeshFactor*maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0*MeshFactor*maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int *ix=std::lower_bound(allowed,allowed+nallowed,MeshSize[0]);
    int *iy=std::lower_bound(allowed,allowed+nallowed,MeshSize[1]);
    int *iz=std::lower_bound(allowed,allowed+nallowed,MeshSize[2]);
    MeshSize[0]=(MeshSize[0]>128)? *ix:(MeshSize[0]+MeshSize[0]%2);
    MeshSize[1]=(MeshSize[1]>128)? *iy:(MeshSize[1]+MeshSize[1]%2);
    MeshSize[2]=(MeshSize[2]>128)? *iz:(MeshSize[2]+MeshSize[2]%2);
    if(Version[0]<2)
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE." << std::endl;
      for(int k=0; k<DistinctTwists.size(); ++k)
      {
        int ik=DistinctTwists[k];
        if(ik==0)
          continue; //already done
        int numg=0;
        if(root)
        {
          std::ostringstream Gpath;
          Gpath    << "/electrons/kpoint_" << ik << "/gvectors";
          HDFAttribIO<std::vector<TinyVector<int,3> > > h_Gvecs(Gvecs[ik]);
          h_Gvecs.read (H5FileID, Gpath.str().c_str());
          numg=Gvecs[ik].size();
        }
        myComm->bcast(numg);
        if(numg==0)
        {
          //copy kpoint_0, default
          Gvecs[ik]=Gvecs[0];
        }
        else
        {
          if(numg !=  MaxNumGvecs)
          {
            std::ostringstream o;
            o<< "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
             << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if(!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}

bool EinsplineSetBuilderInterface::ReadGvectors_PWSCF()
{
  update_token(__FILE__,__LINE__,"ReadGvectors_PWSCF");
  bool root=myComm->rank() ==0;
  //this is always ugly
  MeshSize = 0;
  int hasPsig=1;
/*
#if defined(__bgp__)||(__bgq__)
  if(root)
  {
    hid_t gid=H5Dopen(H5FileID,"/electrons/kpoint_0/spin_0/state_0/psi_g");
    if(gid<0)
      hasPsig=0;
    H5Dclose(gid);
  }
  myComm->bcast(hasPsig);
#else
  if(root)
  {
    HDFAttribIO<TinyVector<int,3> > h_mesh(MeshSize);
    h_mesh.read (H5FileID, "/electrons/psi_r_mesh");
    h_mesh.read (H5FileID, "/electrons/mesh");
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
#endif
*/
  if(hasPsig)
  {
    int nallowed=257;
    int allowed[] =
    {
      72,75,80,81,90,96,100,108,120,125,
      128,135,144,150,160,162,180,192,200,216,
      225,240,243,250,256,270,288,300,320,324,
      360,375,384,400,405,432,450,480,486,500,
      512,540,576,600,625,640,648,675,720,729,
      750,768,800,810,864,900,960,972,1000,1024,
      1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,
      1458,1500,1536,1600,1620,1728,1800,1875,1920,1944,
      2000,2025,2048,2160,2187,2250,2304,2400,2430,2500,
      2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,
      3375,3456,3600,3645,3750,3840,3888,4000,4050,4096,
      4320,4374,4500,4608,4800,4860,5000,5120,5184,5400,
      5625,5760,5832,6000,6075,6144,6250,6400,6480,6561,
      6750,6912,7200,7290,7500,7680,7776,8000,8100,8192,
      8640,8748,9000,9216,9375,9600,9720,10000,10125,10240,
      10368,10800,10935,11250,11520,11664,12000,12150,12288,12500,
      12800,12960,13122,13500,13824,14400,14580,15000,15360,15552,
      15625,16000,16200,16384,16875,17280,17496,18000,18225,18432,
      18750,19200,19440,19683,20000,20250,20480,20736,21600,21870,
      22500,23040,23328,24000,24300,24576,25000,25600,25920,26244,
      27000,27648,28125,28800,29160,30000,30375,30720,31104,31250,
      32000,32400,32768,32805,33750,34560,34992,36000,36450,36864,
      37500,38400,38880,39366,40000,40500,40960,41472,43200,43740,
      45000,46080,46656,46875,48000,48600,49152,50000,50625,51200,
      51840,52488,54000,54675,55296,56250,57600,58320,59049,60000,
      60750,61440,62208,62500,64000,64800,65536
    };
    MaxNumGvecs=0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int,3> maxIndex(0,0,0);
    Gvecs.resize(NumTwists);
    {
      int numg=0;
      if(root)
      {
        esinterface->getReducedGVecs(Gvecs,0);
        numg=Gvecs[0].size();
      }
      myComm->bcast(numg);
      if(!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs=Gvecs[0].size();
      for (int ig=0; ig<Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0*MeshFactor*maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0*MeshFactor*maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0*MeshFactor*maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int *ix=std::lower_bound(allowed,allowed+nallowed,MeshSize[0]);
    int *iy=std::lower_bound(allowed,allowed+nallowed,MeshSize[1]);
    int *iz=std::lower_bound(allowed,allowed+nallowed,MeshSize[2]);
    MeshSize[0]=(MeshSize[0]>128)? *ix:(MeshSize[0]+MeshSize[0]%2);
    MeshSize[1]=(MeshSize[1]>128)? *iy:(MeshSize[1]+MeshSize[1]%2);
    MeshSize[2]=(MeshSize[2]>128)? *iz:(MeshSize[2]+MeshSize[2]%2);
    if(Version[0]<2)
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  PWSCF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE." << std::endl;
      for(int k=0; k<DistinctTwists.size(); ++k)
      {
        int ik=DistinctTwists[k];
        if(ik==0)
          continue; //already done
        int numg=0;
        if(root)
        {
          esinterface->getReducedGVecs(Gvecs,0);
          numg=Gvecs[ik].size();
        }
        myComm->bcast(numg);
        if(numg==0)
        {
          //copy kpoint_0, default
          Gvecs[ik]=Gvecs[0];
        }
        else
        {
          if(numg !=  MaxNumGvecs)
          {
            std::ostringstream o;
            o<< "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
             << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if(!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}
}

