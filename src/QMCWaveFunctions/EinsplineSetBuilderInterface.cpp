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
    
    

#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <vector>
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "ParticleIO/ESHDFParticleParser.h"
#include "ParticleBase/RandomSeqGenerator.h"

#include "qmc_common.h"

namespace qmcplusplus
{
extern bool sortByIndex(BandInfo leftB, BandInfo rightB);

 //!!!!!!!!!!!!!!!!
 /*
 Here there is all the experimental stuff for the hdf5 interface for the current version of qmcpack;
 most of it is taken from the old ion_mover, especially from the ReadOrbitalInfo in EinsplineSetBuilder_Interface.cpp;
 for the moment everything is dumbly locally defined; when (if?) everthing works all the inclusions & references 
 are to be fixed!
 */
bool
 EinsplineSetBuilder::ReadOrbitalInfo_Interface()
 {
   //int NumCoreStates , NumMuffinTins , NumTwists , NumSpins , NumBands , NumAtomicOrbitals, NumElectrons;
   //Tensor<double,OHMMS_DIM> Lattice, RecipLattice, LatticeInv, GGt;
   //Tensor<int,OHMMS_DIM> TileMatrix;

   std::cerr << "Declare the pointer to ESHDF5interface...";
   ESHDF5Interface esinterface;
   std::cerr << "  Done!\n";

   std::cerr << "Initialize it... \n";
   esinterface.initialize();
   std::cerr << "  Done!\n";

   esinterface.getVersion();

   std::cerr << "Getting primitive vectors...";
   esinterface.getPrimVecs(Lattice);
   std::cerr << "  Done!\n";
 
   OHMMS::Controller->barrier();
   std::cerr << "Barrier: seems to be working\n";

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
   OHMMS::Controller->barrier();
   app_log() << buff;
   snprintf (buff, 1000,
             "  SuperLattice = \n    [ %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f\n"
             "      %9.6f %9.6f %9.6f ]\n",
             SuperLattice(0,0), SuperLattice(0,1), SuperLattice(0,2),
             SuperLattice(1,0), SuperLattice(1,1), SuperLattice(1,2),
             SuperLattice(2,0), SuperLattice(2,1), SuperLattice(2,2));
   CheckLattice();
   OHMMS::Controller->barrier();

   app_log() << buff;
   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
       LatticeInv(i,j) = RecipLattice(i,j)/(2.0*M_PI);
   int have_dpsi = false;
   int NumAtomicOrbitals = 0;
   NumCoreStates = NumMuffinTins = NumTwists = NumSpins = NumBands = NumAtomicOrbitals = 0;
   NumElectrons=TargetPtcl.getTotalNum();

   OHMMS::Controller->barrier();

   NumBands          = esinterface.getNumBands();
   NumSpins          = esinterface.getNumSpins();
   NumTwists         = esinterface.getNumTwists();
   NumCoreStates     = esinterface.getNumCoreStates();
   NumMuffinTins     = esinterface.getNumMuffinTins();
   have_dpsi         = esinterface.getHaveDPsi();
   NumAtomicOrbitals = esinterface.getNumAtomicOrbitals();

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
  esinterface.getSpeciesIDs(species_ids);
  int num_species = species_ids.size();
  app_log() << "#Species: " << num_species << std::endl;
  for(int i=0;i<num_species;i++)
    app_log()<<"speciesids: "<<species_ids[i] << std::endl;
  Vector<int> atomic_numbers(num_species);
  for (int isp=0; isp<num_species; isp++)
  {
    std::ostringstream name;
    name << "/atoms/species_" << isp << "/atomic_number";
    HDFAttribIO<int> h_atomic_number (atomic_numbers[isp]);
    h_atomic_number.read(H5FileID, name.str().c_str());
  }
  app_log()<<"Get Atomic Numbers"<<std::endl;

  esinterface.getAtomicNumbers(atomic_numbers);

  for (int isp=0; isp<num_species; isp++)
    app_log()<<"speciesids: "<<species_ids[isp] << "\t"<< atomic_numbers[isp] <<std::endl;
  IonTypes.resize(species_ids.size());
  for (int i=0; i<species_ids.size(); i++)
    IonTypes[i] = atomic_numbers[species_ids[i]];
  //HDFAttribIO<Vector<TinyVector<double,3> > > h_IonPos(IonPos);
  //h_IonPos.read   (H5FileID, "/atoms/positions");

   for (int i=0; i<IonTypes.size(); i++)
     app_log() << "Atom type(" << i << ") = " << IonTypes[i] << std::endl;
   app_log()<<"get Ion Positions"<<std::endl;
   esinterface.getIonPositions(IonPos);
   app_log()<<"got teh Positions! I think"<<std::endl;

   esinterface.getAtomicOrbitals(AtomicOrbitals);
   std::vector<double> dummy(1);
   esinterface.getTwistData(TwistAngles, dummy, TwistSymmetry);
   std::cerr << TwistAngles[0] << std::endl;

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


//   exit(0);
   return true;
 }

void EinsplineSetBuilder::OccupyBands_Interface(int spin, int sortBands, int numOrbs)
{
  update_token(__FILE__,__LINE__,"OccupyBands_Interface");
  if (myComm->rank() != 0)
    return;

  std::vector<BandInfo>& SortBands(*FullBands[spin]);
  SortBands.clear(); //??? can exit if SortBands is already made?
  int maxOrbs(0);
  app_log()<<"DistincTwists.size()="<<DistinctTwists.size()<<" in OccupyBands_Interface\n";
  app_log()<<" spin="<<spin<<" sortBands="<<sortBands<<" numOrbs="<<numOrbs<<std::endl;
  for (int ti=0; ti<DistinctTwists.size(); ti++)
  {
    app_log()<<"Doop\n";
    app_log()<<"DistinctTwist[ti]="<<ti<<" "<<DistinctTwists[ti]<<std::endl;

    int tindex = DistinctTwists[ti];
//     First, read valence states
    std::ostringstream ePath;
    ePath << "/electrons/kpoint_" << tindex << "/spin_"
          << spin << "/eigenvalues";
    std::vector<double> eigvals;
  H5FileID = H5Fopen(H5FileName.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    std::cerr << H5FileName << std::endl;
    std::cerr << ePath.str() << std::endl;
    HDFAttribIO<std::vector<double> > h_eigvals(eigvals);
    h_eigvals.read(H5FileID, ePath.str().c_str());
//    ESHDF5Interface esinterface; /// ???
    app_log()<<"WEEEEEEEEE"<<std::endl;
//    esinterface.getOrbEigenvals(spin,ti,eigvals);  // HDF5 is complaining here
    app_log()<<"WEEEEEEEEE"<<std::endl;
    for (int bi=0; bi<NumBands; bi++)
    {
      BandInfo band;
      band.IsCoreState = false;
      band.TwistIndex = tindex;
      band.BandIndex  = bi;
      band.MakeTwoCopies = MakeTwoCopies[ti];
    app_log()<<"WEEE"<<bi << std::endl;
      band.Energy = eigvals[bi];  // It's you!!
    app_log()<<"WEEE"<<bi << std::endl;
      if (band.Energy > -1.0e100)
        SortBands.push_back(band);
      if (MakeTwoCopies[ti])
        maxOrbs+=2;
      else
        maxOrbs++;
    }
    // Now, read core states
    for (int cs=0; cs<NumCoreStates; cs++)
    {
      APP_ABORT("Core states not supported with interface yet")
      BandInfo band;
      band.IsCoreState = true;
      band.TwistIndex = tindex;
      band.BandIndex  = cs;
      band.MakeTwoCopies = MakeTwoCopies[ti];
      HDFAttribIO<double> h_energy(band.Energy);
      h_energy.read   (H5FileID, (CoreStatePath(ti,cs)+"eigenvalue").c_str());
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


}

