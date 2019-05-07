#ifndef QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H

#include <Message/MPIObjectBase.h>
#include "AFQMC/config.h"
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"
#include "Utilities/Timer.h"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class BackPropagatedEstimator: public EstimatorBase
{

  // allocators
  using Allocator = device_allocator<ComplexType>;

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;
  using stdCVector = boost::multi::array<ComplexType,1>;
  using stdCMatrix = boost::multi::array<ComplexType,2>;
  using stdCTensor = boost::multi::array<ComplexType,3>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;

  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, WALKER_TYPES wlk, WalkerSet& wset, 
        Wavefunction& wfn, Propagator& prop, bool impsamp_=true) :
                                          EstimatorBase(info),TG(tg_), walker_type(wlk),
                                          Refs({0,0,0},shared_allocator<ComplexType>{TG.TG_local()}),
                                          wfn0(wfn), prop0(prop),
                                          greens_function(false),max_nback_prop(10),
                                          nStabalize(10), path_restoration(false), block_size(1),
                                          writer(false), importanceSampling(impsamp_)
  {

    int nave(1);
    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(max_nback_prop, "nsteps", "int");
      m_param.add(nave, "naverages", "int");
      m_param.add(restore_paths, "path_restoration", "std::string");
      m_param.add(block_size, "block_size", "int");
      m_param.add(nblocks_skip, "nskip", "int");
      m_param.put(cur);
      if(restore_paths == "true") {
        path_restoration = true;
      } else {
        path_restoration = false;
      }
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    nback_prop_steps.reserve(nave);
    for(int i=1; i<nave; i++)
      nback_prop_steps.push_back(i*max_nback_prop/nave);
    nback_prop_steps.push_back(max_nback_prop);

    // sort the requested number of steps
    std::sort(nback_prop_steps.begin(),nback_prop_steps.end());

    if(max_nback_prop <= 0)
      APP_ABORT("max_nback_prop <= 0 is not allowed.\n");

    int ncv(prop0.global_number_of_cholesky_vectors());
    int nref(wfn0.number_of_references_for_back_propagation());
    wset.resize_bp(max_nback_prop,ncv,nref);
    wset.setBPPos(0);
    // set SMN in case BP begins right away
    if(nblocks_skip==0)
      for(auto it=wset.begin(); it<wset.end(); ++it)
        it->setSlaterMatrixN(); 

    ncores_per_TG = TG.getNCoresPerTG();
    if(ncores_per_TG > 1)
      APP_ABORT("ncores > 1 is broken with back propagation. Fix this.");
    core_rank = TG.getLocalTGRank();
    writer = (TG.getGlobalRank()==0);
    if(walker_type == CLOSED) {
      dm_size = NMO*NMO;
      dm_dims = {NMO,NMO};
    } else if(walker_type == COLLINEAR) {
      dm_size = 2*NMO*NMO;
      dm_dims = {2*NMO,NMO};
    } else if(walker_type == NONCOLLINEAR) {
      dm_size = 4*NMO*NMO;
      dm_dims = {2*NMO,2*NMO};
    }
    if(DMBuffer.size(0) != nave || DMBuffer.size(1) != dm_size) 
      DMBuffer.reextent({nave,dm_size});
    if(DMAverage.size(0) != nave || DMAverage.size(1) != dm_size) 
      DMAverage.reextent({nave,dm_size});
    using std::fill_n;
    fill_n(DMBuffer.origin(), DMBuffer.num_elements(), ComplexType(0.0,0.0));
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
    denom.reextent({nave,1});
    denom_average.reextent({nave,1});
    fill_n(denom_average.origin(), nave, ComplexType(0.0)); 
  }

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    accumulated_in_last_block=false;
    int bp_step = wset.getBPPos();
    if(bp_step <=0)
      APP_ABORT(" Error: Found bp_step <=0 in BackPropagate::accumulate_block. \n"); 
    if(bp_step > max_nback_prop)
      APP_ABORT(" Error: max_nback_prop in back propagation estimator must be conmensuate with nStep*nSubStep.\n");
    if(max_nback_prop > wset.NumBackProp())
      APP_ABORT(" Error: max_nback_prop > wset.NumBackProp() \n");

    // check if measurement is needed 
    int iav(-1);
    for(int i=0; i<nback_prop_steps.size(); i++) {
      if(bp_step == nback_prop_steps[i]) {
        iav = i; 
        break;
      }    
    }
    if(iav < 0) return;

    using std::fill_n;
    // 0. skip if requested 
    if(bp_step == max_nback_prop && iblock < nblocks_skip) {
      if( iblock+1 == nblocks_skip )
        for(auto it=wset.begin(); it<wset.end(); ++it)
          it->setSlaterMatrixN(); 
      iblock++;
      wset.setBPPos(0);
      return;
    }

    AFQMCTimers[back_propagate_timer]->start();
    int nrefs = wfn0.number_of_references_for_back_propagation();
    int nrow(NMO*((walker_type==NONCOLLINEAR)?2:1));
    int ncol(NAEA+((walker_type==CLOSED)?0:NAEB));
    int nx((walker_type==COLLINEAR)?2:1);

    // 1. check structures
    if(Refs.size(0) != wset.size() || Refs.size(1) != nrefs || Refs.size(2) !=nrow*ncol) 
      Refs = std::move(mpi3CTensor({wset.size(),nrefs,nrow*ncol},Refs.get_allocator()));
    if(detR.size(0) != wset.size() || detR.size(1) != nx*nrefs)
      detR.reextent({wset.size(),nrefs*nx}); 

    // temporary!
    int wpop = wset.GlobalPopulation();
    int nw = wset.size();  // assuming all groups have the same size
    int iw0 = wset.getTG().getTGNumber()*nw;
    if(wdetR.size(0) != nback_prop_steps.size() ||  wdetR.size(1) !=  wpop || wdetR.size(2) != nx*nrefs)
      wdetR.reextent({nback_prop_steps.size(),wpop,nrefs*nx}); 
    if(wDMsum.size(0) != nback_prop_steps.size() ||  wDMsum.size(1) !=  wpop || wDMsum.size(2) != nx*nrefs)
      wDMsum.reextent({nback_prop_steps.size(),wpop,nrefs*nx}); 
    if(wOvlp.size(0) != nback_prop_steps.size() ||  wOvlp.size(1) !=  wpop || wOvlp.size(2) != nx*nrefs)
      wOvlp.reextent({nback_prop_steps.size(),wpop,nrefs*nx}); 

    int n0,n1;
    std::tie(n0,n1) = FairDivideBoundary(TG.getLocalTGRank(),int(Refs.size(2)),TG.getNCoresPerTG());
    boost::multi::array_ref<ComplexType,3> Refs_(to_address(Refs.origin()),Refs.extensions());
    fill_n(denom[iav].origin(),denom[iav].num_elements(),ComplexType(0.0,0.0));
    fill_n(DMBuffer[iav].origin(), DMBuffer[iav].num_elements(), ComplexType(0.0,0.0));

    // 2. setup back propagated references
    wfn0.getReferencesForBackPropagation(Refs_[0]);
    for(int iw=1; iw<wset.size(); ++iw)
      for(int ref=0; ref<nrefs; ++ref)
       copy_n(Refs_[0][ref].origin()+n0,n1-n0,Refs_[iw][ref].origin()+n0);
    TG.TG_local().barrier();

    //3. propagate backwards the references
    prop0.BackPropagate(bp_step,nStabalize,wset,Refs_,detR);

    // MAM:
    using std::copy_n;
    copy_n(detR.origin(),detR.num_elements(),wdetR[iav][iw0].origin());

    //4. calculate properties, only rdm now but make a list or properties later
    // adjust weights here is path restoration
    stdCVector wgt(iextensions<1u>{wset.size()});
    wset.getProperty(WEIGHT,wgt);
    if(path_restoration) {
      auto&& factors(wset.getWeightFactors());
      for(int k=0; k<bp_step; k++)  
        for(int i=0; i<wgt.size(); i++) 
          wgt[i] *= factors[k][i];
    } else if(!importanceSampling) {
      stdCVector phase(iextensions<1u>{wset.size()});
      wset.getProperty(PHASE,phase);
      for(int i=0; i<wgt.size(); i++) wgt[i] *= phase[i];
    }
    CMatrix_ref BackPropDM(DMBuffer[iav].origin(), {dm_dims.first,dm_dims.second});
    stdCVector_ref denom_(denom[iav].origin(), iextensions<1u>{denom.size(1)});
    
    wfn0.WalkerAveragedDensityMatrix(wset, wgt, BackPropDM, denom_, wOvlp[iav].sliced(iw0,iw0+nw), wDMsum[iav].sliced(iw0,iw0+nw), 
                                   !importanceSampling, std::addressof(Refs_),
                                   std::addressof(detR));

    if(bp_step == max_nback_prop) {
      // 5. setup for next block 
      for(auto it=wset.begin(); it<wset.end(); ++it)
        it->setSlaterMatrixN(); 
      wset.setBPPos(0);

      // 6. increase block counter
      iblock++;
      accumulated_in_last_block=true;
    }
    AFQMCTimers[back_propagate_timer]->stop();
  }

  void tags(std::ofstream& out) {
    if(writer) 
      out<<"BP_timer ";
  }

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset)
  {
    // I doubt we will ever collect a billion blocks of data.
    int n_zero = 9;
    if(writer) {
      out<<std::setprecision(5)
                    <<AFQMCTimers[back_propagate_timer]->get_total() <<" ";
      AFQMCTimers[back_propagate_timer]->reset();  
    }
    if(accumulated_in_last_block) {
      TG.Global().reduce_in_place_n(wDMsum.origin(),wDMsum.num_elements(),std::plus<>());
      TG.Global().reduce_in_place_n(wOvlp.origin(),wOvlp.num_elements(),std::plus<>());
      // not correct with ncores>1  
      TG.Global().reduce_in_place_n(wdetR.origin(),wdetR.num_elements(),std::plus<>());
    }
    if(writer && accumulated_in_last_block) {
      int nave(nback_prop_steps.size());
      int nref(detR.size(1));
      if(write_metadata) {
        dump.push("Metadata");
        dump.write(max_nback_prop, "NumBackProp");
        dump.write(nave, "NumAverages");
        dump.write(nref, "NumReferences");
        dump.pop();
        write_metadata = false;
      }
#ifdef ENABLE_CUDA
      stdCMatrix buff(DMBuffer);
#else
      CMatrix& buff(DMBuffer);
#endif
      ma::axpy(ComplexType(1.0),buff,DMAverage);
      for(int i=0; i<nave; ++i)  
        denom_average[i][0] += denom[i][0];
      if(iblock%block_size == 0) {
//        for(int i = 0; i < DMAverage.size(); i++)
//          DMAverage[i] /= block_size;
//        denom_average[0] /= block_size;
        ma::scal(ComplexType(1.0/block_size),DMAverage);
        ma::scal(ComplexType(1.0/block_size),denom_average);
        dump.push("BackPropagated");
        for(int i=0; i<nave; ++i) {
          dump.push(std::string("NumBackProp_")+std::to_string(nback_prop_steps[i]));
          std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          stdCVector_ref DMAverage_( DMAverage[i].origin(), {DMAverage.size(1)}); 
          stdCVector_ref denom_average_( denom_average[i].origin(), {denom_average.size(1)}); 
          stdCVector_ref wOvlp_( wOvlp[i].origin(), {wOvlp.size(1)*wOvlp.size(2)}); 
          stdCVector_ref wDMsum_( wDMsum[i].origin(), {wDMsum.size(1)*wDMsum.size(2)}); 
          stdCVector_ref wdetR_( wdetR[i].origin(), {wdetR.size(1)*wdetR.size(2)}); 
          dump.write(DMAverage_, "one_rdm_"+padded_iblock);
          dump.write(denom_average_, "one_rdm_denom_"+padded_iblock);
          dump.write(wOvlp_, "one_rdm_walker_overlaps_"+padded_iblock);
          dump.write(wDMsum_, "one_rdm_walker_dm_sums_"+padded_iblock);
          dump.write(wdetR_, "one_rdm_detR_"+padded_iblock);
          dump.pop();
        }
        dump.pop();
        using std::fill_n;  
        fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
        fill_n(denom_average.origin(), denom_average.num_elements(), ComplexType(0.0,0.0));
        fill_n(wdetR.origin(),wdetR.num_elements(),ComplexType(0.0));
        fill_n(wOvlp.origin(),wOvlp.num_elements(),ComplexType(0.0));
        fill_n(wDMsum.origin(),wDMsum.num_elements(),ComplexType(0.0));
      }
    }
  }

  private:

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  bool writer;
  bool accumulated_in_last_block;

  int max_nback_prop;
  std::vector<int> nback_prop_steps;  

  Wavefunction& wfn0;

  Propagator& prop0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  CMatrix DMBuffer;
  stdCMatrix DMAverage;
  mpi3CTensor Refs;
  boost::multi::array<ComplexType,2> detR;
  boost::multi::array<ComplexType,3> wdetR;
  boost::multi::array<ComplexType,3> wDMsum;
  boost::multi::array<ComplexType,3> wOvlp;

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;
  int iblock = 0;
  int nblocks_skip = 0;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  // Print the full mixed estimator for the one-particle reduced density matrix.
  bool greens_function;
  // Frequency of reorthogonalisation.
  int nStabalize;
  // Block size over which RDM will be averaged.
  int block_size;
  // Whether to restore cosine projection and real local energy apprximation for weights
  // along back propagation path.
  bool path_restoration, importanceSampling;
  std::vector<ComplexType> weights;
  int dm_size;
  std::pair<int,int> dm_dims;
  stdCMatrix denom;
  stdCMatrix denom_average;
  bool write_metadata = true;

};
}
}

#endif
