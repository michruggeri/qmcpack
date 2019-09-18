//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/Crowd.h"

namespace qmcplusplus
{
void Crowd::clearResults()
{
  // These were cleared to 1.0 each loop by VMCUpdatePbyP advance walker
  // refactored code may depend on this initial value.
  std::fill(log_gf_.begin(), log_gf_.end(), 1.0);
  std::fill(log_gb_.begin(), log_gb_.end(), 1.0);
}


void Crowd::reserve(int crowd_size)
{
  auto reserveCS = [crowd_size](auto& avector) { avector.reserve(crowd_size); };
  reserveCS(mcp_walkers_);
  reserveCS(walker_elecs_);
  reserveCS(walker_twfs_);
  reserveCS(walker_hamiltonians_);

  auto resizeCS = [crowd_size](auto& avector) { avector.resize(crowd_size); };
  resizeCS(grads_now_);
  resizeCS(grads_new_);
  resizeCS(ratios_);
  resizeCS(log_gf_);
  resizeCS(log_gb_);
  resizeCS(prob_);
}

void Crowd::addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian)
{
  mcp_walkers_.push_back(walker);
  mcp_wfbuffers_.push_back(walker.DataSet);
  walker_elecs_.push_back(elecs);
  walker_twfs_.push_back(twf);
  walker_hamiltonians_.push_back(hamiltonian);
};

void Crowd::loadWalkers()
{
  auto it_walker       = mcp_walkers_.begin();
  auto it_walker_elecs = walker_elecs_.begin();
  //flex walkers here
  while (it_walker != mcp_walkers_.end())
  {
    (*it_walker_elecs).get().loadWalker(*it_walker, true);
    ++it_walker;
    ++it_walker_elecs;
  }
}

void Crowd::startBlock(int num_steps)
{
  n_accept = 0;
  n_reject = 0;
  estimator_manager_crowd_.startBlock(num_steps);
}

} // namespace qmcplusplus
