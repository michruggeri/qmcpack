//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/VMC/VMCFactoryNew.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "Message/OpenMP.h"

#ifdef QMC_CUDA
#include "QMCDrivers/VMC/VMC_CUDA.h"
#endif

namespace qmcplusplus
{
QMCDriverInterface* VMCFactoryNew::create(MCPopulation& pop,
                                       TrialWaveFunction& psi,
                                       QMCHamiltonian& h,
                                       ParticleSetPool& ptclpool,
                                       HamiltonianPool& hpool,
                                       WaveFunctionPool& ppool,
                                       Communicate* comm)
{
  int np = omp_get_max_threads();
  //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
  QMCDriverInterface* qmc = nullptr;
  if (VMCMode == 0 || VMCMode == 1) //(0,0,0) (0,0,1)
  {
    qmc = new VMCBatched(pop, psi, h, ppool, comm);
  }
  else
  {
    APP_ABORT("VMCBatch driver not yet supported");
  }
  qmc->setUpdateMode(VMCMode & 1);
  return qmc;
}
} // namespace qmcplusplus
