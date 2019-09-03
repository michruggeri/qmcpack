//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VMCBATCHED_H
#define QMCPLUSPLUS_VMCBATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "Particle/MCPopulation.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCBatched : public QMCDriverNew
{
public:
  /// Constructor.
  VMCBatched(QMCDriverInput&& qmcdriver_input,
             VMCDriverInput&& input,
             MCPopulation& pop,
             TrialWaveFunction& psi,
             QMCHamiltonian& h,
             WaveFunctionPool& ppool,
             Communicate* comm);

  bool run();
  void setup();
  //inline std::vector<RandomGenerator_t*>& getRng() { return Rng;}
  IndexType calc_default_local_walkers();

private:
  int prevSteps;
  int prevStepsBetweenSamples;
  VMCDriverInput vmcdriver_input_;
  QMCRunType getRunType() { return QMCRunType::VMC_BATCH; }
  ///Ways to set rn constant
  RealType logoffset, logepsilon;
  ///option to enable/disable drift equation or RN for VMC
  std::string UseDrift;
  ///check the run-time environments
  void resetRun();
  ///copy constructor
  VMCBatched(const VMCBatched&) = delete;
  /// Copy operator (disabled).
  VMCBatched& operator=(const VMCBatched&) = delete;
};

extern std::ostream& operator<<(std::ostream& o_stream, const VMCBatched& vmc_batched);
} // namespace qmcplusplus

#endif
