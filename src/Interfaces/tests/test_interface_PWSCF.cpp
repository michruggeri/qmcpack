//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSetBuilderInterface.h"


#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Einspline SPO from PWSCF interface", "[interfaces]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(1);
  ions_.R[0][0] = 9.44863;
  ions_.R[0][1] = 9.44863;
  ions_.R[0][2] = 9.44863;


  elec_.setName("elec");
  elec_.create(2);
  elec_.R[0][0] = 0.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  // monoO
  /*
  elec_.Lattice.R(0,0) = 5.10509515;
  elec_.Lattice.R(0,1) = -3.23993545;
  elec_.Lattice.R(0,2) = 0.0;
  elec_.Lattice.R(1,0) = 5.10509515;
  elec_.Lattice.R(1,1) = 3.23993545;
  elec_.Lattice.R(1,2) = 0.0;
  elec_.Lattice.R(2,0) = -6.49690625;
  elec_.Lattice.R(2,1) = 0.0;
  elec_.Lattice.R(2,2) = 7.08268015;
 */

  // diamondC_1x1x1
  elec_.Lattice.R(0, 0) = 1.889726e+01;
  elec_.Lattice.R(0, 1) = 0.0;
  elec_.Lattice.R(0, 2) = 0.0;
  elec_.Lattice.R(1, 0) = 0.0;
  elec_.Lattice.R(1, 1) = 1.889726e+01;
  elec_.Lattice.R(1, 2) = 0.0;
  elec_.Lattice.R(2, 0) = 0.0;
  elec_.Lattice.R(2, 1) = 0.0;
  elec_.Lattice.R(2, 2) = 1.889726e+01;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
//  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
//  tspecies(chargeIdx, downIdx) = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_, DT_SOA);
#else
  elec_.addTable(ions_, DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();


  TrialWaveFunction psi(c);
  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"qmcqepack\" href=\"pwscf.in\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";


  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
//  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilderInterface einSet(elec_, ptcl.getPool(), c, ein1,"PWSCF");
  SPOSet* spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != NULL);

#if defined(QMC_COMPLEX)
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // for vgl
  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // value
  REQUIRE(psiM[1][0] == ComplexApprox(std::complex<double>(-0.00400826,-0.00165963)));
  REQUIRE(psiM[1][1] == ComplexApprox(std::complex<double>( 0.00570092,-0.084762)));
  // grad
  REQUIRE(dpsiM[1][0][0] == ComplexApprox(std::complex<double> (-0.04942   , -0.00672706)));
  REQUIRE(dpsiM[1][0][1] == ComplexApprox(std::complex<double> ( 0.0524425 ,  0.00218184)));
  REQUIRE(dpsiM[1][0][2] == ComplexApprox(std::complex<double> (-0.0494052 , -0.00673416)));
  REQUIRE(dpsiM[1][1][0] == ComplexApprox(std::complex<double> ( 0.0157443 , -0.00742662)));
  REQUIRE(dpsiM[1][1][1] == ComplexApprox(std::complex<double> (-0.0683265 ,  0.00343797)));
  REQUIRE(dpsiM[1][1][2] == ComplexApprox(std::complex<double> ( 0.00937362, -0.00479787)));
  // lapl
  REQUIRE(d2psiM[1][0] == ComplexApprox(std::complex<double> ( 0.0514297,0.0306235)));
  REQUIRE(d2psiM[1][1] == ComplexApprox(std::complex<double> (-0.0325769,0.948198 )));
  // for vgh
//  SPOSet::ValueVector_t psiV(psiM[1], spo->getOrbitalSetSize());
//  SPOSet::GradVector_t dpsiV(dpsiM[1], spo->getOrbitalSetSize());
//  SPOSet::HessVector_t ddpsiV(spo->getOrbitalSetSize());
//  spo->evaluate(elec_, 1, psiV, dpsiV, ddpsiV);

// for vgh
  SPOSet::ValueVector_t psiV(psiM[1], spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsiV(dpsiM[1], spo->getOrbitalSetSize());
  SPOSet::HessVector_t ddpsiV(spo->getOrbitalSetSize());
  spo->evaluateVGH(elec_, 1, psiV, dpsiV, ddpsiV);

  // Catch default is 100*(float epsilson)
  double eps = 2000 * std::numeric_limits<float>::epsilon();
  //hess
  REQUIRE(ddpsiV[1](0, 0) == ComplexApprox(std::complex<double>(-0.011268,0.300841)));
  REQUIRE(ddpsiV[1](0, 1) == ComplexApprox(std::complex<double>(0.0523478,-0.0266171)));
  REQUIRE(ddpsiV[1](0, 2) == ComplexApprox(std::complex<double>(0.0654908,0.123955)));
  REQUIRE(ddpsiV[1](1, 0) == ComplexApprox(std::complex<double>(0.0532509,-0.0196289)));
  REQUIRE(ddpsiV[1](1, 1) == ComplexApprox(std::complex<double>(-0.00465563,0.326825)).epsilon(eps));
  REQUIRE(ddpsiV[1](1, 2) == ComplexApprox(std::complex<double>(0.0577793,-0.0239127)));
  REQUIRE(ddpsiV[1](2, 0) == ComplexApprox(std::complex<double>(0.0657093,0.124484)));
  REQUIRE(ddpsiV[1](2, 1) == ComplexApprox(std::complex<double>(0.0570947,-0.0303713)));
  REQUIRE(ddpsiV[1](2, 2) == ComplexApprox(std::complex<double>(-0.0166532,0.320532)));

  SPOSet::HessMatrix_t hesspsiV(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GGGMatrix_t d3psiV(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, hesspsiV, d3psiV);

#endif
}
} // namespace qmcplusplus
