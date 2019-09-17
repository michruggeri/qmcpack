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
TEST_CASE("Einspline SPO from PWSCF interface", "[wavefunction]")
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
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

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
<determinantset type=\"qmcqepack\" href=\"pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

  // monoO
  //<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"6\"/> \

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilderInterface einSet(elec_, ptcl.getPool(), c, ein1,"PWSCF");
  SPOSet* spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != NULL);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
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
  SPOSet::ValueVector_t psiV(psiM[1], spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsiV(dpsiM[1], spo->getOrbitalSetSize());
  SPOSet::HessVector_t ddpsiV(spo->getOrbitalSetSize());
  spo->evaluate(elec_, 1, psiV, dpsiV, ddpsiV);

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
/*
  REQUIRE(ddpsiV[1](0, 1) == ComplexApprox(1.8089479397).compare_real_only());
  REQUIRE(ddpsiV[1](0, 2) == ComplexApprox(0.5608575749).compare_real_only());
  REQUIRE(ddpsiV[1](1, 0) == ComplexApprox(1.8089479397).compare_real_only());
  REQUIRE(ddpsiV[1](1, 1) == ComplexApprox(-0.07996207476).epsilon(eps).compare_real_only());
  REQUIRE(ddpsiV[1](1, 2) == ComplexApprox(0.5237969314).compare_real_only());
  REQUIRE(ddpsiV[1](2, 0) == ComplexApprox(0.5608575749).compare_real_only());
  REQUIRE(ddpsiV[1](2, 1) == ComplexApprox(0.5237969314).compare_real_only());
  REQUIRE(ddpsiV[1](2, 2) == ComplexApprox(-2.316497764).compare_real_only());
*/
  SPOSet::HessMatrix_t hesspsiV(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GGGMatrix_t d3psiV(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, hesspsiV, d3psiV);

  //The reference values for grad_grad_grad_psi.
  /*
  d3psiV(1,0)[0][0]=(0.046337127685546875,-0.046337127685546875)
  d3psiV(1,0)[0][1]=(1.1755813360214233,-1.1755813360214233)
  d3psiV(1,0)[0][2]=(0.066015571355819702,-0.066015541553497314)
  d3psiV(1,0)[0][4]=(0.041470438241958618,-0.041470438241958618)
  d3psiV(1,0)[0][5]=(-0.51674127578735352,0.51674127578735352)
  d3psiV(1,0)[0][8]=(0.065953642129898071,-0.065953642129898071)
  d3psiV(1,0)[1][4]=(-4.8771157264709473,4.8771157264709473)
  d3psiV(1,0)[1][5]=(0.041532635688781738,-0.041532605886459351)
  d3psiV(1,0)[1][8]=(1.1755810976028442,-1.1755810976028442)
  d3psiV(1,0)[2][8]=(0.046399354934692383,-0.046399354934692383)
 
  d3psiV(1,1)[0][0]=(6.7155771255493164,-7.5906991958618164)
  d3psiV(1,1)[0][1]=(5.545051097869873,-5.0280308723449707)
  d3psiV(1,1)[0][2]=(0.98297119140625,-0.50021600723266602)
  d3psiV(1,1)[0][4]=(-3.1704092025756836,3.8900821208953857)
  d3psiV(1,1)[0][5]=(-1.9537661075592041,1.7758266925811768)
  d3psiV(1,1)[0][8]=(1.9305641651153564,-2.1480715274810791)
  d3psiV(1,1)[1][4]=(3.605137825012207,-3.2767453193664551)
  d3psiV(1,1)[1][5]=(-0.73825764656066895,-0.33745908737182617)
  d3psiV(1,1)[1][8]=(5.5741839408874512,-5.0784988403320312)
  d3psiV(1,1)[2][8]=(3.131234884262085,-1.3596141338348389)
*/

#if 0 //Enable when finite precision issue on Rhea is found.
  REQUIRE(d3psiV(1,0)[0][0] ==ComplexApprox(0.0463371276).compare_real_only());
  REQUIRE(d3psiV(1,0)[0][1] ==ComplexApprox(1.1755813360).compare_real_only());
  REQUIRE(d3psiV(1,0)[0][2] ==ComplexApprox(0.0660155713).compare_real_only());
  REQUIRE(d3psiV(1,0)[0][4] ==ComplexApprox(0.0414704382).compare_real_only());
  REQUIRE(d3psiV(1,0)[0][5] ==ComplexApprox(-0.5167412758).compare_real_only());
  REQUIRE(d3psiV(1,0)[0][8] ==ComplexApprox(0.0659536421).compare_real_only());
  REQUIRE(d3psiV(1,0)[1][4] ==ComplexApprox(-4.8771157264).compare_real_only());
  REQUIRE(d3psiV(1,0)[1][5] ==ComplexApprox(0.0415326356).compare_real_only());
  REQUIRE(d3psiV(1,0)[1][8] ==ComplexApprox(1.1755810976).compare_real_only());
  REQUIRE(d3psiV(1,0)[2][8] ==ComplexApprox(0.0463993549).compare_real_only());
  
  REQUIRE(d3psiV(1,1)[0][0] ==ComplexApprox(6.7155771255).compare_real_only());
  REQUIRE(d3psiV(1,1)[0][1] ==ComplexApprox(5.5450510978).compare_real_only());
  REQUIRE(d3psiV(1,1)[0][2] ==ComplexApprox(0.9829711914).compare_real_only());
  REQUIRE(d3psiV(1,1)[0][4] ==ComplexApprox(-3.1704092025).compare_real_only());
  REQUIRE(d3psiV(1,1)[0][5] ==ComplexApprox(-1.9537661075).compare_real_only());
  REQUIRE(d3psiV(1,1)[0][8] ==ComplexApprox(1.9305641651).compare_real_only());
  REQUIRE(d3psiV(1,1)[1][4] ==ComplexApprox(3.6051378250).compare_real_only());
  REQUIRE(d3psiV(1,1)[1][5] ==ComplexApprox(-0.7382576465).compare_real_only());
  REQUIRE(d3psiV(1,1)[1][8] ==ComplexApprox(5.5741839408).compare_real_only());
  REQUIRE(d3psiV(1,1)[2][8] ==ComplexApprox(3.1312348842).compare_real_only());
#endif

#endif

#if 0
  // Dump values of the orbitals
  int orbSize= spo->getOrbitalSetSize();
  int basisSize= spo->getBasisSetSize();
  printf("orb size = %d basis set size = %d\n",orbSize, basisSize);

  FILE *fspo = fopen("spo.dat", "w");
  for (int ix = 0; ix < 30; ix++) {
    for (int iy = 0; iy < 30; iy++) {
      for (int iz = 0; iz < 30; iz++) {
        double x = 0.1*ix - 1.5;
        double y = 0.1*iy - 1.5;
        double z = 0.1*iz - 1.5;
        elec_.R[0][0] = x;
        elec_.R[0][1] = y;
        elec_.R[0][2] = z;
        elec_.update();
        SPOSet::ValueVector_t orbs(orbSize);
        spo->evaluate(elec_, 0, orbs);
        fprintf(fspo, "%g %g %g",x,y,z);
        for (int j = 0; j < orbSize; j++) {
          fprintf(fspo, " %g ",orbs[j]);
        }
        fprintf(fspo, "\n");
      }
    }
  }
  fclose(fspo);
#endif
}
} // namespace qmcplusplus
