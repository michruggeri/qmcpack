//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSetBuilderInterface.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2R.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepReal.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
//<<<<<<< HEAD
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderInterface.h"
//#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"
//#include "QMCWaveFunctions/BsplineFactory/SplineHybridAdoptorReaderP.h"
//=======
//#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderPInterface.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReaderInterface.h"
//#include "QMCWaveFunctions/BsplineFactory/HybridRepSetReader.h"
//>>>>>>> 7d0e44ea47f573c19d9f2c33afe1ad7b896c8d55

namespace qmcplusplus
{
BsplineReaderBase* createBsplineRealSingle(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU)
{
  BsplineReaderBase* aReader = nullptr;

  if (hybrid_rep)
    aReader = nullptr;// new HybridRepSetReader<HybridRepReal<SplineR2R<float>>>(e);
  else
    aReader = new SplineSetReader<SplineR2R<float>>(e);
  return aReader;
}
BsplineReaderInterface* createBsplineRealSingle(EinsplineSetBuilderInterface* e, bool hybrid_rep, const std::string& useGPU)
{
  BsplineReaderInterface* aReader = nullptr;

//  aReader = new SplineSetReaderInterface<SplineR2RSoA<float, OHMMS_PRECISION>>(e);
  //aReader = new SplineAdoptorReaderInterface<SplineR2RSoA<float, OHMMS_PRECISION>>(e);
  return aReader;
}
} // namespace qmcplusplus
