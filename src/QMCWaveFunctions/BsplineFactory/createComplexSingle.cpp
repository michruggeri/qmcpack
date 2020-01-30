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
#include "Numerics/e2iphi.h"
#include "simd/vmath.hpp"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSetBuilderInterface.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2R.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2C.h"
#if defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/BsplineFactory/SplineC2ROMP.h"
#endif
#include "QMCWaveFunctions/BsplineFactory/HybridRepCplx.h"
#include <fftw3.h>
#include <QMCWaveFunctions/einspline_helper.hpp>
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
//<<<<<<< HEAD
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderInterface.h"
//#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"
//#include "QMCWaveFunctions/BsplineFactory/SplineHybridAdoptorReaderP.h"
//=======
//#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderPInterface.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReader.h"
#include "QMCWaveFunctions/BsplineFactory/SplineSetReaderInterface.h"
#include "QMCWaveFunctions/BsplineFactory/HybridRepSetReader.h"
//>>>>>>> 7d0e44ea47f573c19d9f2c33afe1ad7b896c8d55

namespace qmcplusplus
{
BsplineReaderBase* createBsplineComplexSingle(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU)
{
  typedef OHMMS_PRECISION RealType;
  BsplineReaderBase* aReader = nullptr;

#if defined(QMC_COMPLEX)
  if (hybrid_rep)
    aReader = nullptr;//new HybridRepSetReader<HybridRepCplx<SplineC2C<float>>>(e);
  else
    aReader = new SplineSetReader<SplineC2C<float>>(e);
#else //QMC_COMPLEX
#if defined(ENABLE_OFFLOAD)
  if (useGPU == "yes")
  {
    if (hybrid_rep)
    {
      APP_ABORT("OpenMP offload has not been enabled with hybrid orbital representation!");
    }
    else
      aReader = new SplineSetReader<SplineC2ROMP<float>>(e);
  }
  else
#endif
  {
    if (hybrid_rep)
      aReader = new HybridRepSetReader<HybridRepCplx<SplineC2R<float>>>(e);
    else
      aReader = new SplineSetReader<SplineC2R<float>>(e);
  }
#endif
  return aReader;
}

BsplineReaderInterface* createBsplineComplexSingle(EinsplineSetBuilderInterface* e, bool hybrid_rep, const std::string& useGPU)
{
  typedef OHMMS_PRECISION RealType;
  BsplineReaderInterface* aReader = nullptr;

#if defined(QMC_COMPLEX)
    aReader = new SplineSetReaderInterface<SplineC2C<float>>(e);
#else //QMC_COMPLEX
    aReader = new SplineSetReaderInterface<SplineC2R<float>>(e);
#endif
  return aReader;
}
} // namespace qmcplusplus
