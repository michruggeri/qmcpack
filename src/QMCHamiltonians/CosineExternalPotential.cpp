//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <QMCHamiltonians/CosineExternalPotential.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{
bool CosineExternalPotential::put(xmlNodePtr cur)
{
  using std::sqrt;
  intensity = -1.0;
  qx = -1.0;
  qy = -1.0;
  qz = -1.0;

  OhmmsAttributeSet attrib;
  attrib.add(qx, "qx");
  attrib.add(qy, "qy");
  attrib.add(qz, "qz");
  attrib.add(intensity, "intensity");
  attrib.put(cur);

  return true;
}


bool CosineExternalPotential::get(std::ostream& os) const
{
  os << "External Cosine potential" << std::endl;
  return true;
}


OperatorBase* CosineExternalPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return new CosineExternalPotential(*this);
}


CosineExternalPotential::Return_t CosineExternalPotential::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles)
    Value = evaluate_sp(P);
  else
  {
#endif
    Value              = 0.0;
    for (int i = 0; i < P.getTotalNum(); ++i)
    {
      PosType r   = P.R[i];
      P.Lattice.applyMinimumImage(r);
      double rr = r[0]*qx + r[1]*qy + r[2]*qz;
      RealType v1 = intensity*std::cos(rr);
      Value += v1;
    }
#if !defined(REMOVE_TRACEMANAGER)
  }
#endif
  return Value;
}


#if !defined(REMOVE_TRACEMANAGER)
CosineExternalPotential::Return_t CosineExternalPotential::evaluate_sp(ParticleSet& P)
{
  Array<TraceReal, 1>& V_samp = *V_sample;
  Value                       = 0.0;
  for (int i = 0; i < P.getTotalNum(); ++i)
  {
    PosType r   = P.R[i];
    P.Lattice.applyMinimumImage(r);
    double rr = r[0]*qx + r[1]*qy + r[2]*qz;
    RealType v1 = intensity*std::cos(rr);
    V_samp(i)   = v1;
    Value += v1;
  }
  return Value;
}
#endif

} // namespace qmcplusplus
