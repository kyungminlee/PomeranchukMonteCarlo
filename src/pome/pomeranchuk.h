#pragma once

#include "geometry.h"
#include "wavefunction.h"
#include "measurement.h"

namespace pome {
class PomeranchukProblem
{
public:
  using SystemType = CflWavefunction;
  using StateType = ParticleCoordinates;
  using UpdateType = ParticleMove;
  using UpdateOptType = ParticleMoveOpt;
  using CacheType = CflCache;
  using MeasurementType = Measurement;
  using RealType = double;
};
}