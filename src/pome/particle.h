#pragma once

#include <random>
#include "basics.h"

namespace pome
{
class CflWavefunction;
struct CflCache;

struct ParticleMove
{
  Integer i_electron;
  Complex z_new;
};

struct ParticleMoveOpt
{
  Real radius;
  ParticleMoveOpt(Real r) : radius(r) { assert(r > 0); }
};

struct ParticleCoordinates
{
  kore::array::Array<Complex, 1> coordinates;

  ParticleCoordinates(Integer n_electron) : coordinates(n_electron) { }
  Complex operator()(Integer i_electron) const { return coordinates(i_electron); }
  Complex& operator()(Integer i_electron) { return coordinates(i_electron); }
};

} // namespace pome