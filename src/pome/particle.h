#pragma once

#include "basics.h"

namespace pome
{
class CflWavefunction;
struct CflCache;

struct ParticleMove
{
  Integer i_electron;
  Complex z_new;
  Complex z_old;
};

struct ParticleMoveOpt
{
  Real radius;
  explicit ParticleMoveOpt(Real r) : radius(r) { assert(r > 0); }
};

struct ParticleCoordinates
{
  kore::array::Array<Complex, 1> coordinates;

  explicit ParticleCoordinates(Integer n_electron) : coordinates(n_electron) { }
  Complex operator()(Integer i_electron) const { return coordinates(i_electron); }
  Complex& operator()(Integer i_electron) { return coordinates(i_electron); }
};

} // namespace pome