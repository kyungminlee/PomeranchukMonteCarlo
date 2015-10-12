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

enum class ParticleMoveStrategy { UNSET, CIRCLE, SQUARE };

struct ParticleMoveOpt
{
  Real radius;
  ParticleMoveStrategy strategy;
};

inline
std::ostream& operator<<(std::ostream& os, ParticleMoveOpt opt)
{
  os << "ParticleMoveOpt[";
  os << "radius=" << opt.radius;
  os << ",strategy=";
  switch(opt.strategy) {
   case ParticleMoveStrategy::UNSET:  os << "unset"; break;
   case ParticleMoveStrategy::CIRCLE: os << "circle"; break;
   case ParticleMoveStrategy::SQUARE: os << "square"; break;
   default: os << "error";
  }
  os << "]";
  return os;
}

struct ParticleCoordinates
{
  kore::array::Array<Complex, 1> coordinates;

  explicit ParticleCoordinates(Integer n_electron) : coordinates(n_electron) { }
  Complex operator()(Integer i_electron) const { return coordinates(i_electron); }
  Complex& operator()(Integer i_electron) { return coordinates(i_electron); }
};

} // namespace pome