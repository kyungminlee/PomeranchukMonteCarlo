#pragma once

#include <random>
#include <limits>
#include "basics.h"
#include "particle.h"
//#include "naive_theta.h"

namespace pome
{
class CflWavefunction;
struct CflCache;

// Thin container? or a proper object?
struct CflCache
{
  kore::array::Array<Complex, 3> theta, theta_old;
  kore::array::Array<Complex, 2> slater;
  Real log_determinant_abs_squared;
  Real log_determinant_abs_squared_old; // with theta_old
  bool dirty; // true if theta different from theta_old or anything like that

  // NOTE: slater gets destroyed everytime determinant is called, so don't even bother.

  explicit CflCache(Integer n_electron)
    // array
    : theta(n_electron, n_electron, n_electron)
    , theta_old(n_electron, n_electron, n_electron)
    , slater(n_electron, n_electron)
    // value cache
    , log_determinant_abs_squared(-std::numeric_limits<Real>::infinity())
    , log_determinant_abs_squared_old(-std::numeric_limits<Real>::infinity())
    // flags
    , dirty(true)
  {
    theta.fill(std::numeric_limits<double>::quiet_NaN());
    theta_old.fill(std::numeric_limits<double>::quiet_NaN());
    slater.fill(std::numeric_limits<double>::quiet_NaN());
  }

  void show(std::ostream& os = std::cout) const
  {
#define SHOWVAR(x) { os << #x << " : " << (x) << std::endl; }
    SHOWVAR(theta);
    SHOWVAR(theta_old);
    SHOWVAR(slater);
    SHOWVAR(dirty);
    SHOWVAR(log_determinant_abs_squared);
    SHOWVAR(log_determinant_abs_squared_old);
#undef SHOWVAR
  }
}; // class CflCache


class CflWavefunction
{
public:

  CflWavefunction(Geometry const & geometry,
                  Integer n_electron,
                  const Complex* displacements)
    : geometry_(geometry), n_electron_(n_electron)

    , displacements_(n_electron)
    , displacement_differences_(n_electron, n_electron)
    , displacement_difference_magnitudes_(n_electron, n_electron)
    , maximum_displacement_difference_(4.0) // TODO

    , displacement_kappas_(n_electron)
    , displacement_kappa_differences_(n_electron, n_electron)
    , displacement_kappa_difference_magnitudes_(n_electron, n_electron)
    , maximum_displacement_kappa_difference_(4.0) // TODO

    , displacement_inverse_kappas_(n_electron)
    , displacement_inverse_kappa_differences_(n_electron, n_electron)
    , displacement_inverse_kappa_difference_magnitudes_(n_electron, n_electron)
    , maximum_displacement_inverse_kappa_difference_(4.0) // TODO

    , momentums_(n_electron)
    , momentum_differences_(n_electron, n_electron)
    , momentum_difference_magnitudes_(n_electron, n_electron)
    , maximum_momentum_difference_(4.0) // TODO

    , distribution_electron_number_(0, n_electron - 1)
    , distribution_01_(0.0, 1.0)
    , distribution_angle_(0.0, 2*M_PI)
  {
    assert(n_electron > 0);
    init_displacements(displacements);
  }


  void init_displacements(Complex const * displacements)
  {
    using namespace kore::array;
    auto ne = n_electron_;
    auto kappa = geometry_.kappa();

    momentum_total_ = 0.0;
    displacement_total_ = 0.0;
    displacement_kappa_total_ = 0.0;
    displacement_inverse_kappa_total_ = 0.0;

    for (Integer i = 0; i < ne; ++i) {
      auto dvec = displacements[i];
      //Complex kvec(dvec.imag() * 2.0, -dvec.real() * 2); //!< d=ik/2,  k=-2id   if ell_B = 1/sqrt(2)
      Complex kvec(dvec.imag(), -dvec.real()); //!< d=ik,  k=-id   (k1, k2) = (d2, -d1)

      auto dvec_kappa = dvec * kappa;
      auto dvec_inverse_kappa = dvec / kappa;

      assert(geometry_.is_valid_momentum(kvec));
      assert(geometry_.is_valid_displacement(dvec));

      momentums_(i) = kvec;
      displacements_(i) = dvec;
      displacement_kappas_(i) = dvec_kappa;
      displacement_inverse_kappas_(i) = dvec_inverse_kappa;

      momentum_total_ += kvec;
      displacement_total_ += dvec;
      displacement_kappa_total_ += dvec_kappa;
      displacement_inverse_kappa_total_ += dvec_inverse_kappa;
    }

    {
      auto inv_ne = 1.0 / static_cast<Real>(n_electron_);
      momentum_mean_ = momentum_total_ * inv_ne;
      displacement_mean_ = displacement_total_ * inv_ne;
      displacement_kappa_mean_ = displacement_kappa_total_ * inv_ne;
      displacement_inverse_kappa_mean_ = displacement_inverse_kappa_total_ * inv_ne;
    }

    auto map_difference = [ne](Array<Complex, 1> const & src, 
                               Array<Complex, 2>       & tgt,
                               Array<Real ,2>          & tgt_mag) {
      for (Integer i = 0; i < ne; ++i) {
        for (Integer j = 0; j < ne; ++j) {
          auto q = src(i) - src(j);
          tgt(i, j) = q;
          tgt_mag(i, j) = std::abs(q);
        }
      }
    };

    map_difference(momentums_,
                   momentum_differences_,
                   momentum_difference_magnitudes_);
    map_difference(displacements_,
                   displacement_differences_,
                   displacement_difference_magnitudes_);
    map_difference(displacement_kappas_,
                   displacement_kappa_differences_,
                   displacement_kappa_difference_magnitudes_);
    map_difference(displacement_inverse_kappas_,
                   displacement_inverse_kappa_differences_,
                   displacement_inverse_kappa_difference_magnitudes_);
  }


  //! \name MonteCarloEngine conformant functions
  //@{
  template <typename RandGen>
  ParticleMove candidate(ParticleCoordinates const & coordinates,
                         CflCache const & cache,
                         ParticleMoveOpt const & opt,
                         RandGen& gen) const
  {
    auto ie = distribution_electron_number_(gen);
    auto z_old = coordinates(ie);
#if 1
    auto radius = distribution_01_(gen) * opt.radius;
    auto angle = distribution_angle_(gen);
    auto dz = Complex(::cos(angle), ::sin(angle)) * radius;
#endif
#if 0
    auto dx = (2.0 * distribution_01_(gen) - 1.0) * opt.radius;
    auto dy = (2.0 * distribution_01_(gen) - 1.0) * opt.radius;
    auto dz = Complex(dx, dy);
#endif
    auto z_new = geometry_.mod_kappa_z(z_old + dz);
    return{ ie, z_new, z_old };
  }

  void update(ParticleMove const & move,
              ParticleCoordinates & coordinates,
              CflCache & cache) const
  {
    update_theta(move.i_electron, move.z_new, coordinates, cache);
    coordinates.coordinates(move.i_electron) = move.z_new;
    generate_slater(coordinates, cache);
  }

  void commit(ParticleMove const & move,
              ParticleCoordinates & coordinates,
              CflCache & cache) const
  {
    coordinates(move.i_electron) = move.z_new;
    commit_theta(move.i_electron, cache);
  }

  void revert(ParticleMove const & move,
              ParticleCoordinates & coordinates,
              CflCache & cache) const
  {
    coordinates(move.i_electron) = move.z_old;
    revert_theta(move.i_electron, cache);
  }


  //! Probability 
  //!
  //! \f[ P(\{z_i\}) \propto \left\vert \det_{i,j} \psi_{i} (\{z_{j}\}) \right\vert^2 \f]
  Real acceptance(ParticleMove const & move,
                  ParticleCoordinates const & coords,
                  CflWavefunction const & wavefunction,
                  CflCache const & cache) const
  {
    Real logprob_old = cache.log_determinant_abs_squared_old;
    Real logprob_new = cache.log_determinant_abs_squared;
    if (logprob_new >= logprob_old) { return 1.0; }
    return std::exp(logprob_new - logprob_old);
  }
  //@}

  void generate_theta(ParticleCoordinates const & coordinates,
                      CflCache & cache) const
  {
    //static const Complex I(0.0, 1.0);
    //auto kappa = geometry_.kappa();
    auto nel = n_electron_;
    auto const & zs = coordinates;
    auto const & kd = displacement_kappas_;
    auto kd0 = displacement_kappa_mean_;
    auto inv_tp1 = geometry_.inverse_theta1_prime();
    for (Integer i = 0; i < nel; ++i) {
      for (Integer j = 0; j < nel; ++j) {
        for (Integer k = 0; k < nel; ++k) {
          if (k != i) {
            auto w = zs(i) - zs(k) + 2.0 * (kd(j) - kd0);
            //w = geometry_.mod_kappa_z(w);
            auto phi = geometry_.theta1(w) * inv_tp1;
            cache.theta(i, k, j) = phi;
            cache.theta_old(i, k, j) = phi;
          } // if k != i
        } // for k
      } // for j
    } // for i
    cache.dirty = false;
  }

  void update_theta(Integer i_electron, Complex z_new,
                    ParticleCoordinates const & coordinates,
                    CflCache & cache) const
  {
    auto nel = n_electron_;
    auto const & zs = coordinates;
    auto const & kd = displacement_kappas_;
    auto kd0 = displacement_kappa_mean_;
    auto inv_tp1 = geometry_.inverse_theta1_prime();
#if 0 
    {
      Integer k = i_electron;
      for (Integer i = 0; i < nel; ++i) {
        if (k != i) {
          for (Integer j = 0; j < nel; ++j) {
            auto w = zs(i) - z_new + 2.0 * (kd(j) - kd0);
            //w = geometry_.mod_kappa_z(w);
            auto phi = geometry_.theta1(w) * inv_tp1;
            cache.theta(i, k, j) = phi;
          } // for j
        }  // if k != i
      } // for i
    }

    {
      Integer i = i_electron;
      for (Integer k = 0; k < nel; ++k) {
        if (k != i) {
          for (Integer j = 0; j < nel; ++j) {
            auto w = z_new - zs(k) + 2.0 * (kd(j) - kd0);
            //w = geometry_.mod_kappa_z(w);
            auto phi = geometry_.theta1(w) * inv_tp1;
            cache.theta(i, k, j) = phi;
          }  // for j
        } // if k != i
      } // for k
    } // for i
#endif

    {
      for (Integer i = 0; i < nel; ++i) {
        if (i != i_electron) {
          for (Integer j = 0; j < nel; ++j) {
            {
              auto w = zs(i) - z_new + 2.0 * (kd(j) - kd0);
              //w = geometry_.mod_kappa_z(w);
              auto phi = geometry_.theta1(w) * inv_tp1;
              cache.theta(i, i_electron, j) = phi;
            }
            {
              auto w = z_new - zs(i) + 2.0 * (kd(j) - kd0);
              //w = geometry_.mod_kappa_z(w);
              auto phi = geometry_.theta1(w) * inv_tp1;
              cache.theta(i_electron, i, j) = phi;
            }
          } // for j
        }  // if i != i_electron
      } // for i
    }
    cache.dirty = true;
  }


  void commit_theta(Integer i_electron,
                    CflCache& cache) const
  {
    //auto const& zs = coordinates;
    //zs(i_electron) = z_new;
    auto const nel = n_electron_;
    for (Integer i = 0; i < nel; ++i) {
      if (i != i_electron) {
        for (Integer j = 0; j < nel; ++j) {
          cache.theta_old(i, i_electron, j) = cache.theta(i, i_electron, j);
          cache.theta_old(i_electron, i, j) = cache.theta(i_electron, i, j);
        }
      }
    }
    cache.log_determinant_abs_squared_old = cache.log_determinant_abs_squared;
    cache.dirty = false;
  }


  void revert_theta(Integer i_electron,
                    CflCache& cache) const
  {
    auto const nel = n_electron_;
    for (Integer i = 0; i < nel; ++i) {
      if (i != i_electron) {
        for (Integer j = 0; j < nel; ++j) {
          cache.theta(i, i_electron, j) = cache.theta_old(i, i_electron, j);
          cache.theta(i_electron, i, j) = cache.theta_old(i_electron, i, j);
        }
      }
    }
    //cache.log_determinant_abs_squared = cache.log_determinant_abs_squared_old;
    //<- unnecessary, but let's keep it here.
    cache.dirty = false;
  }


  void generate_slater(ParticleCoordinates const & coordinates,
                       CflCache& cache) const
  {
    static const Complex I(0.0, 1.0);
    const Complex inv_kappa = geometry_.inverse_kappa();
    auto const nel = n_electron_;
    auto const & zs = coordinates; // kappa z
    auto const inv_nel = 1.0 / static_cast<Real>(nel);

    auto y_total = 0.0;
    auto y_squared_total = 0.0;
    for (Integer i = 0; i < nel; ++i) {
      auto y = zs(i).imag();
      y_total += y;
      y_squared_total += y * y;
    }
    y_total *= inv_kappa.real();
    y_squared_total *= inv_kappa.real() * inv_kappa.real();

    auto xp = std::exp((y_total * y_total * inv_nel
                        - y_squared_total
                        + 2.0 * displacement_total_.imag() * y_total * inv_nel)
                       * 0.5 * inv_nel / (nel - 1));
    // k = -id   (k1, k2) = (d2, -d1)

    for (Integer i = 0; i < nel; ++i) {
      for (Integer j = 0; j < nel; ++j) {
        Complex factor(1.0, 0.0);
        for (Integer k = 0; k < nel; ++k) {
          if (k != i) {
            factor *= cache.theta(i, k, j) * xp;
          } // if (k != i)
        } // for k
        cache.slater(i, j)
          = factor * std::exp( I * zs(i) * displacement_inverse_kappas_(j).imag() );
      } // for j
    } // for i
    auto ldsq = kore::array::linalg::log_determinant_abs_squared(cache.slater);
    cache.log_determinant_abs_squared = ldsq;
  }

  //! \name Getter
  //@{
  Complex momentum_difference(Integer i, Integer j) const
  {
    assert(0 <= i && i < n_electron_);
    assert(0 <= j && j < n_electron_);
    return momentum_differences_(i, j);
  }

  Real momentum_difference_magnitude(Integer i, Integer j) const
  {
    assert(0 <= i && i < n_electron_);
    assert(0 <= j && j < n_electron_);
    return momentum_difference_magnitudes_(i, j);
  }

  Integer n_electron() const { return n_electron_; }
  //@}

  kore::array::Array<Complex const, 1> const & momentums() const { return momentums_; }
  kore::array::Array<Complex const, 2> const & momentum_differences() const { return momentum_differences_; }
  kore::array::Array<Real const, 2> const & momentum_difference_magnitudes() const { return momentum_difference_magnitudes_; }
  Complex momentum_total() const { return momentum_total_; }
  Complex momentum_mean() const { return momentum_mean_; }

  kore::array::Array<Complex const, 1> const & displacements() const { return displacements_; }
  kore::array::Array<Complex const, 2> const & displacement_differences() const { return displacement_differences_; }
  kore::array::Array<Real const, 2> const & displacement_difference_magnitudes() const { return displacement_difference_magnitudes_; }
  Complex displacement_total() const { return displacement_total_; }
  Complex displacement_mean() const { return displacement_mean_; }

  kore::array::Array<Complex const, 1> const & displacement_kappas() const { return displacement_kappas_; }
  kore::array::Array<Complex const, 2> const & displacement_kappa_differences() const { return displacement_kappa_differences_; }
  kore::array::Array<Real const, 2> const & displacement_kappa_difference_magnitudes() const { return displacement_kappa_difference_magnitudes_; }
  Complex displacement_kappa_total() const { return displacement_kappa_total_; }
  Complex displacement_kappa_mean() const { return displacement_kappa_mean_; }

  kore::array::Array<Complex const, 1> const & displacement_inverse_kappas() const { return displacement_inverse_kappas_; }
  kore::array::Array<Complex const , 2> const & displacement_inverse_kappa_differences() const { return displacement_inverse_kappa_differences_; }
  kore::array::Array<Real const, 2> const & displacement_inverse_kappa_difference_magnitudes() const { return displacement_inverse_kappa_difference_magnitudes_; }
  Complex displacement_inverse_kappa_total() const { return displacement_inverse_kappa_total_; }
  Complex displacement_inverse_kappa_mean() const { return displacement_inverse_kappa_mean_; }


  void show(std::ostream& os = std::cout) const
  {
#define SHOWVAR(x) { os << #x << " : " << (x) << std::endl; }
    SHOWVAR(n_electron_);
    SHOWVAR(momentums_);
    SHOWVAR(momentum_differences_);
    SHOWVAR(momentum_difference_magnitudes_);
    SHOWVAR(momentum_total_);
    SHOWVAR(momentum_mean_);
    SHOWVAR(maximum_momentum_difference_);
#undef SHOWVAR
  }

private:
  Geometry const & geometry_;
  Integer n_electron_;


  kore::array::Array<Complex, 1> displacements_;
  kore::array::Array<Complex, 2> displacement_differences_;
  kore::array::Array<Real, 2> displacement_difference_magnitudes_;
  Complex displacement_total_;
  Complex displacement_mean_;
  Real maximum_displacement_difference_;

  kore::array::Array<Complex, 1> displacement_kappas_;
  kore::array::Array<Complex, 2> displacement_kappa_differences_;
  kore::array::Array<Real, 2> displacement_kappa_difference_magnitudes_;
  Complex displacement_kappa_total_;
  Complex displacement_kappa_mean_;
  Real maximum_displacement_kappa_difference_;

  kore::array::Array<Complex, 1> displacement_inverse_kappas_;
  kore::array::Array<Complex, 2> displacement_inverse_kappa_differences_;
  kore::array::Array<Real, 2> displacement_inverse_kappa_difference_magnitudes_;
  Complex displacement_inverse_kappa_total_;
  Complex displacement_inverse_kappa_mean_;
  Real maximum_displacement_inverse_kappa_difference_;

  kore::array::Array<Complex, 1> momentums_;
  kore::array::Array<Complex, 2> momentum_differences_;
  kore::array::Array<Real, 2> momentum_difference_magnitudes_;
  Complex momentum_total_;
  Complex momentum_mean_;
  Real maximum_momentum_difference_;


  // TODO: WHY not const?
  mutable std::uniform_int_distribution<Integer> distribution_electron_number_;
  mutable std::uniform_real_distribution<Real> distribution_01_;
  mutable std::uniform_real_distribution<Real> distribution_angle_;
};


} // namespace pome

