#pragma once

#include <cmath>
#include "geometry.h"

namespace pome {

inline
Geometry::Geometry(Real ell1, Real ell2, Real angle)
  : ell1_(ell1), ell2_(ell2), angle_(angle)
  , aspect_ratio_(ell1 / ell2)
  , tau_(aspect_ratio_ * ::cos(angle), aspect_ratio_ * ::sin(angle))
  , i_pi_tau_(-tau_.imag() * M_PI, tau_.real() * M_PI)
  , kappa_(M_PI / ell1)
  , inverse_kappa_(1.0 / kappa_)
  , inverse_lattice_basis_(2)
  , theta1_function_(tau_)
  , theta2_function_(tau_)
  , theta3_function_(tau_)
  , theta4_function_(tau_)
  , theta1_prime_(0.0, 0.0)
  , displacements_(0)
{
  assert(ell1 > 0);
  assert(ell2 > 0);
  assert(angle > 0);
  assert(angle <= M_PI / 2);
  assert(std::imag(tau_) > 0);

  Real area = ell1 * ell2 * ::sin(angle);
  Real n_flux_real = area / (2 * M_PI);
  Integer n_flux_integer = static_cast<Integer>(std::round(n_flux_real));

  assert(::fabs(static_cast<Real>(n_flux_integer) - n_flux_real) < 1E-13);
  n_flux_ = n_flux_integer;
  assert(n_flux_ % 2 == 0);

  Complex zero(0.0, 0.0);
  theta1_prime_ = theta2(zero) * theta3(zero) * theta4(zero);

  // L = | L1x  L2x |
  //     | L1y  L2y |
  kore::array::Array<Real, 2> lattice(2, 2);
  lattice(0, 0) = ell1_;
  lattice(1, 0) = 0.0;
  lattice(0, 1) = ell2_ * ::cos(angle);
  lattice(1, 1) = ell2_ * ::sin(angle);

  // G = | g1x  g1y |
  //     | g2x  g2y |
  // such that G . L = L . G = identity(2)
  auto inverse_lattice = kore::array::linalg::inverse(lattice);

  inverse_lattice_basis_(0) = Complex(inverse_lattice(0, 0), inverse_lattice(0, 1));
  inverse_lattice_basis_(1) = Complex(inverse_lattice(1, 0), inverse_lattice(1, 1));

  // TODO: fix this
  reference_k_ = 0.0;
  reference_z0_ = 0.0;
  init_displacements();
}


inline bool
Geometry::is_valid_displacement(Complex dvec) const
{
  // Given a system of size L1 x L2 (at angle theta),  valid momenta are
  // 
  //
  //
  // $N_{\Phi} \bfd_{i} = m \bfL_1 + n \bfL_2
  dvec *= static_cast<Real>(n_flux_);
  Real m1_real = inverse_lattice_basis_(0).real() * dvec.real()
    + inverse_lattice_basis_(0).imag() * dvec.imag();
  Real m2_real = inverse_lattice_basis_(1).real() * dvec.real()
    + inverse_lattice_basis_(1).imag() * dvec.imag();
  Integer m1_int = static_cast<Integer>(std::round(m1_real));
  Integer m2_int = static_cast<Integer>(std::round(m2_real));
  if (std::fabs(m1_real - m1_int) < 1E-13
    && std::fabs(m2_real - m2_int) < 1E-13) {
    return true;
  }
  else {
    return false;
  }
}


// TODO: k vector generation and validation
inline bool
Geometry::is_valid_momentum(Complex k) const
{
  Complex dvec(k.imag(), -k.real());
  return is_valid_displacement(dvec);
}


//!       +-----------+
//!      /           /
//!     /     z     / L1 tau
//!    /           /
//!   +-----------+
//!        L1
inline Complex
Geometry::mod_z(Complex z) const {
  z /= ell1_;
  Real y_cell = ::floor(z.imag() / tau_.imag() + 1E-15);
  z -= y_cell * tau_;
  Real offset = z.imag() * tau_.real() / tau_.imag();
  Real x_cell = ::floor(z.real() - offset + 1E-15);
  return Complex(ell1_ * (z.real() - x_cell), ell1_ * z.imag());
}


//!       +-----------+
//!      /           /
//!     /  kappa z  / pi tau
//!    /           /
//!   +-----------+
//!        pi
inline Complex
Geometry::mod_kappa_z(Complex z) const {
  z /= M_PI;
  Real y_cell = ::floor(z.imag() / tau_.imag() + 1E-15);
  z -= y_cell * tau_;
  Real offset = z.imag() * tau_.real() / tau_.imag();
  Real x_cell = ::floor(z.real() - offset + 1E-15);
  return Complex(M_PI * (z.real() - x_cell), M_PI * z.imag());
}


inline void
Geometry::init_displacements()
{
  Complex a1(ell1_ / n_flux_, 0.0);
  Complex a2(ell2_ * ::cos(angle_) / n_flux_, ell2_ * ::sin(angle_) / n_flux_);

  kore::array::Array<Complex, 2u> dvectors(n_flux_, n_flux_);
  for (Integer i1 = 0; i1 < n_flux_; ++i1) {
    for (Integer i2 = 0; i2 < n_flux_; ++i2) {
      dvectors(i1, i2) = a1 * static_cast<Real>(i1) + a2 * static_cast<Real>(i2);
      assert(is_valid_momentum(dvectors(i1, i2)));
    }
  }
  displacements_ = std::move(dvectors);
}



} // namespace pome
