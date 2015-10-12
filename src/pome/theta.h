#pragma once
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES 1
#endif
#include <cassert>
#include <cmath>
#include <complex>
#include <array>

//! \f$ \vartheta_1(z | \tau) \f$
template <typename _RealType = double, size_t CacheSize = 128>
class EllipticTheta1
{
public:
  typedef _RealType RealType;
  typedef std::complex<RealType> ComplexType;
  static_assert(std::is_same<RealType, double>::value, "RealType should be double");

  explicit EllipticTheta1(ComplexType tau)
    : tau_(tau)
  {
    static const ComplexType complex_eye(0.0, 1.0);
    RealType rtau = std::real(tau);
    //rtau = std::fmod(rtau + ::ceil(fabs(rtau))*2.0 + 1.0, 2.0) - 1.0;
    //tau_mod2_ = ComplexType(rtau, std::imag(tau_));

    ComplexType foo = 0;
    for (size_t i = 0; i < CacheSize; ++i) {
      cache_two_pi_i_tau_[i] = foo;
      foo += 2.0 * M_PI * complex_eye * tau_;
    }
    overall_factor_ = complex_eye * std::exp(0.25 * M_PI * complex_eye * tau_);
  }

  ComplexType operator()(ComplexType z) const
  {
    //static const ComplexType I(0.0, 1.0);
    ComplexType theta = 0.0;
    ComplexType two_iz(-2.0 * z.imag(), 2.0 * z.real()); // = 2 I z

    ComplexType z1(std::imag(z), -std::real(z)); // = -I * z
    ComplexType z2 = -z1;

    for (size_t i = 0; i < CacheSize; ++i) {
      ComplexType foo = cache_two_pi_i_tau_[i];
      z1 += foo + two_iz;
      z2 += foo - two_iz;
      //ComplexType next_term = (std::exp(z1) - std::exp(z2));
      ComplexType next_term =
        ComplexType(::cos(z1.imag()), ::sin(z1.imag())) * ::exp(z1.real()) -
        ComplexType(::cos(z2.imag()), ::sin(z2.imag())) * ::exp(z2.real());
      //theta += (i % 2 == 0) ? -next_term : next_term;
      theta += (i % 2) ? next_term : -next_term;
      if (::fabs(next_term.real()) + ::fabs(next_term.imag()) < 1E-15) { break; } 
      // faster than std::abs                                                                           
      //if (std::abs(next_term) < 1E-15) { break; }
    }
    theta *= overall_factor_;
    return theta;
  }

private:
  ComplexType tau_, overall_factor_;
  std::array<ComplexType, CacheSize> cache_two_pi_i_tau_;
  std::array<ComplexType, CacheSize> cache_exponential_tau_;
  // <- (-1)^n exp(n (n+1) i pi tau)
};

//! \f$ \vartheta_{1}'(z | \tau) \f$
template <typename _RealType = double, size_t CacheSize = 64>
class EllipticThetaPrime1
{
public:
  typedef _RealType RealType;
  typedef std::complex<RealType> ComplexType;

  explicit EllipticThetaPrime1(ComplexType tau)
    : tau_(tau)
  {
    static const ComplexType I(0, 1);

    for (size_t i = 0; i < CacheSize; ++i) {
      ComplexType val = (i % 2 == 0) ? 1 : -1;
      val *= std::exp(I * M_PI * tau * RealType(i*(i + 1)));
      cache_exponential_tau_[i] = val;
    }
  }

  ComplexType operator()(int i_pow, ComplexType z) const {
    assert(i_pow >= 0);
    static const ComplexType I(0, 1);
    ComplexType theta = 0.0;
    for (size_t i = 0; i < CacheSize; ++i) {
      RealType alpha = RealType(2 * i + 1);
      RealType alpha_pow = std::pow(alpha, i_pow);
      ComplexType next_term;
      switch (i_pow % 4) {
      case 0:
        next_term = cache_exponential_tau_[i] * std::sin(z * alpha) * alpha_pow;
        break;
      case 1:
        next_term = cache_exponential_tau_[i] * std::cos(z * alpha) * alpha_pow;
        break;
      case 2:
        next_term = -cache_exponential_tau_[i] * std::sin(z * alpha) * alpha_pow;
        break;
      case 3:
        next_term = -cache_exponential_tau_[i] * std::cos(z * alpha) * alpha_pow;
        break;
      }
      theta += next_term;
      if (std::abs(next_term) < 1E-15) { break; }
    }
    theta *= 2.0 * std::exp(0.25 * M_PI * I);
    return theta;
  }

private:
  ComplexType tau_;
  std::array<ComplexType, CacheSize> cache_exponential_tau_;
  // <- (-1)^n exp(n (n+1) i pi tau)
};

//! \f$ \vartheta_2(z | \tau) \f$
template <typename _RealType = double, size_t CacheSize = 64>
class EllipticTheta2
{
public:
  typedef _RealType RealType;
  typedef std::complex<RealType> ComplexType;
  static_assert(std::is_same<RealType, double>::value, "RealType should be double");

  explicit EllipticTheta2(ComplexType tau)
    : tau_(tau)
  {
    static const ComplexType complex_eye(0.0, 1.0);
    ComplexType foo = 0;
    for (size_t i = 0; i < CacheSize; ++i) {
      cache_two_pi_i_tau_[i] = foo;
      foo += 2.0 * M_PI * complex_eye * tau_;
    }
    overall_factor_ = std::exp(0.25 * M_PI * complex_eye * tau_);
  }

  ComplexType operator()(ComplexType z) const
  {
    ComplexType theta = 0.0;
    ComplexType two_iz(-2.0 * z.imag(), 2.0 * z.real()); // = 2 I z

    ComplexType z1(std::imag(z), -std::real(z)); // = -I * z
    ComplexType z2 = -z1;

    for (size_t i = 0; i < CacheSize; ++i) {
      ComplexType foo = cache_two_pi_i_tau_[i];
      z1 += foo + two_iz;
      z2 += foo - two_iz;
      ComplexType next_term =
        ComplexType(::cos(z1.imag()), ::sin(z1.imag())) * ::exp(z1.real()) +
        ComplexType(::cos(z2.imag()), ::sin(z2.imag())) * ::exp(z2.real());
      theta += next_term;
      if (::fabs(next_term.real()) + ::fabs(next_term.imag()) < 1E-15) { break; } // faster than std::abs
                                                                                  //if (std::abs(next_term) < 1E-15) { break; }
    }
    theta *= overall_factor_;
    return theta;
  }

private:
  ComplexType tau_, overall_factor_;
  std::array<ComplexType, CacheSize> cache_two_pi_i_tau_;
  std::array<ComplexType, CacheSize> cache_exponential_tau_;
  // <- (-1)^n exp(n (n+1) i pi tau)
};

//! \f$ \vartheta_3(z | \tau) \f$
template <typename _RealType = double, size_t CacheSize = 64>
class EllipticTheta3
{
public:
  typedef _RealType RealType;
  typedef std::complex<RealType> ComplexType;
  static_assert(std::is_same<RealType, double>::value, "RealType should be double");

  explicit EllipticTheta3(ComplexType tau)
    : tau_(tau)
  {
    static const ComplexType complex_eye(0.0, 1.0);
    ComplexType foo = -M_PI * complex_eye * tau_;
    for (size_t i = 0; i < CacheSize; ++i) {
      cache_two_pi_i_tau_[i] = foo;
      foo += 2.0 * M_PI * complex_eye * tau_;
    }
  }

  ComplexType operator()(ComplexType z) const
  {
    ComplexType theta = 1.0;
    ComplexType two_iz(-2.0 * z.imag(), 2.0 * z.real()); // = 2 I z
    ComplexType z1 = 0, z2 = 0;

    for (int i = 1; i < CacheSize; ++i) {
      ComplexType foo = cache_two_pi_i_tau_[i];
      z1 += foo + two_iz;
      z2 += foo - two_iz;
      ComplexType next_term =
        ComplexType(::cos(z1.imag()), ::sin(z1.imag())) * ::exp(z1.real()) +
        ComplexType(::cos(z2.imag()), ::sin(z2.imag())) * ::exp(z2.real());
      theta += next_term;
      if (::fabs(next_term.real()) + ::fabs(next_term.imag()) < 1E-15) { break; } // faster than std::abs
                                                                                  //if (std::abs(next_term) < 1E-15) { break; }
    }
    return theta;
  }

private:
  ComplexType tau_;
  std::array<ComplexType, CacheSize> cache_two_pi_i_tau_;
  std::array<ComplexType, CacheSize> cache_exponential_tau_;
  // <- (-1)^n exp(n (n+1) i pi tau)
};

//! \f$ \vartheta_4(z | \tau) \f$
template <typename _RealType = double, size_t CacheSize = 64>
class EllipticTheta4
{
public:
  typedef _RealType RealType;
  typedef std::complex<RealType> ComplexType;
  static_assert(std::is_same<RealType, double>::value, "RealType should be double");

  explicit EllipticTheta4(ComplexType tau)
    : tau_(tau)
  {
    static const ComplexType complex_eye(0.0, 1.0);
    ComplexType foo = -M_PI * complex_eye * tau_;
    for (size_t i = 0; i < CacheSize; ++i) {
      cache_two_pi_i_tau_[i] = foo;
      foo += 2.0 * M_PI * complex_eye * tau_;
    }
  }

  ComplexType operator()(ComplexType z) const
  {
    ComplexType theta = 1.0;
    ComplexType two_iz(-2.0 * z.imag(), 2.0 * z.real()); // = 2 I z
    ComplexType z1 = 0, z2 = 0;
    for (int i = 1; i < CacheSize; ++i) {
      ComplexType foo = cache_two_pi_i_tau_[i];
      z1 += foo + two_iz;
      z2 += foo - two_iz;
      ComplexType next_term =
        ComplexType(::cos(z1.imag()), ::sin(z1.imag())) * ::exp(z1.real()) +
        ComplexType(::cos(z2.imag()), ::sin(z2.imag())) * ::exp(z2.real());
      theta += (i % 2) ? -next_term : next_term;
      if (::fabs(next_term.real()) + ::fabs(next_term.imag()) < 1E-15) { break; } // faster than std::abs
                                                                                  //if (std::abs(next_term) < 1E-15) { break; }
    }
    return theta;
  }

private:
  ComplexType tau_;
  std::array<ComplexType, CacheSize> cache_two_pi_i_tau_;
  std::array<ComplexType, CacheSize> cache_exponential_tau_;
  // <- (-1)^n exp(n (n+1) i pi tau)
};
