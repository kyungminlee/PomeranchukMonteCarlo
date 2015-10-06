#pragma once

#define _USE_MATH_DEFINES
#include <cassert>
#include <cmath>
#include <complex>
#include <array>

inline
std::complex<double> naive_theta1(std::complex<double> tau, std::complex<double> z)
{
  std::complex<double> I(0, 1.0), theta(0.0, 0.0);
  for (int k = 0; k < 64; ++k) {
    int isign = (k % 2) ? -1 : 1;
    std::complex<double> zpr = static_cast<double>(k*(k + 1))*M_PI*(I*tau)
      + std::log(std::sin((static_cast<double>(2 * k + 1)*z)));
    std::complex<double> next_term = static_cast<double>(isign)*std::exp(zpr);
    theta += next_term;
    if (std::abs(next_term) < 1E-15) break;
  }
  theta *= 2.0*(std::exp((M_PI*I*tau) / 4.0));
  return theta;
}

inline
std::complex<double> naive_theta2(std::complex<double> tau, std::complex<double> z)
{
  std::complex<double> I(0, 1.0), theta(0.0, 0.0);
  for (int k = 0; k < 64; ++k) {
    std::complex<double> zpr = static_cast<double>(k*(k + 1))*M_PI*(I*tau)
      + std::log(std::cos((static_cast<double>(2 * k + 1)*z)));
    std::complex<double> next_term = std::exp(zpr);
    theta += next_term;
    if (std::abs(next_term) < 1E-15) break;
  }
  theta *= 2.0*(std::exp((M_PI*I*tau) / 4.0));
  return theta;
}

inline
std::complex<double> naive_theta3(std::complex<double> tau, std::complex<double> z)
{
  std::complex<double> I(0, 1.0), theta(0.0, 0.0);
  for (int k = 1; k < 64; ++k) {
    std::complex<double> zpr = static_cast<double>(k*k)*M_PI*(I*tau)
      + std::log(std::cos((static_cast<double>(2 * k)*z)));
    std::complex<double> next_term = std::exp(zpr);
    theta += next_term;
    if (std::abs(next_term) < 1E-15) break;
  }
  theta = 1.0 + 2.0*theta;
  return theta;
}

inline
std::complex<double> naive_theta4(std::complex<double> tau, std::complex<double> z)
{
  std::complex<double> I(0, 1.0), theta(0.0, 0.0);
  for (int k = 1; k < 64; ++k) {
    int isign = (k % 2) ? -1 : 1;
    std::complex<double> zpr = static_cast<double>(k*k)*M_PI*(I*tau)
      + std::log(std::cos((static_cast<double>(2 * k)*z)));
    std::complex<double> next_term = static_cast<double>(isign)*std::exp(zpr);
    theta += next_term;
    if (std::abs(next_term) < 1E-15) break;
  }
  theta = 1.0 + 2.0*theta;
  return theta;
}

