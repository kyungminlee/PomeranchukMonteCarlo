#pragma once

#include "basics.h"
#include "theta.h"
#include <kore/array/array.h>
#include <kore/array/array_linalg.h>

namespace pome {

/*! \verbatim
        +------------+
       /            /
      /            / L2
     / angle      /
    +------------+
         L1
\endverbatim

Use \f$ \ell_{B} == 1/\sqrt{2} \f$
*/
class Geometry {
private:
  Real ell1_, ell2_, angle_;

  // Useful variables
  Real aspect_ratio_;
  Complex tau_, i_pi_tau_;
  Real kappa_;
  Real inverse_kappa_;
  Integer n_flux_;
  Complex theta1_prime_; // theta1'(0) = theta2(0) * theta3(0) * theta4(0)
  Complex inverse_theta1_prime_;

  kore::array::Array<Complex, 1> inverse_lattice_basis_;
  kore::array::Array<Complex, 2> displacements_;

  //Real reference_k_;  //!< ala Haldane and Rezayi, PRB 31, 2529 (1985)
  //Real reference_z0_; //!< ala Haldane and Rezayi, PRB 31, 2529 (1985)
                      // Cache Objects
  EllipticTheta1<Real> theta1_function_;
  EllipticTheta2<Real> theta2_function_;
  EllipticTheta3<Real> theta3_function_;
  EllipticTheta4<Real> theta4_function_;

public:
  //! ell1, ell2 > 0 and 0 < angle <= PI/2
  Geometry(Real ell1, Real ell2, Real angle);

  bool is_valid_displacement(Complex dvec) const;

  // TODO: k vector generation and validation
  bool is_valid_momentum(Complex k) const;

  //! \name modular
  //@{
  Complex mod_z(Complex z) const;
  Complex mod_kappa_z(Complex z) const;
  //@}

  //! \name theta functions
  //@{
  Complex theta1(Complex z) const { return theta1_function_(z); } //!< \f$\vartheta_{1}(z|\tau)\f$
  Complex theta2(Complex z) const { return theta2_function_(z); } //!< \f$\vartheta_{2}(z|\tau)\f$
  Complex theta3(Complex z) const { return theta3_function_(z); } //!< \f$\vartheta_{3}(z|\tau)\f$
  Complex theta4(Complex z) const { return theta4_function_(z); } //!< \f$\vartheta_{4}(z|\tau)\f$
  Complex theta1_prime() const { return theta1_prime_; }          //!< \f$\vartheta_{1}'(0|\tau)\f$
  Complex inverse_theta1_prime() const { return inverse_theta1_prime_; }          //!< \f$1/\vartheta_{1}'(0|\tau)\f$
  //@}

public:
  //! \name Getter.
  //@{
  auto displacements() const -> const decltype(displacements_) & {
    return displacements_;
  }

  Real ell1() const { return ell1_; }                   //!< \f$L_{1}\f$
  Real ell2() const { return ell2_; }                   //!< \f$L_{2}\f$
  Real L1() const { return ell1_; }                     //!< \f$L_{1}\f$
  Real L2() const { return ell2_; }                     //!< \f$L_{2}\f$
  Integer n_flux() const { return n_flux_; }            //!< \f$N_{\Phi} \equiv L_1 L_2 \sin \theta /(2\pi)\f$
  Real aspect_ratio() const { return ell2_ / ell1_; }   //!< \f$L_{2}/L_{1}\f$
  Complex tau() const { return tau_; }                  //!< \f$\tau \equiv L_2 e^{i\theta}/L_1 \f$
  Complex i_pi_tau() const { return i_pi_tau_; }        //!< \f$i \pi \tau\f$
  Real kappa() const { return kappa_; }                 //!< \f$\kappa\f$
  Real inverse_kappa() const { return inverse_kappa_; } //!< \f$1/\kappa\f$
  //@}

public:
  void debug_show() const;

private:

  //! initialize displacements 
  //!
  //! \f$ \mathbf{d}_i \f$
  void init_displacements();

};

} // namespace pome

#include "geometry_impl.h"