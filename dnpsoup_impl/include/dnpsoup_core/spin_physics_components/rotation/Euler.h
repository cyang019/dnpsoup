#ifndef DNPSOUP_EULER_H
#define DNPSOUP_EULER_H

#include "dnpsoup_core/spin_physics_components/rotation/RotationType.h"
#include "dnpsoup_core/spin_physics_components/rotation/Quaternion.h"
#include <type_traits>
#include <iostream>
#include <iomanip>		// setprecision


namespace dnpsoup {

  // ZYZ rotations
  template<typename T = ActiveRotation>
  class Euler {
  public:
      Euler() : m_alpha(0.0), m_beta(0.0), m_gamma(0.0) {}
      Euler(double a, double b, double g) : m_alpha(a), m_beta(b), m_gamma(g) {}
      Euler(const Euler<T> &);
      Euler(Euler<T> &&) noexcept;
      Euler<T>& operator=(const Euler<T> &);
      Euler<T>& operator=(Euler<T> &&) noexcept;
      ~Euler() {}

      const double& alpha() const { return m_alpha; }
      const double& beta() const { return m_beta; }
      const double& gamma() const { return m_gamma; }

      Euler<T>& alpha(double a) { m_alpha = a; return *this; }
      Euler<T>& beta(double b) { m_beta = b; return *this; }
      Euler<T>& gamma(double g) { m_gamma = g; return *this; }
      bool equals(
          const Euler<T> &rhs,
          double eps=std::numeric_limits<double>::epsilon()) const
      {
        const bool same_a = (m_alpha - rhs.alpha()) > -eps && (m_alpha - rhs.alpha()) < eps;
        const bool same_b = (m_beta - rhs.beta()) > -eps && (m_beta - rhs.beta()) < eps;
        const bool same_g = (m_gamma - rhs.gamma()) > -eps && (m_alpha - rhs.gamma()) < eps;
        return same_a && same_b && same_g;
      }
  private:
      double m_alpha;
      double m_beta;
      double m_gamma;
  };  // class Euler

  // using quaternion to calculate
  template<typename T>
  Euler<T> operator*(const Euler<T> &e1, const Euler<T> &e2);

  template<typename R>
  Euler<R> toEuler(const Quaternion &q);

  template<typename R>
  Quaternion toQuaternion(const Euler<R> &);

  template<typename T>
  std::ostream& operator<<(std::ostream &os, const Euler<T> &euler);
} // namespace dnpsoup


#include "dnpsoup_core/spin_physics_components/rotation/EulerImpl.hpp"

#endif
