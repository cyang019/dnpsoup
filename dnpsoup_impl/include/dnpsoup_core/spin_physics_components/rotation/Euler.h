#ifndef DNPSOUP_EULER_H
#define DNPSOUP_EULER_H

#include <type_traits>


namespace dnpsoup {
  class Quaternion;

  class ActiveRotation {};

  class PassiveRotation {};

  template<typename T> struct is_rotation_type : std::false_type {};
  template<> struct is_rotation_type<ActiveRotation> : std::true_type {};
  template<> struct is_rotation_type<PassiveRotation> : std::true_type {};

  // ZYZ rotations
  template<typename T = ActiveRotation>
  class Euler {
  public:
      Euler() : m_alpha(0.0), m_beta(0.0), m_gamma(0.0) {}
      Euler(double a, double b, double g) : m_alpha(a), m_beta(b), m_gamma(g) {}
      Euler(const Quaternion &);
      Euler(const Euler<T> &) = default;
      Euler(Euler<T> &&) noexcept = default;
      Euler<T>& operator=(const Euler<T> &) = default;
      Euler<T>& operator=(Euler<T> &&) noexcept = default;
      ~Euler() {}

      double alpha() const { return m_alpha; }
      double beta() const { return m_beta; }
      double gamma() const { return m_gamma; }

      Euler<T>& alpha(double a) { m_alpha = a; return *this; }
      Euler<T>& beta(double b) { m_beta = b; return *this; }
      Euler<T>& gamma(double g) { m_gamma = g; return *this; }
  private:
      double m_alpha;
      double m_beta;
      double m_gamma;
  };  // class Euler

  // using quaternion to calculate
  template<typename T>
  Euler<T> operator*(const Euler<T> &e1, const Euler<T> &e2);
} // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/rotation/EulerImpl.hpp"

#endif
