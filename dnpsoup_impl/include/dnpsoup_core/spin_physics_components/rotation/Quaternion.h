#ifndef DNPSOUP_QUATERNION_H
#define DNPSOUP_QUATERNION_H

#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/common.h"

namespace dnpsoup {
  Quaternion operator+(const Quaternion &, const Quaternion &);
  Quaternion operator-(const Quaternion &, const Quaternion &);
  Quaternion operator*(const Quaternion &, const Quaternion &);

  class Quaternion {
  public:
    friend Quaternion operator+(const Quaternion &, const Quaternion &);
    friend Quaternion operator-(const Quaternion &, const Quaternion &);
    friend Quaternion operator*(const Quaternion &, const Quaternion &);
    Quaternion();
    Quaternion(double, double, double, double);
    Quaternion(const Euler &e);
    Quaternion(const Quaternion &) = default;
    Quaternion(Quaternion &&) noexcept = default;
    Quaternion& operator=(const Quaternion &) = default;
    Quaternion& operator=(Quaternion &&) noexcept = default;
    ~Quaternion();

    double r() const { return m_r; }
    double i() const { return m_i; }
    double j() const { return m_j; }
    double k() const { return m_k; }

    Quaternion& inv();

    MatrixDbl toRotationMatrix() const;
  private:
    double m_r;
    double m_i;
    double m_j;
    double m_k;
  };  // class Quaternion
} // namespace dnpsoup


#endif
