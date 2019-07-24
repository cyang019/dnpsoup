#include "dnpsoup_core/spin_physics_components/rotation/Quaternion.h"
#include <cmath>

using namespace std;


namespace dnpsoup {
  Quaternion::Quaternion()
    : m_r(1.0), m_i(0.0), m_j(0.0), m_k(0.0)
  {}

  Quaternion::Quaternion(double r, double i, double j, double k)
    : m_r(r), m_i(i), m_j(j), m_k(k)
  {}

  Quaternion::~Quaternion() {}

  Quaternion& Quaternion::inv()
  {
    const double sq_norm = m_r * m_r + m_i * m_i + m_j * m_j + m_k * m_k;
    m_r /= sq_norm;
    m_i /= -sq_norm;
    m_j /= -sq_norm;
    m_k /= -sq_norm;
    return *this;
  }

  MatrixDbl Quaternion::toRotationMatrix() const
  {
    MatrixDbl mat(3,3);
    mat(0,0) = m_r * m_r + m_i * m_i - m_j * m_j - m_k * m_k;
    mat(0,1) = 2.0 * (m_i * m_j - m_r * m_k);
    mat(0,2) = 2.0 * (m_i * m_k + m_r * m_j);
    mat(1,0) = 2.0 * (m_i * m_j + m_r * m_k);
    mat(1,1) = m_r * m_r - m_i * m_i + m_j * m_j - m_k * m_k;
    mat(1,2) = 2.0 * (m_j * m_k - m_r * m_i);
    mat(2,0) = 2.0 * (m_i * m_k - m_r * m_j);
    mat(2,1) = 2.0 * (m_j * m_k + m_r * m_i);
    mat(2,2) = m_r * m_r - m_i * m_i - m_j * m_j + m_k * m_k;
    return mat;
  }

  // operators
  Quaternion operator+(const Quaternion &q1, const Quaternion &q2)
  {
    Quaternion q;
    q.m_r = q1.m_r + q2.m_r;
    q.m_i = q1.m_i + q2.m_i;
    q.m_j = q1.m_j + q2.m_j;
    q.m_k = q1.m_k + q2.m_k;
    return q;
  }

  Quaternion operator-(const Quaternion &q1, const Quaternion &q2)
  {
    Quaternion q;
    q.m_r = q1.m_r - q2.m_r;
    q.m_i = q1.m_i - q2.m_i;
    q.m_j = q1.m_j - q2.m_j;
    q.m_k = q1.m_k - q2.m_k;
    return q;
  }

  Quaternion operator*(const Quaternion &q1, const Quaternion &q2)
  {
    Quaternion q;
    q.m_r = q1.m_r * q2.m_r - q1.m_i * q2.m_i - q1.m_j * q2.m_j - q1.m_k * q2.m_k;
    q.m_i = q1.m_r * q2.m_i + q1.m_i * q2.m_r + q1.m_j * q2.m_k - q1.m_k * q2.m_j;
    q.m_j = q1.m_r * q2.m_j - q1.m_i * q2.m_k + q1.m_j * q2.m_r + q1.m_k * q2.m_i;
    q.m_k = q1.m_r * q2.m_k + q1.m_i * q2.m_j - q1.m_j * q2.m_i + q1.m_k * q2.m_r;
    return q;
  }
} // dnpsoup
