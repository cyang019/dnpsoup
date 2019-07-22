#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/Quaternion.h"
#include <cmath>

using namespace std;

namespace dnpsoup {
  Euler::Euler(const Quaternion &q)
  {
    double mat02 = 2.0 * (q.i() * q.k() + q.r() * q.j());
    double mat12 = 2.0 * (q.j() * q.k() - q.r() * q.i());
    double mat20 = 2.0 * (q.i() * q.k() - q.r() * q.j());
    double mat21 = 2.0 * (q.j() * q.k() + q.r() * q.i());
    double mat22 = q.r() * q.r() - q.i() * q.i() - q.j() * q.j() + q.k() * q.k();

    m_alpha = dnpsoup::atan(mat02,mat12);
    m_beta = dnpsoup::atan(sqrt(1.0 - mat22),mat22);
    m_gamma = dnpsoup::atan(mat21,-mat20);
  }

  Euler operator*(const Euler &e1, const Euler &e2)
  {
    Quaternion q1(e1);
    Quaternion q2(e2);
    Quaternion res = q1 * q2;
    return Euler(res);
  }
} // namespace dnpsoup
