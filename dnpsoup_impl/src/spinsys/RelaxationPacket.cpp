#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/constants.h"
#include <limits>

using namespace std;


namespace dnpsoup {
    RelaxationPacket::RelaxationPacket(
        const SpinId &sid, const SpinEntity &sinfo,
        std::size_t nbefore, std::size_t nafter)
      : m_id(sid), m_t1(sinfo.getT1()), m_t2(sinfo.getT2())
    {
      auto n = getMatrixDimension(sinfo.getSpinType());
      m_x = expandMatrix(spin<X>(n), nbefore, nafter);
      m_y = expandMatrix(spin<Y>(n), nbefore, nafter);
      m_z = expandMatrix(spin<Z>(n), nbefore, nafter);
    }

    MatrixCxDbl RelaxationPacket::genSuperOpT1(const MatrixCxDbl &eigenvec) const
    {
      double t1_inv = 0.0;
      if(m_t1 < numeric_limits<double>::max()){
        t1_inv = 1.0/m_t1;
      }
      return t1SuperOp(t1_inv, eigenvec, m_x, m_y);
    }

    MatrixCxDbl RelaxationPacket::genSuperOpT2(const MatrixCxDbl &eigenvec) const
    {
      if(m_t2 >= 2 * m_t1 - eps){
        throw RelaxationValueError("T2 cannot be larger than 2T1");
      }
      double t2_prime_inv = 0;
      if(m_t2 < numeric_limits<double>::max() && m_t1 < numeric_limits<double>::max()){
        t2_prime_inv = 1.0/m_t2 - 0.5/m_t1;
      }
      return t2SuperOp(t2_prime_inv, eigenvec, m_z);
    }
} // namespace dnpsoup
