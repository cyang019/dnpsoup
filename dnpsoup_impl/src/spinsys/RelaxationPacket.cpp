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

    MatrixCxDbl RelaxationPacket::genSuperOpT1() const
    {
      double t1_inv = 0.0;
      if(m_t1 < numeric_limits<double>::max()){
        t1_inv = 1.0/m_t1;
      }
      return t1SuperOp(t1_inv, m_x, m_y);
    }

    MatrixCxDbl RelaxationPacket::genSuperOpT2() const
    {
      if(m_t2 >= 2 * m_t1 - eps){
        throw RelaxationValueError("T2 cannot be larger than 2T1");
      }
      double t2_prime_inv = 0;
      if(m_t2 < numeric_limits<double>::max() && m_t1 < numeric_limits<double>::max()){
        t2_prime_inv = 1.0/m_t2 - 0.5/m_t1;
      }
      return t2SuperOp(t2_prime_inv, m_z);
    }

    CustomRelaxationPacket::CustomRelaxationPacket(
        const vector<pair<vector<pair<SpinType, OperatorType>>, double>> &ops,
        double t)
      : m_t(t)
    {
      vector<MatrixCxDbl> mats;
      for(const auto &[vec, scale] : ops) {
        vector<MatrixCxDbl> sub_mats;
        for(const auto &[s_type, op_type] : vec) {
          const auto sz = getMatrixDimension(s_type);
          switch(op_type) {
            case OperatorType::Identity:
              sub_mats.emplace_back(spin<OperatorType::Identity>(sz));
              break;
            case OperatorType::Minus:
              sub_mats.emplace_back(spin<OperatorType::Minus>(sz));
              break;
            case OperatorType::Plus:
              sub_mats.emplace_back(spin<OperatorType::Plus>(sz));
              break;
            case OperatorType::X:
              sub_mats.emplace_back(spin<OperatorType::X>(sz));
              break;
            case OperatorType::Y:
              sub_mats.emplace_back(spin<OperatorType::Y>(sz));
              break;
            case OperatorType::Z:
              sub_mats.emplace_back(spin<OperatorType::Z>(sz));
              break;
          }
        }
        MatrixCxDbl sub_mat = scale * kron(sub_mats);
        mats.push_back(sub_mat);
      }
      if(mats.size() > 0) {
        m_mat = mats[0];
        for(size_t i = 1; i < mats.size(); ++i) {
          m_mat += mats[i];
        }
      }
    }

    MatrixCxDbl CustomRelaxationPacket::genSuperOp() const
    {
      double t_inv = 0.0;
      if(m_t < numeric_limits<double>::max()){
        t_inv = 1.0/m_t;
      }
      return t_inv * secularRelaxationSuperOp(m_mat);
    }

    MatrixCxDbl RelaxationPacketCollection::genMatrix() const
    {
      if(m_rpackets.size() == 0) return MatrixCxDbl(0);
      auto t1_superop = m_rpackets[0].genSuperOpT1();
      auto t2_superop = m_rpackets[0].genSuperOpT2();
      for(size_t i = 1; i < m_rpackets.size(); ++i) {
        t1_superop += m_rpackets[i].genSuperOpT1();
        t2_superop += m_rpackets[i].genSuperOpT2();
      }
      if(m_crpackets.size() == 0) {
        return t1_superop + t2_superop;
      }
      auto t_custom_superop = m_crpackets[0].genSuperOp();
      for(size_t i = 1; i < m_crpackets.size(); ++i) {
        t_custom_superop += m_crpackets[i].genSuperOp();
      }
      return t1_superop + t2_superop + t_custom_superop;
    }

    void RelaxationPacketCollection::addRelaxationPacket(
        const SpinId &sid, const SpinEntity &sinfo,
        std::size_t nbefore, std::size_t nafter)
    {
      m_rpackets.emplace_back(sid, sinfo, nbefore, nafter);
    }

    void RelaxationPacketCollection::addCustomRelaxationPacket(
        const vector<pair<vector<pair<SpinType, OperatorType>>, double>> &ops,
        double t)
    {
      m_crpackets.emplace_back(ops, t);
    }
} // namespace dnpsoup
