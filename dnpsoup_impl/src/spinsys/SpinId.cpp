#include "dnpsoup_core/spinsys/SpinId.h"
#include <utility>

using namespace std;

namespace dnpsoup {
  SpinId::SpinId(int id_val)
    : m_id(id_val)
  {}

  SpinId::~SpinId()
  {}


  bool SpinId::operator==(const SpinId &rhs)
  {
    return rhs.m_id == m_id;
  }

  int SpinId::get() const 
  { return m_id; }

  SpinId& SpinId::set(int value)
  { m_id = value; return *this; }

  std::size_t SpinIdHash::operator()(const SpinId &id_val) const
  {
    return std::hash<int>()(id_val.get());
  }
} // namespace dnpsoup
