#include "dnpsoup_core/spinsys/SpinId.h"
#include <utility>

using namespace std;

namespace dnpsoup {
  SpinId::SpinId(int id_val)
    : m_id(id_val)
  {}

  SpinId::~SpinId()
  {}


  bool SpinId::operator==(const SpinId &rhs) const
  {
    return rhs.m_id == m_id;
  }

  bool SpinId::operator!=(const SpinId &rhs) const
  {
    return rhs.m_id != m_id;
  }

  int SpinId::get() const 
  { return m_id; }

  SpinId& SpinId::set(int value)
  { m_id = value; return *this; }

  std::size_t SpinIdHash::operator()(const SpinId &s) const
  {
    return std::hash<int>{}(s.get());
  }
} // namespace dnpsoup
