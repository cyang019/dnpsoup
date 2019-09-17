#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/common.h"
#include <string>
#include <algorithm>


namespace dnpsoup {
  ObservableId::ObservableId(const InteractionType &t, const SpinId &id1)
  {
    m_id = interactionName(t) + "_" + std::to_string(id1.get()); 
  }

  ObservableId::ObservableId(const InteractionType &t, const SpinId &id1, const SpinId &id2)
  {
    if(id1 < id2){
      m_id = interactionName(t) + "_" 
        + std::to_string(id1.get()) 
        + "_" + std::to_string(id2.get());
    } else {
      m_id = interactionName(t) + "_" 
        + std::to_string(id2.get()) 
        + "_" + std::to_string(id1.get());
    }
  }

  ObservableId::ObservableId(
      const InteractionType &i_type, const SpinType &s_type)
  {
    m_id = interactionName(i_type) + "_" + toString(s_type);
  }

  // class Observable
  Observable::Observable()
    : m_type(InteractionType::Null)
  {}

  Observable::Observable(const InteractionType &t, const SpinId &id)
    : m_type(t), m_spin_ids({id})
  {}

  Observable::Observable(const InteractionType &t, const SpinId &id1, const SpinId &id2)
    : m_type(t), m_spin_ids({id1, id2})
  {}

  Observable::Observable(
      const InteractionType &t, const std::vector<SpinId> &spin_ids)
    : m_type(t), m_spin_ids(spin_ids)
  {}

  Observable::Observable(Observable &&rhs) noexcept
    : m_type(std::move(rhs.m_type)), m_spin_ids(std::move(rhs.m_spin_ids)),
    m_e(std::move(rhs.m_e)), m_p(std::move(rhs.m_p))
  {}

  Observable& Observable::operator=(Observable &&rhs) noexcept
  {
    m_type = std::move(rhs.m_type);
    m_spin_ids = std::move(rhs.m_spin_ids);
    m_e = std::move(rhs.m_e);
    m_p = std::move(rhs.m_p);
    return *this;
  }

  Observable& Observable::setPropertyValue(const ValueName &vname, double val)
  {
    this->m_p.set(vname, val);
    return *this;
  }

  double Observable::getPropertyValue(const ValueName &vname) const
  {
    return this->m_p.get(vname);
  }

  bool Observable::hasSpinId(const SpinId &sid) const
  {
    if(std::find(m_spin_ids.begin(), m_spin_ids.end(), sid) == m_spin_ids.end()){
      return false;
    }
    return true;
  }
} // namespace dnpsoup
