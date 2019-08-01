#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/common.h"
#include <string>


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

  // class Observable
  Observable::Observable(const InteractionType &t, const SpinId &id)
    : m_type(t), m_spin_ids({id})
  {}

  Observable::Observable(const InteractionType &t, const SpinId &id1, const SpinId &id2)
    : m_type(t), m_spin_ids({id1, id2})
  {}

  Observable::Observable(Observable &&rhs) noexcept
    : m_type(std::move(rhs.m_type)), m_spin_ids(std::move(rhs.m_spin_ids)),
    m_e(std::move(rhs.m_e)), m_p(std::move(rhs.m_p)),
    m_t1(std::move(rhs.m_t1)), m_t2(std::move(rhs.m_t1))
  {}

  Observable& Observable::operator=(Observable &&rhs) noexcept
  {
    m_type = std::move(rhs.m_type);
    m_spin_ids = std::move(rhs.m_spin_ids);
    m_e = std::move(rhs.m_e);
    m_p = std::move(rhs.m_p);
    m_t1 = std::move(rhs.m_t1);
    m_t2 = std::move(rhs.m_t2);
    return *this;
  }
} // namespace dnpsoup
