#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/common.h"
#include <string>


namespace dnpsoup {
  ObservableId::ObservableId(const InteractionType &t, const SpinId &id1, const SpinId &id2)
  {
    if(id1.get() < id2.get()){
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

  Observable(const InteractionType &t, const SpinId &id1, const SpinId &id2)
    : m_type(t), m_spin_ids({id1, id2})
  {}
} // namespace dnpsoup
