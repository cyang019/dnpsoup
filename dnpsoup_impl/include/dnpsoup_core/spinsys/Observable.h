#ifndef DNPSOUP_OBSERVABLE_H
#define DNPSOUP_OBSERVABLE_H

#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spin_physics_components/Rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include <cstdint>

namespace dnpsoup {
  class ObservableId {
  public:
    ObservableId(const SpinId &id1, const SpinId &id2);
    ObservableId(const ObservableId &) = default;
    ObservableId(ObservableId &&) noexcept = default;
    ObservableId& operator=(const ObservableId &) = default;
    ObservableId& operator=(ObservableId &&) noexcept = default;
    ~ObservableId() {}

    std::int64_t get() const { return m_id; }
  private:
    std::int64_t m_id; 
  };

  class Observable {
  public:
    
  private:
    double m_t1;
    double m_t2;
    Euler<> m_e;
    Property m_p;
  };
} // namespace dnpsoup

#endif
