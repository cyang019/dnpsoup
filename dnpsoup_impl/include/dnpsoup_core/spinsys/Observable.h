#ifndef DNPSOUP_OBSERVABLE_H
#define DNPSOUP_OBSERVABLE_H

#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/interaction_types.h"
#include <vector>
#include <cstdint>
#include <string>
#include <functional>     // hash


namespace dnpsoup {
  class ObservableId {
  public:
    ObservableId(const InteractionType&, const SpinId &id1);
    ObservableId(const InteractionType&, const SpinId &id1, const SpinId &id2);
    ObservableId(const ObservableId &) = default;
    ObservableId(ObservableId &&) noexcept = default;
    ObservableId& operator=(const ObservableId &) = default;
    ObservableId& operator=(ObservableId &&) noexcept = default;
    ~ObservableId() {}

    std::string get() const { return m_id; }
  private:
    std::string m_id; 
  };

  inline bool operator==(const ObservableId &oid1, const ObservableId &oid2)
  { return oid1.get() == oid2.get(); }

  inline bool operator<(const ObservableId &oid1, const ObservableId &oid2)
  { return oid1.get() < oid2.get(); }

  inline bool operator>(const ObservableId &oid1, const ObservableId &oid2)
  { return oid1.get() > oid2.get(); }

  class Observable {
  public:
    Observable(const InteractionType &, const SpinId &);
    Observable(const InteractionType &, const SpinId &, const SpinId &);
    Observable(const Observable &) = default;
    Observable(Observable &&) noexcept;
    Observable& operator=(const Observable&) = default;
    Observable& operator=(Observable &&) noexcept;
    ~Observable() {}

    std::vector<SpinId> getSpinIds() const { return m_spin_ids; }
    InteractionType getType() const { return m_type; }

    double T1() const { return m_t1; }
    double T2() const { return m_t2; }
    Observable& setT1(double t1) 
    { m_t1 = t1; return *this; }

    Observable& setT2(double t2) 
    { m_t2 = t2; return *this; }

    Euler<> getEuler() const { return m_e; }
    Observable& setEuler(const Euler<> &e) 
    { m_e = e; return *this; }

    Property getProperty() const { return m_p; }
    Observable& setProperty(const Property &p) 
    { m_p = p; return *this; }
  private:
    InteractionType m_type;
    std::vector<SpinId> m_spin_ids;

    Euler<> m_e;
    Property m_p;
    double m_t1;    // T1 relaxation, default infinity
    double m_t2;    // T2 relaxation, default infinity
  };

  class ObservableIdHash {
  public:
    std::size_t operator()(const ObservableId &o_id) const
    {
      return std::hash<std::string>{}(o_id.get());
    }
  };
} // namespace dnpsoup

#endif
