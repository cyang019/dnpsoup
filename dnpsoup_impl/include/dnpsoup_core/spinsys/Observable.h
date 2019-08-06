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
    ObservableId(const InteractionType&, const SpinType &t);
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
    Observable(const InteractionType&, const std::vector<SpinId> &sids); 
    Observable(const Observable &) = default;
    Observable(Observable &&) noexcept;
    Observable& operator=(const Observable&) = default;
    Observable& operator=(Observable &&) noexcept;
    ~Observable() {}

    const std::vector<SpinId>& getSpinIds() const { return m_spin_ids; }
    const InteractionType& getType() const { return m_type; }

    bool hasSpinId(const SpinId &) const;

    const Euler<>& getEuler() const { return m_e; }
    Observable& setEuler(const Euler<> &e) 
    { m_e = e; return *this; }

    const Property& getProperty() const { return m_p; }
    Observable& setProperty(const Property &p) 
    { m_p = p; return *this; }

    Observable& setPropertyValue(const ValueName &vname, double val);
    double getPropertyValue(const ValueName &vname) const;
  private:
    InteractionType m_type;
    std::vector<SpinId> m_spin_ids;

    Euler<> m_e;
    Property m_p;
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
