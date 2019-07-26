#ifndef DNPSOUP_OBSERVABLE_H
#define DNPSOUP_OBSERVABLE_H

#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spin_physics_components/Rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include <vector>
#include <cstdint>
#include <string>
#include <hash>


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

  class Observable {
  public:
    Observable(const InteractionType &, const SpinId &);
    Observable(const InteractionType &, const SpinId &, const SpinId &);
    Observable(const Observable &) = default;
    Observable(Observable &&) noexcept = default;
    Observable& operator=(const Observable&) = default;
    Observable& operator=(Observable &&) noexcept = default;
    ~Observable() {}

    std::vector<SpinId> getSpinIds() const { return m_spin_ids; }
    InteractionType getType() const { return m_type; }

    double T1() const { return m_t1; }
    double T2() const { return m_t2; }
    Observable& setT1(double t1) { m_t1 = t1; return *this; }
    Observable& setT2(double t2) { m_t2 = t2; return *this; }

    Euler<> getEuler() const { return m_e; }
    Observable& setEuler(const Euler<> &e) const { m_e = e; return *this; }

    Property getProperty() const { return m_p; }
    Observable& setProperty(const Property &p) const { m_p = p; return *this; }
  private:
    InteractionType m_type;
    std::vector<SpinId> m_spin_ids;
    double m_t1;    // T1 relaxation, default infinity
    double m_t2;    // T2 relaxation, default infinity
    Euler<> m_e;
    Property m_p;
  };

  class ObservableIdHash {
    std::size_t operator(const ObservableId &o_id) const
    {
      return std::hash<std::string>{}(o_id.get());
    }
  };
} // namespace dnpsoup

#endif
