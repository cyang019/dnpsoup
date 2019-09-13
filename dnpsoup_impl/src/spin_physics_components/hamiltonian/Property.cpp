#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/errors.h"
#include <sstream>
#include <utility>

using namespace std;

namespace dnpsoup {
  std::ostream& operator<<(std::ostream &os, const ValueName &vname)
  {
    switch(vname){
      case ValueName::scalar:
        os << "scalar";
        break;
      case ValueName::distance:
        os << "distance";
        break;
      case ValueName::offset:
        os << "offset";
        break;
      case ValueName::xx:
        os << "xx";
        break;
      case ValueName::yy:
        os << "yy";
        break;
      case ValueName::zz:
        os << "zz";
        break;
      case ValueName::b0:
        os << "b0";
        break;
      case ValueName::freq:
        os << "freq";
        break;
      case ValueName::phase:
        os << "phase";
        break;
      default:
        break;
    }
    return os;
  }

  std::string toString(const ValueName &vname)
  {
    switch(vname){
      case ValueName::scalar:
        return "scalar";
        break;
      case ValueName::distance:
        return "distance";
        break;
      case ValueName::offset:
        return "offset";
        break;
      case ValueName::xx:
        return "xx";
        break;
      case ValueName::yy:
        return "yy";
        break;
      case ValueName::zz:
        return "zz";
        break;
      case ValueName::b0:
        return "b0";
        break;
      case ValueName::freq:
        return "freq";
        break;
      case ValueName::phase:
        return "phase";
        break;
      default:
        return "Unknown";
        break;
    }
    return "Unknown";
  }

  Property::Property()
  {
  }
  
  const double& Property::get(const ValueName &name) const
  {
#ifndef NDEBUG
    if(m_values.find(name) == m_values.end()){
      const std::string err_str = toString(name) + " not found in property.";
      throw NameNotFoundInProperty(err_str);
    }
#endif
    return m_values.at(name);
  }
  
  Property& Property::set(const ValueName &name, double value)
  {
    m_values[name] = value;
    return *this;
  }

  Property genCsaProperty(double b0, double sxx, double syy, double szz)
  {
    Property p;
    p.set(ValueName::xx, sxx);
    p.set(ValueName::yy, syy);
    p.set(ValueName::zz, szz);
    p.set(ValueName::b0, b0);
    return p;
  }

  // distances are in Anstrom
  Property genDipoleProperty(double distance)
  {
    Property p;
    p.set(ValueName::distance, distance);
    return p;
  }

  Property genScalarProperty(double val)
  {
    Property p;
    p.set(ValueName::scalar, val);
    return p;
  }

  Property genWaveProperty(double freq, double phase, double offset)
  {
    Property p;
    p.set(ValueName::freq, freq);
    p.set(ValueName::phase, phase);
    p.set(ValueName::offset, offset);
    return p;
  }
} // namespace dnpsoup
