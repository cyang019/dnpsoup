#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/common.h"


namespace dnpsoup {
  ObservableId::ObservableId(const SpinId &id1, const SpinId &id2)
  {
    m_id = genUniqueInt(
        static_cast<std::int64_t>(id1.get()),
        static_cast<std::int64_t>(id2.get()));
  }
} // namespace dnpsoup
