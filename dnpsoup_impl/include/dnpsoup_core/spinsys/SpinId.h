#ifndef DNPSOUP_SPINID_H
#define DNPSOUP_SPINID_H

#include <utility>
#include <string>


namespace dnpsoup {
  class SpinId {
  public:
    SpinId(int);
    SpinId(const SpinId &) = default;
    SpinId(SpinId &&) noexcept = default;
    SpinId& operator=(const SpinId &) = default;
    SpinId& operator=(SpinId &&) noexcept = default;
    ~SpinId();

    bool operator==(const SpinId &) const;
    bool operator!=(const SpinId &) const;

    int get() const;
    SpinId& set(int);
  private:
    int m_id;
  };

  inline std::string toString(const SpinId &t){
    return std::to_string(t.get());
  }

  inline bool operator<(const SpinId &id1, const SpinId &id2)
  {
    return id1.get() < id2.get();
  }

  inline bool operator>(const SpinId &id1, const SpinId &id2)
  {
    return id1.get() > id2.get();
  }

  class SpinIdHash {
    std::size_t operator()(const SpinId &) const;
  };
} // namespace dnpsoup

#endif
