#ifndef DNPSOUP_SPINID_H
#define DNPSOUP_SPINID_H

#include <utility>


namespace dnpsoup {
  class SpinId {
  public:
    SpinId(int);
    SpinId(const SpinId &) = default;
    SpinId(SpinId &&) noexcept = default;
    SpinId& operator=(const SpinId &) = default;
    SpinId& operator=(SpinId &&) noexcept = default;
    ~SpinId();

    bool operator==(const SpinId &);

    int get() const;
    SpinId& set(int);
  private:
    int m_id;
  };

  class SpinIdHash {
    std::size_t operator()(const SpinId &) const;
  };
} // namespace dnpsoup

#endif
