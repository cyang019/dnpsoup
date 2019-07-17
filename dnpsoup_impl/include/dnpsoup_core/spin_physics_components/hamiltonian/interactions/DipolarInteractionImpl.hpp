namespace dnpsoup {
  template<FrameType T1, FrameType T2>
  DipolarInteraction<T1, T2>::DipolarInteraction(size_t n1, size_t n2)
    : m_n1(n1), m_n2(n2), m_nbefore(0), m_nbetween(0), m_nafter(0)
  {
    if constexpr(std::is_same<T1, FrameType::Lab>::value
        && std::is_same<T2, FrameType::Lab>::value) {
    } else if constexpr(std::is_same<T1, FrameType::Lab>::value
        && std::is_same<T2, FrameType::Rotate>::value) {
    } else if constexpr(std::is_same<T1, FrameType::Rotate>::value
        && std::is_same<T2, FrameType::Lab>::value) {
    } else { // Rotate & Rotate
    }
  }
} // namespace dnpsoup

