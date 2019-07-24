namespace dnpsoup {
  template<typename T>
  Quaternion::Quaternion(const Euler<T> &e)
  {
    static_assert(is_rotation_type<T>::value, 
        "Euler can either be ActiveRotation or PassiveRotation");
    if constexpr (std::is_same<T, ActiveRotation>::value){
      m_r = cos(0.5 * e.beta()) * cos(0.5 * (e.alpha() + e.gamma()));
      m_i = -sin(0.5 * e.beta()) * sin(0.5 * (e.alpha() - e.gamma()));
      m_j = sin(0.5 * e.beta()) * cos(0.5 * (e.alpha() - e.gamma()));
      m_k = cos(0.5 * e.beta()) * sin(0.5 * (e.alpha() + e.gamma()));
    } else {
      m_r = cos(0.5 * e.beta()) * cos(0.5 * (e.alpha() + e.gamma()));
      m_i = sin(0.5 * e.beta()) * sin(0.5 * (e.alpha() - e.gamma()));
      m_j = -sin(0.5 * e.beta()) * cos(0.5 * (e.alpha() - e.gamma()));
      m_k = -cos(0.5 * e.beta()) * sin(0.5 * (e.alpha() + e.gamma()));
    }
  }
} // namespace dnpsoup
