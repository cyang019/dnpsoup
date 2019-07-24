namespace dnpsoup {
  template<typename T>
  Euler<T>::Euler(const Quaternion &q)
  {
    static_assert(is_rotation_type<T>::value, 
        "Euler can either be ActiveRotation or PassiveRotation");

    double mat02 = 2.0 * (q.i() * q.k() + q.r() * q.j());
    double mat12 = 2.0 * (q.j() * q.k() - q.r() * q.i());
    double mat20 = 2.0 * (q.i() * q.k() - q.r() * q.j());
    double mat21 = 2.0 * (q.j() * q.k() + q.r() * q.i());
    double mat22 = q.r() * q.r() - q.i() * q.i() - q.j() * q.j() + q.k() * q.k();

    if constexpr(std::is_same<T, ActiveRotation>::value){
      m_alpha = dnpsoup::atan(mat02,mat12);
      m_beta = dnpsoup::atan(sqrt(1.0 - mat22),mat22);
      m_gamma = dnpsoup::atan(mat21,-mat20);
    } else {
      m_gamma = -dnpsoup::atan(mat02,mat12);
      m_beta = -dnpsoup::atan(sqrt(1.0 - mat22),mat22);
      m_alpha = -dnpsoup::atan(mat21,-mat20);
    }
  }

  template<typename T>
  inline
  Euler<T> operator*(const Euler<T> &e1, const Euler<T> &e2)
  {
    Quaternion q1(e1);
    Quaternion q2(e2);
    Quaternion res = q1 * q2;
    return Euler<T>(res);
  }
} // namespace dnpsoup

