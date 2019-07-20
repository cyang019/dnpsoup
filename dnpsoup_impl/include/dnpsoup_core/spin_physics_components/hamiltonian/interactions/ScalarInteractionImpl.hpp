namespace dnpsoup {
  template<typename T1, typename T2>
  ScalarInteraction<T1, T2>::ScalarInteraction(double g1, double g2, size_t n1, size_t n2)
    : m_n1(n1), m_n2(n2),  
    m_nbefore(0), m_nbetween(0), m_nafter(0),
    m_gamma1(g1), m_gamma2(g2)
  {
    static_assert(is_frame_type<T1>::value && is_frame_type<T2>::value,
        "T1 and T2 need to be either LabFrame or RotatingFrame respectively.");
    m_zz = kron(spin<Z>(n1), spin<Z>(n2));
  }

  template<typename T1, typename T2>
  ScalarInteraction<T1, T2>::ScalarInteraction(double g1, double g2, 
      size_t n1, size_t n2, size_t nbefore, size_t nbetween, size_t nafter)
    : m_n1(n1), m_n2(n2),  
    m_nbefore(nbefore), m_nbetween(nbetween), m_nafter(nafter),
    m_gamma1(g1), m_gamma2(g2)
  {
    static_assert(is_frame_type<T1>::value && is_frame_type<T2>::value,
        "T1 and T2 need to be either LabFrame or RotatingFrame respectively.");
    // 1. disregard nbefore & nafter to do arithmetics
    // 2. expand the matrix
    const MatrixCxDbl eye_1 = identity<cxdbl>(nbefore);
    const MatrixCxDbl eye_2 = identity<cxdbl>(nbetween);
    const MatrixCxDbl eye_3 = identity<cxdbl>(nafter);
    MatrixCxDbl zz = kron( kron(spin<Z>(n1), eye_2), spin<Z>(n2));
    m_zz = kron(kron(eye_1, zz), eye_3);
  }

  template<typename T1, typename T2>
  MatrixCxDbl ScalarInteraction<T1, T2>::genMatrix(
      const Property &p, const Euler &e) const
  {
    return p.get(ValueName::scalar) * m_zz;
  }
} // namespace dnpsoup


