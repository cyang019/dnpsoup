namespace dnpsoup {
  template<typename T1, typename T2>
  template<typename U1, typename U2,
           std::enable_if_t<is_frame_type<U1>::value, int>,
           std::enable_if_t<is_frame_type<U2>::value, int> >
  DipolarInteraction<T1, T2>::DipolarInteraction(double g1, double g2, size_t n1, size_t n2)
    : m_gamma1(g1), m_gamma2(g2), m_n1(n1), m_n2(n2),
    m_nbefore(0), m_nbetween(0), m_nafter(0)
  {
    if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2))
        - 0.25 * (kron(spin<P>(n1), spin<M>(n2)) + kron(spin<M>(n1) + spin<P>(n2)));
      m_a21 = -1.5 * (kron(spin<Z>(n1), spin<P>(n2)) + kron(spin<P>(n1), spin<Z>(n2)));
      m_a2n1 = -1.5 * (kron(spin<Z>(n1), spin<M>(n2)) + kron(spin<M>(n1), spin<Z>(n2)));
      m_a22 = -0.75 * kron(spin<P>(n1), spin<P>(n2));
      m_a2n2 = -0.75 * kron(spin<M>(n1), spin<M>(n2));
    } else if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, RotatingFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      m_a21 = -1.5 * kron(spin<P>(n1), spin<Z>(n2));
      m_a2n1 = -1.5 * kron(spin<Z>(n1), spin<M>(n2));
      m_a22 = MatrixCxDbl(0,0);
      m_a2n2 = MatrixCxDbl(0,0);
    } else if constexpr(std::is_same<T1, RotatingFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      m_a21 = -1.5 * kron(spin<Z>(n1), spin<P>(n2));
      m_a2n1 = -1.5 * kron(spin<Z>(n1), spin<M>(n2));
      m_a22 = MatrixCxDbl(0,0);
      m_a2n2 = MatrixCxDbl(0,0);
    } else { // Rotate & Rotate
      // homonuclear dipolar secular approx.
      if(approxEqual(m_gamma1, m_gamma2, eps)){ 
        m_a20 = kron(spin<Z>(n1), spin<Z>(n2))
        - 0.25 * (kron(spin<P>(n1), spin<M>(n2)) + kron(spin<M>(n1) + spin<P>(n2)));
      } else {
        // heteronuclear dipolar secular approx.
        m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      }
      m_a21 = MatrixCxDbl(0,0);
      m_a2n1 = MatrixCxDbl(0,0);
      m_a22 = MatrixCxDbl(0,0);
      m_a2n2 = MatrixCxDbl(0,0);
    }
  }

  template<typename T1, typename T2>
  template<typename U1, typename U2,
           std::enable_if_t<is_frame_type<U1>::value, int>,
           std::enable_if_t<is_frame_type<U2>::value, int> >
  DipolarInteraction<T1, T2>::DipolarInteraction(double g1, double g2, 
      size_t n1, size_t n2, size_t nbefore, size_t nbetween, size_t nafter)
    : m_gamma1(g1), m_gamma2(g2), m_n1(n1), m_n2(n2),
    m_nbefore(nbefore), m_nbetween(nbetween), m_nafter(nafter)
  {
    // 1. disregard nbefore & nafter to do arithmetics
    // 2. expand the matrix
    eye_1 = identity<cxdbl>(nbefore);
    eye_2 = identity<cxdbl>(nbetween);
    eye_3 = identity<cxdbl>(nafter);
    MatrixCxDbl izsz = 
      kron( kron(spin<Z>(n1), eye_2), spin<Z>(n2));
    if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      MatrixCxDbl a20_pre = izsz - 0.25 * (
          kron( kron(spin<P>(n1), eye_2), spin<M>(n2))
        + kron( kron(spin<M>(n1), eye_2), spin<P>(n2)));
      MatrixCxDbl a21_pre = -1.5 * (
          kron( kron(spin<Z>(n1), eye_2), spin<P>(n2))
        + kron( kron(spin<P>(n1), eye_2), spin<Z>(n2)));
      MatrixCxDbl a2n1_pre = -1.5 * (
          kron( kron(spin<Z>(n1), eye_2), spin<M>(n2))
        + kron( kron(spin<M>(n1), eye_2), spin<Z>(n2)));
      MatrixCxDbl a22_pre = -0.75 * kron( kron(spin<P>(n1), eye_2), spin<P>(n2));
      MatrixCxDbl a2n2_pre = -0.75 * kron( kron(spin<M>(n1), eye_2), spin<M>(n2));
      m_a20 = kron(kron(eye_1, a20_pre), eye_3);
      m_a21 = kron(kron(eye_1, a21_pre), eye_3);
      m_a2n1 = kron(kron(eye_1, a2n1_pre), eye_3);
      m_a22 = kron(kron(eye_1, a22_pre), eye_3);
      m_a2n2 = kron(kron(eye_1, a2n2_pre), eye_3);
    } else if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, RotatingFrame>::value) {
      MatrixCxDbl a20_pre = izsz;
      MatrixCxDbl a21_pre = -1.5 * kron( kron(spin<P>(n1), eye_2), spin<Z>(n2));
      MatrixCxDbl a2n1_pre = -1.5 * kron( kron(spin<M>(n1), eye_2), spin<Z>(n2));
      m_a20 = kron(kron(eye_1, a20_pre), eye_3);
      m_a21 = kron(kron(eye_1, a21_pre), eye_3);
      m_a2n1 = kron(kron(eye_1, a2n1_pre), eye_3);
      m_a22 = MatrixCxDbl(0,0);
      m_a2n2 = MatrixCxDbl(0,0);
    }
} // namespace dnpsoup

