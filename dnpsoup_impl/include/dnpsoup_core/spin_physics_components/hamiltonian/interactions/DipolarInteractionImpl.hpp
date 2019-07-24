namespace dnpsoup {
  template<typename T1, typename T2>
  DipolarInteraction<T1, T2>::DipolarInteraction(double g1, double g2, size_t n1, size_t n2)
    : m_n1(n1), m_n2(n2),  
    m_nbefore(0), m_nbetween(0), m_nafter(0),
    m_gamma1(g1), m_gamma2(g2)
  {
    static_assert(is_frame_type<T1>::value && is_frame_type<T2>::value,
        "T1 and T2 need to be either LabFrame or RotatingFrame respectively.");
    if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2))
        - 0.25 * (kron(spin<P>(n1), spin<M>(n2)) + kron(spin<M>(n1), spin<P>(n2)));
      m_a21 = -1.5 * (kron(spin<Z>(n1), spin<P>(n2)) + kron(spin<P>(n1), spin<Z>(n2)));
      m_a2n1 = -1.5 * (kron(spin<Z>(n1), spin<M>(n2)) + kron(spin<M>(n1), spin<Z>(n2)));
      m_a22 = -0.75 * kron(spin<P>(n1), spin<P>(n2));
      m_a2n2 = -0.75 * kron(spin<M>(n1), spin<M>(n2));
    } else if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, RotatingFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      m_a21 = -1.5 * kron(spin<P>(n1), spin<Z>(n2));
      m_a2n1 = -1.5 * kron(spin<Z>(n1), spin<M>(n2));
      m_a22 = MatrixCxDbl();
      m_a2n2 = MatrixCxDbl();
    } else if constexpr(std::is_same<T1, RotatingFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      m_a21 = -1.5 * kron(spin<Z>(n1), spin<P>(n2));
      m_a2n1 = -1.5 * kron(spin<Z>(n1), spin<M>(n2));
      m_a22 = MatrixCxDbl();
      m_a2n2 = MatrixCxDbl();
    } else { // Rotate & Rotate
      // homonuclear dipolar secular approx.
      if(approxEqual(m_gamma1, m_gamma2, eps)){ 
        m_a20 = kron(spin<Z>(n1), spin<Z>(n2))
        - 0.25 * (kron(spin<P>(n1), spin<M>(n2)) + kron(spin<M>(n1), spin<P>(n2)));
      } else {
        // heteronuclear dipolar secular approx.
        m_a20 = kron(spin<Z>(n1), spin<Z>(n2));
      }
      m_a21 = MatrixCxDbl();
      m_a2n1 = MatrixCxDbl();
      m_a22 = MatrixCxDbl();
      m_a2n2 = MatrixCxDbl();
    }
  }

  template<typename T1, typename T2>
  DipolarInteraction<T1, T2>::DipolarInteraction(double g1, double g2, 
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
      m_a22 = MatrixCxDbl();
      m_a2n2 = MatrixCxDbl();
    } else if constexpr(std::is_same<T1, RotatingFrame>::value
        && std::is_same<T2, LabFrame>::value){
      MatrixCxDbl a20_pre = izsz;
      MatrixCxDbl a21_pre = -1.5 * kron( kron(spin<Z>(n1), eye_2), spin<P>(n2));
      MatrixCxDbl a2n1_pre = -1.5 * kron( kron(spin<Z>(n1), eye_2), spin<M>(n2));
      m_a20 = kron(kron(eye_1, a20_pre), eye_3);
      m_a21 = kron(kron(eye_1, a21_pre), eye_3);
      m_a2n1 = kron(kron(eye_1, a2n1_pre), eye_3);
      m_a22 = MatrixCxDbl();
      m_a2n2 = MatrixCxDbl();
    } else {  // Rotating & Rotating
      if (approxEqual(m_gamma1, m_gamma2, eps)){  // homonuclear dipole
        MatrixCxDbl a20_pre = izsz - 0.25 * (
            kron( kron(spin<P>(n1), eye_2), spin<M>(n2))
          + kron( kron(spin<M>(n1), eye_2), spin<P>(n2)));
        m_a20 = kron(kron(eye_1, a20_pre), eye_3);
        m_a21 = MatrixCxDbl();
        m_a2n1 = MatrixCxDbl();
        m_a22 = MatrixCxDbl();
        m_a2n2 = MatrixCxDbl();
      } else { // heteronuclear dipole
        MatrixCxDbl a20_pre = izsz;
        m_a20 = kron(kron(eye_1, a20_pre), eye_3);
        m_a21 = MatrixCxDbl();
        m_a2n1 = MatrixCxDbl();
        m_a22 = MatrixCxDbl();
        m_a2n2 = MatrixCxDbl();
      }
    }
  }

  template<typename T1, typename T2>
  template<typename R>
  MatrixCxDbl DipolarInteraction<T1, T2>::genMatrix(
      const Property &p, const Euler<R> &e) const
  {
    const double distance = p.get(ValueName::distance);
    const double d = 1.0e-7 * m_gamma1 * m_gamma2 * dnpsoup::h / (distance * distance * distance) * 1.0e30;
    if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, LabFrame>::value){
      if constexpr(std::is_same<R, ActiveRotation>::value){
        return d * (m_a20 * calcF20(e.beta()) 
            + m_a21 * calcF21(e.alpha(), e.beta())
            + m_a2n1 * calcF2n1(e.alpha(), e.beta())
            + m_a22 * calcF22(e.alpha(), e.beta())
            + m_a2n2 * calcF2n2(e.alpha(), e.beta()));
      } else {
        return d * (m_a20 * calcF20(-e.beta()) 
            + m_a21 * calcF21(-e.gamma(), -e.beta())
            + m_a2n1 * calcF2n1(-e.gamma(), -e.beta())
            + m_a22 * calcF22(-e.gamma(), -e.beta())
            + m_a2n2 * calcF2n2(-e.gamma(), -e.beta()));
      }
    } else if constexpr(std::is_same<T1, LabFrame>::value
        && std::is_same<T2, RotatingFrame>::value){
      if constexpr(std::is_same<R, ActiveRotation>::value){
        return d * (m_a20 * calcF20(e.beta()) 
            + m_a21 * calcF21(e.alpha(), e.beta())
            + m_a2n1 * calcF2n1(e.alpha(), e.beta()));
      } else {
        return d * (m_a20 * calcF20(-e.beta()) 
            + m_a21 * calcF21(-e.gamma(), -e.beta())
            + m_a2n1 * calcF2n1(-e.gamma(), -e.beta()));
      }
    } else if constexpr(std::is_same<T1, RotatingFrame>::value
        && std::is_same<T2, LabFrame>::value) {
      if constexpr(std::is_same<R, ActiveRotation>::value){
        return d * (m_a20 * calcF20(e.beta()) 
            + m_a21 * calcF21(e.alpha(), e.beta())
            + m_a2n1 * calcF2n1(e.alpha(), e.beta()));
      } else {
        return d * (m_a20 * calcF20(-e.beta()) 
            + m_a21 * calcF21(-e.gamma(), -e.beta())
            + m_a2n1 * calcF2n1(-e.gamma(), -e.beta()));
      }
    } else {  // Rotating Frame for Both
      if constexpr(std::is_same<R, ActiveRotation>::value){
        return d * (m_a20 * calcF20(e.beta()));
      } else {
      }
        return d * (m_a20 * calcF20(-e.beta()));
    }
  }
} // namespace dnpsoup

