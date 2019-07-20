namespace dnpsoup {

  template<typename T>
  ChemicalShiftInteraction<T>::ChemicalShiftInteraction(double gamma, size_t n)
    : InteractionInterface(), m_gamma(gamma), m_n(n), m_nbefore(0), m_nafter(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    if constexpr(std::is_same<T, LabFrame>::value){
      m_x = spin<X>(m_n);
      m_y = spin<Y>(m_n);
      m_z = spin<Z>(m_n);
    } else {  // Rotating Frame
      m_z = spin<Z>(m_n);
    }
  }

  template<typename T>
  ChemicalShiftInteraction<T>::ChemicalShiftInteraction(double gamma, size_t n, size_t nbefore, size_t nafter)
    : InteractionInterface(), m_gamma(gamma), m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    if constexpr(std::is_same<T, LabFrame>::value){
      m_x = kron(kron(identity<cxdbl>(nbefore), spin<X>(n)), identity<cxdbl>(nafter));
      m_y = kron(kron(identity<cxdbl>(nbefore), spin<Y>(n)), identity<cxdbl>(nafter));
      m_z = kron(kron(identity<cxdbl>(nbefore), spin<Z>(n)), identity<cxdbl>(nafter));
    } else {  // Rotating Frame
      m_z = kron(kron(identity<cxdbl>(nbefore), spin<Z>(n)), identity<cxdbl>(nafter));
    }
  }

  template<typename T>
  MatrixCxDbl ChemicalShiftInteraction<T>::genMatrix(
      const Property &csa,
      const Euler &e) const
  {
    const double szz = csa.get(ValueName::zz);
    const double sxx = csa.get(ValueName::xx);
    const double syy = csa.get(ValueName::yy);

    const double sb = sin(e.beta());
    const double cb = cos(e.beta());
    const double sg = sin(e.gamma());
    const double cg = cos(e.gamma());
    double coeff3 = szz * cb * cb + sb * sb * (sxx * cg * cg + syy * sg * sg);

    if constexpr(std::is_same<T, LabFrame>::value){
      const double sa = sin(e.alpha());
      const double ca = cos(e.alpha());
      const double s2b = sin(2.0*e.beta());
      const double s2g = sin(2.0*e.gamma());

      const double bz = csa.get(ValueName::bz);
      double coeff1 = (sxx - syy) * 0.5 * sa * sb * s2g + ca * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
      double coeff2 = (syy - sxx) * 0.5 * ca * sb * s2g + sa * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
      MatrixCxDbl res = m_gamma * bz * m_z + coeff1 * m_x + coeff2 * m_y + coeff3 * m_z;
      return res;
    } else {  // Rotating Frame
      MatrixCxDbl res = coeff3 * m_z;
      return res;
    }
  }

  template<typename T>
  size_t ChemicalShiftInteraction<T>::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }

} // namespace dnpsoup
