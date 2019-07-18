namespace dnpsoup {

  template<typename T>
  template<typename U,
           std::enable_if_t<is_frame_type<U>::value, int>>
  ChemicalShiftInteraction<T>::ChemicalShiftInteraction(double gamma, size_t n)
    : InteractionInterface(), m_gamma(gamma), m_n(n), m_nbefore(0), m_nafter(0)
  { 
    if constexpr(std::is_same<T, LabFrame>::value){
      m_ix = spin<X>(m_n);
      m_iy = spin<Y>(m_n);
      m_iz = spin<Z>(m_n);
    } else {  // Rotating Frame
      m_iz = spin<Z>(m_n);
    }
  }

  template<typename T>
  template<typename U,
           std::enable_if_t<is_frame_type<U>::value, int>>
  ChemicalShiftInteraction<T>::ChemicalShiftInteraction(double gamma, size_t n, size_t nbefore, size_t nafter)
    : InteractionInterface(), m_gamma(gamma), m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  { 
    if constexpr(std::is_same<T, LabFrame>::value){
      m_ix = kron(kron(identity<cxdbl>(nbefore), spin<X>(n)), identity<cxdbl>(nafter));
      m_iy = kron(kron(identity<cxdbl>(nbefore), spin<Y>(n)), identity<cxdbl>(nafter));
      m_iz = kron(kron(identity<cxdbl>(nbefore), spin<Z>(n)), identity<cxdbl>(nafter));
    } else {  // Rotating Frame
      m_iz = kron(kron(identity<cxdbl>(nbefore), spin<Z>(n)), identity<cxdbl>(nafter));
    }
  }

  template<typename T>
  template<typename U,
           std::enable_if_t<is_frame_type<U>::value, int>>
  MatrixCxDbl ChemicalShiftInteraction<T>::genMatrix(
      const Property *ptr_csa,
      const Euler &e) const
  {
    const double szz = ptr_csa->get(ValueName::zz);

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
      const double sxx = ptr_csa->get(ValueName::xx);
      const double syy = ptr_csa->get(ValueName::yy);

      const double bz = ptr_csa->get(ValueName::bz);
      double coeff1 = (sxx - syy) * 0.5 * sa * sb * s2g + ca * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
      double coeff2 = (syy - sxx) * 0.5 * ca * sb * s2g + sa * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
      MatrixCxDbl res = m_gamma * bz * m_iz + coeff1 * m_ix + coeff2 * m_iy + coeff3 * m_iz;
      return res;
    } else {  // Rotating Frame
      MatrixCxDbl res = coeff3 * m_iz;
      return res;
    }
  }

  template<typename T>
  template<typename U,
           std::enable_if_t<is_frame_type<U>::value, int>>
  size_t ChemicalShiftInteraction<T>::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }

} // namespace dnpsoup
