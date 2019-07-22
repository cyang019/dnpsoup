namespace dnpsoup {

  template<typename T>
  ShieldingInteraction<T>::ShieldingInteraction(double beta, size_t n)
    : InteractionInterface(), m_beta(beta), m_n(n), m_nbefore(0), m_nafter(0)
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
  ShieldingInteraction<T>::ShieldingInteraction(double beta, size_t n, size_t nbefore, size_t nafter)
    : InteractionInterface(), m_beta(beta), m_n(n), m_nbefore(nbefore), m_nafter(nafter)
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
  MatrixCxDbl ShieldingInteraction<T>::genMatrix(
      const Property &g,
      const Euler &e) const
  {
    const double gzz = g.get(ValueName::zz);
    const double gxx = g.get(ValueName::xx);
    const double gyy = g.get(ValueName::yy);
    const double bz = g.get(ValueName::bz);
    const double offset = g.get(ValueName::offset);

    const double sb = sin(e.beta());
    const double cb = cos(e.beta());
    const double sg = sin(e.gamma());
    const double cg = cos(e.gamma());
    double coeff3 = gzz * cb * cb + sb * sb * (gxx * cg * cg + gyy * sg * sg);

    if constexpr(std::is_same<T, LabFrame>::value){
      const double sa = sin(e.alpha());
      const double ca = cos(e.alpha());
      const double s2b = sin(2.0*e.beta());
      const double s2g = sin(2.0*e.gamma());

      double coeff1 = (gxx - gyy) * 0.5 * sa * sb * s2g + ca * s2b * (gzz - gxx * (cg * cg) - gyy * sg * sg);
      double coeff2 = (gyy - gxx) * 0.5 * ca * sb * s2g + sa * s2b * (gzz - gxx * (cg * cg) - gyy * sg * sg);
      MatrixCxDbl res = -m_beta * bz * (coeff1 * m_x + coeff2 * m_y + coeff3 * m_z);
      return res;
    } else {  // Rotating Frame
      MatrixCxDbl res = (-m_beta * bz * coeff3 - offset) * m_z;
      return res;
    }
  }

  template<typename T>
  size_t ShieldingInteraction<T>::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }

} // namespace dnpsoup

