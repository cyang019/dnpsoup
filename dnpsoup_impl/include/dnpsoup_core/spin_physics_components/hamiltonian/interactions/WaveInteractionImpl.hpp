namespace dnpsoup {

  template<typename T>
  WaveInteraction<T>::WaveInteraction(const std::vector<SpinType> &spins, const SpinType &irradiated)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> wave_operators;

    if constexpr(std::is_same<T, LabFrame>::value){
      m_x = spin<X>(m_n);
      m_y = spin<Y>(m_n);
      m_z = spin<Z>(m_n);
    } else {  // Rotating Frame
      m_z = spin<Z>(m_n);
    }
  }

  template<typename T>
  MatrixCxDbl WaveInteraction<T>::genMatrix(
      const Property &p,
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
      MatrixCxDbl res = - m_gamma * bz * m_z + coeff1 * m_x + coeff2 * m_y + coeff3 * m_z;
      return 2 * pi * res;
    } else {  // Rotating Frame
      MatrixCxDbl res = coeff3 * m_z;
      return 2 * pi * res;
    }
  }

  template<typename T>
  size_t WaveInteraction<T>::dimension() const
  { 
    return m_ntotal;
  }

} // namespace dnpsoup

