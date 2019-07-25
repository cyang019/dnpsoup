namespace dnpsoup {

  template<typename T>
  WaveInteraction<T>::WaveInteraction(const std::vector<SpinType> &spin_types, const SpinType &irradiated_type)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_x;
    std::vector<MatrixCxDbl> operators_y;

    for(const auto &t : spin_types){
      std::size_t n = getMatrixDimension(t);
      if(n > 0){
        m_ntotal = (m_ntotal == 0) ? n : m_ntotal * n;
      }
    }

    std::size_t nbefore = 0;
    std::size_t nafter = 0;
    for(const auto &t : spin_types){
      std::size_t n = getMatrixDimension(t);
      if (t == irradiated_type){
        auto temp_x = spin<X>(n);
        nafter = nbefore > 0 ? m_ntotal / (n * nbefore) : m_ntotal / n;
        MatrixCxDbl x_op = kron(std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_x, identity<cxdbl>(nafter)});
        operators_x.emplace_back(std::move(x_op));
        if constexpr(std::is_same<T, RotatingFrame>::value){
          auto temp_y = spin<Y>(n);
          MatrixCxDbl y_op = kron(std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_y, identity<cxdbl>(nafter)});
          operators_y.emplace_back(std::move(y_op));
        }
      }
      if(n > 0){
        nbefore = (nbefore == 0) ? n : nbefore * n;
      }
    }

    m_x = zeros<cxdbl>(m_ntotal, m_ntotal);
    for(const auto &m : operators_x){
      m_x += m;
    }
    if constexpr(std::is_same<T, RotatingFrame>::value){
      m_y = zeros<cxdbl>(m_ntotal, m_ntotal);
      for(const auto &m : operators_y){
        m_y += m;
      }
    } else {  // Lab Frame, just a placeholder
      m_y = MatrixCxDbl();
    }
  }

  template<typename T>
  MatrixCxDbl WaveInteraction<T>::genMatrixInternal(const Property &p) const
  {
    const double freq = p.get(ValueName::freq);
    const double phase = p.get(ValueName::phase);

    if constexpr(std::is_same<T, RotatingFrame>::value){
      MatrixCxDbl res = freq * (m_x * std::cos(phase) + m_y * std::sin(phase));
      return res;
    } else {  // Lab Frame
      const double phase0 = p.get(ValueName::phase0);
      MatrixCxDbl res = freq * m_x * std::cos(phase0 + phase);
      return res;
    }
  }
  
  template<typename T>
  MatrixCxDbl WaveInteraction<T>::genMatrix(
      const Property &p,
      [[maybe_unused]] const Euler<ActiveRotation> &e) const
  {
    return genMatrixInternal(p);
  }

  template<typename T>
  MatrixCxDbl WaveInteraction<T>::genMatrix(
      const Property &p,
      [[maybe_unused]] const Euler<PassiveRotation> &e) const
  {
    return genMatrixInternal(p);
  }

  template<typename T>
  std::size_t WaveInteraction<T>::dimension() const
  { 
    return m_ntotal;
  }

} // namespace dnpsoup

