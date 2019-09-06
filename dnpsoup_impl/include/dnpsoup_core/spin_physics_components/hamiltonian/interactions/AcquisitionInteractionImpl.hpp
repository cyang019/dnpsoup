namespace dnpsoup {

  template<typename T>
  AcquisitionInteraction<T>::AcquisitionInteraction(const std::vector<SpinType> &spin_types, const SpinType &irradiated_type)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_p;

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
        auto temp_p = spin<P>(n);
        if(n > 0){
          nafter = nbefore > 0 ? m_ntotal / (n * nbefore) : m_ntotal / n;
        } else {
          nafter = nbefore > 0 ? m_ntotal / nbefore : m_ntotal;
        }
        MatrixCxDbl p_op = kron(std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_p, identity<cxdbl>(nafter)});
        operators_p.emplace_back(std::move(p_op));
      }
      if(n > 0){
        nbefore = (nbefore == 0) ? n : nbefore * n;
      }
    }

    m_p = zeros<cxdbl>(m_ntotal, m_ntotal);
    for(const auto &m : operators_p){
      m_p += m;
    }
  }

  template<typename T>
  AcquisitionInteraction::AcquisitionInteraction(
      const std::map<SpinId, SpinEntity> &spins,
      const std::vector<SpinId> &irradiated)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_p;

    for(auto [sid, s_info] : irradiated){
      std::size_t n = getMatrixDimension(s_info.getSpinType());
      if(n > 0){
        m_ntotal = (m_ntotal == 0) ? n : m_ntotal * n;
      }
    }

    std::size_t nbefore = 0;
    std::size_t nafter = 0;
    for(const auto &sid : irradiated){
      const std::size_t n = getMatrixDimension(spins[sid].getSpinType());
      auto temp_p = spin<P>(n);
      if(n > 0){
        nafter = nbefore > 0 ? m_ntotal / (n * nbefore) : m_ntotal / n;
      } else {
        nafter = nbefore > 0 ? m_ntotal / nbefore : m_ntotal;
      }
      MatrixCxDbl p_op = kron(
          std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_p, identity<cxdbl>(nafter)});
      operators_p.emplace_back(std::move(p_op));
      if(n > 0){
        nbefore = (nbefore == 0) ? n : nbefore * n;
      }
    }

    m_p = zeros<cxdbl>(m_ntotal, m_ntotal);
    for(const auto &m : operators_p){
      m_p += m;
    }
  }

  template<typename T>
  MatrixCxDbl AcquisitionInteraction<T>::genMatrix(
      [[maybe_unused]] const Property &p,
      [[maybe_unused]] const Euler<ActiveRotation> &e) const
  {
    return m_p;
  }

  template<typename T>
  MatrixCxDbl AcquisitionInteraction<T>::genMatrix(
      [[maybe_unused]] const Property &p,
      [[maybe_unused]] const Euler<PassiveRotation> &e) const
  {
    return m_p;
  }

  template<typename T>
  std::size_t AcquisitionInteraction<T>::dimension() const
  { 
    return m_ntotal;
  }

} // namespace dnpsoup


