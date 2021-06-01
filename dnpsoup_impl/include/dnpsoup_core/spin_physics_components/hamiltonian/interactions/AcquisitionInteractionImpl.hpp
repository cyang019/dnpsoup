namespace dnpsoup {

  template<typename T, typename E>
  AcquisitionInteraction<T, E>::AcquisitionInteraction(const std::vector<SpinType> &spin_types, const SpinType &irradiated_type)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_m;

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
        if(n > 0){
          nafter = nbefore > 0 ? m_ntotal / (n * nbefore) : m_ntotal / n;
        } else {
          nafter = nbefore > 0 ? m_ntotal / nbefore : m_ntotal;
        }
        if constexpr (std::is_same<E, NmrExperiment>::value){
          auto temp_m = spin<M>(n);
          MatrixCxDbl minus_op = kron(std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_m, identity<cxdbl>(nafter)});
          operators_m.emplace_back(std::move(minus_op));
        } else {
          auto temp_z = spin<Z>(n);
          MatrixCxDbl z_op = kron(std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_z, identity<cxdbl>(nafter)});
          operators_m.emplace_back(std::move(z_op));
        }
      }
      if(n > 0){
        nbefore = (nbefore == 0) ? n : nbefore * n;
      }
    }

    m_op = zeros<cxdbl>(m_ntotal, m_ntotal);
    for(const auto &m : operators_m){
      m_op += m;
    }
  }

  template<typename T, typename E>
  AcquisitionInteraction<T, E>::AcquisitionInteraction(
      const std::map<SpinId, SpinEntity> &spins,
      const std::vector<SpinId> &irradiated)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_m;

    for(const auto &sid : irradiated){
      auto s_info = spins.at(sid);
      std::size_t n = getMatrixDimension(s_info.getSpinType());
      if(n > 0){
        m_ntotal = (m_ntotal == 0) ? n : m_ntotal * n;
      }
    }

    std::size_t nbefore = 0;
    std::size_t nafter = 0;
    for(const auto &sid : irradiated){
      const std::size_t n = getMatrixDimension(spins.at(sid).getSpinType());
      if(n > 0){
        nafter = nbefore > 0 ? m_ntotal / (n * nbefore) : m_ntotal / n;
      } else {
        nafter = nbefore > 0 ? m_ntotal / nbefore : m_ntotal;
      }
      if constexpr(std::is_same<E, NmrExperiment>::value){
        auto temp_m = spin<M>(n);
        MatrixCxDbl minus_op = kron(
            std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_m, identity<cxdbl>(nafter)});
        operators_m.emplace_back(std::move(minus_op));
      }
      else {
        auto temp_z = spin<Z>(n);
        MatrixCxDbl z_op = kron(
            std::vector<MatrixCxDbl>{identity<cxdbl>(nbefore), temp_z, identity<cxdbl>(nafter)});
        operators_m.emplace_back(std::move(z_op));
      }
      if(n > 0){
        nbefore = (nbefore == 0) ? n : nbefore * n;
      }
    }

    m_op = zeros<cxdbl>(m_ntotal, m_ntotal);
    for(const auto &m : operators_m){
      m_op += m;
    }
  }

  template<typename T, typename E>
  MatrixCxDbl AcquisitionInteraction<T, E>::genMatrix(
      [[maybe_unused]] const Property &p,
      [[maybe_unused]] const Euler<ActiveRotation> &e,
      [[maybe_unused]] const Euler<ActiveRotation> &e2,
      [[maybe_unused]] const Euler<ActiveRotation> &e3
      ) const
  {
    return m_op;
  }

  template<typename T, typename E>
  MatrixCxDbl AcquisitionInteraction<T, E>::genMatrix(
      [[maybe_unused]] const Property &p,
      [[maybe_unused]] const Euler<PassiveRotation> &e,
      [[maybe_unused]] const Euler<PassiveRotation> &e2,
      [[maybe_unused]] const Euler<PassiveRotation> &e3
      ) const
  {
    return m_op;
  }

  template<typename T, typename E>
  std::size_t AcquisitionInteraction<T, E>::dimension() const
  { 
    return m_ntotal;
  }

} // namespace dnpsoup


