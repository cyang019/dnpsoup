namespace dnpsoup {

  template<typename T>
  WaveInteraction<T>::WaveInteraction(const std::vector<SpinType> &spin_types, const SpinType &irradiated_type)
    : InteractionInterface(), m_ntotal(0)
  { 
    static_assert(is_frame_type<T>::value, "T needs to be either LabFrame or RotatingFrame");
    std::vector<MatrixCxDbl> operators_x;
    std::vector<MatrixCxDbl> operators_y;

    for(const auto &t : spin_types){
      size_t n = getMatrixDimension(t);
      if(n > 0){
        m_ntotal = (m_ntotal == 0) ? n : m_ntotal * n;
      }
      if (t == irradiated_type){
        auto temp_x = spin<X>(n);
        operators_x.emplace_back(std::move(temp_x));
        if constexpr(std::is_same<T, RotatingFrame>::value){
          auto temp_y = spin<Y>(n);
          operators_y.emplace_back(std::move(temp_y));
        }
      } else {
        operators_x.emplace_back(identity<cxdbl>(n));
        if constexpr(std::is_same<T, RotatingFrame>::value){
          operators_y.emplace_back(identity<cxdbl>(n));
        }
      }
    }

    if constexpr(std::is_same<T, RotatingFrame>::value){
      m_x = kron(operators_x);
      m_y = kron(operators_y);
    } else {  // Lab Frame
      m_x = kron(operators_x);
      m_y = MatrixCxDbl();
    }
  }

  template<typename T>
  MatrixCxDbl WaveInteraction<T>::genMatrix(
      const Property &p,
      [[maybe_unused]] const Euler &e) const
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
  size_t WaveInteraction<T>::dimension() const
  { 
    return m_ntotal;
  }

} // namespace dnpsoup

