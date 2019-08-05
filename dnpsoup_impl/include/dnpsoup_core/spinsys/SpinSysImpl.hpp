namespace dnpsoup {
    template<typename T>
    PacketCollection SpinSys::summarize() const
    {
      PacketCollection result;
      for(const auto &obs : m_observables){
        auto ptr_interaction = genInteractionFromObservable<T>(obs.second);
        auto packet = HamiltonianPacket(
            std::move(ptr_interaction), obs.second.getProperty(), obs.second.getEuler());
        result.add(obs.first, std::move(packet));
      }
      return result;
    }

  // private methods
  template<typename T>
  std::unique_ptr<InteractionInterface> 
  SpinSys::genInteractionFromObservable(const Observable &ob) const
  {
    const auto ob_type = ob.getType();
    std::unique_ptr<InteractionInterface> res = nullptr;
    switch(ob_type){
      case InteractionType::Csa:
        {
          SpinId sid = ob.getSpinIds()[0];
          SpinType t = m_spins.at(sid).getSpinType();
          std::size_t n = getMatrixDimension(t);
          std::size_t nbefore = calcDimBeforeId(m_spins, sid);
          std::size_t nafter = calcDimAfterId(m_spins, sid);
          const double gyro = getGyromagneticRatio(t);
          if constexpr (std::is_same<T, DnpExperiment>::value){
            res = std::make_unique<ChemicalShiftInteraction<LabFrame>>(
                gyro, n, nbefore, nafter);
          } else {
            res = std::make_unique<ChemicalShiftInteraction<RotatingFrame>>(
                gyro, n, nbefore, nafter);
          }
        }
        break;
      case InteractionType::Shielding:
        {
          SpinId sid = ob.getSpinIds()[0];
          SpinType t = m_spins.at(sid).getSpinType();
          std::size_t n = getMatrixDimension(t);
          std::size_t nbefore = calcDimBeforeId(m_spins, sid);
          std::size_t nafter = calcDimAfterId(m_spins, sid);
          const double gyro = getGyromagneticRatio(t);
          res = std::make_unique<ShieldingInteraction<RotatingFrame>>(
              gyro, n, nbefore, nafter);
        }
        break;
      case InteractionType::Dipole:
        {
          SpinId sid1 = ob.getSpinIds()[0];
          SpinId sid2 = ob.getSpinIds()[1];
          SpinType t1 = m_spins.at(sid1).getSpinType();
          SpinType t2 = m_spins.at(sid2).getSpinType();
          std::size_t n1 = getMatrixDimension(t1);
          std::size_t n2 = getMatrixDimension(t2);
          std::size_t nbefore = calcDimBeforeId(m_spins, sid1);
          std::size_t nbetween = calcDimBetweenIds(m_spins, sid1, sid2);
          std::size_t nafter = calcDimAfterId(m_spins, sid2);
          const double gyro1 = getGyromagneticRatio(t1);
          const double gyro2 = getGyromagneticRatio(t2);
          if constexpr(std::is_same<T, DnpExperiment>::value){
            if(t1 == SpinType::e && t2 == SpinType::e){
              res = std::make_unique<
                DipolarInteraction<RotatingFrame, RotatingFrame>>(
                    gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);

            } else if(t1 == SpinType::e){
              res = std::make_unique<
                DipolarInteraction<RotatingFrame, LabFrame>>(
                    gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);
            } else if(t2 == SpinType::e){
              res = std::make_unique<
                DipolarInteraction<LabFrame, RotatingFrame>>(
                    gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);
            } else {
              res = std::make_unique<
                DipolarInteraction<LabFrame, LabFrame>>(
                    gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);
            }
          } 
          else {  // nmr experiment
              res = std::make_unique<
                DipolarInteraction<RotatingFrame, RotatingFrame>>(
                    gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);
          }
        }
        break;
      case InteractionType::Scalar:
        {
          SpinId sid1 = ob.getSpinIds()[0];
          SpinId sid2 = ob.getSpinIds()[1];
          SpinType t1 = m_spins.at(sid1).getSpinType();
          SpinType t2 = m_spins.at(sid2).getSpinType();
          std::size_t n1 = getMatrixDimension(t1);
          std::size_t n2 = getMatrixDimension(t2);
          std::size_t nbefore = calcDimBeforeId(m_spins, sid1);
          std::size_t nbetween = calcDimBetweenIds(m_spins, sid1, sid2);
          std::size_t nafter = calcDimAfterId(m_spins, sid2);
          const double gyro1 = getGyromagneticRatio(t1);
          const double gyro2 = getGyromagneticRatio(t2);
          res = std::make_unique<
            ScalarInteraction>(
                gyro1, gyro2, n1, n2, nbefore, nbetween, nafter);
        }
        break;
      case InteractionType::Wave:
        {
          std::vector<SpinType> types = this->getSpinTypes();
          if(ob.getSpinIds().size() == 0){
            /// either frame is fine, since nothing irradiated.
            res = std::make_unique<
              WaveInteraction<RotatingFrame>>(types, SpinType::Null);
          } else {
            /// irradiated on the same type.
            const SpinId sid0 = ob.getSpinIds()[0];
            const SpinType t = m_spins.at(sid0).getSpinType();
            if constexpr(std::is_same<T, DnpExperiment>::value){
              if(t == SpinType::e){
                res = std::make_unique<
                  WaveInteraction<RotatingFrame>>(types, t);
              } else {
                res = std::make_unique<
                  WaveInteraction<LabFrame>>(types, t);
              }
            } else {
                res = std::make_unique<
                  WaveInteraction<RotatingFrame>>(types, t);
            }
          }
        }
        break;
      default:
        throw NotImplementedError("Unknown InteractionType");
        break;
    }

    return res;
  }
} // namespace dnpsoup
