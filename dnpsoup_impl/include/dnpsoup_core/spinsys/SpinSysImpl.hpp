namespace dnpsoup {
    template<typename T>
    SpinPacketCollection SpinSys::Summarize() const
    {
      std::size_t ntotal = calcTotalDimension();
      SpinPacketCollection spin_packets;
      std::vector<std::size_t> dimensions = calcDimensions();

      if constexpr(std::is_same<T, DnpExperiment>::value){
        // e in rotating frame, everything else in lab frame
        for(const auto &ob : m_observables){
        }
      } else if constexpr(std::is_same<T, NmrExperiment>::value){
        // everything in the rotating frame
      } else{
        throw NotImplementedError(
            "Summarize() only defined for DnpExperiment and NmrExperiment types.");
      }
    }

  // private methods
  template<typename T>
  std::unique_ptr<InteractionInterface> 
  SpinSys::genInteractionFromObservable(const Observable &ob) const
  {
    std::vector<std::size_t> dims = self.calcDimensions();
    const auto ob_type = ob.getType();
    std::unique_ptr<InteractionInterface> res = nullptr;
    switch(ob_type){
    case InteractionType::Csa:
      break;
    case InteractionType::Dipole:
      break;
    case InteractionType::Scalar:
      break;
    case InteractionType::Shielding:
      break;
    default:
      break;
    }
  }
} // namespace dnpsoup
