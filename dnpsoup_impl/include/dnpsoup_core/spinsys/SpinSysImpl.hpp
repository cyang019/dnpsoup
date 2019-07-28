namespace dnpsoup {
    template<typename T>
    SpinPacketCollection SpinSys::Summarize() const
    {
      if constexpr(std::is_same<T, DnpExperiment>::value){
        // e in rotating frame, everything else in lab frame
      } else if constexpr(std::is_same<T, NmrExperiment>::value){
        // everything in the rotating frame
      } else{
        throw NotImplementedError("Summarize() only defined for DnpExperiment and NmrExperiment types.");
      }
    }
} // namespace dnpsoup
