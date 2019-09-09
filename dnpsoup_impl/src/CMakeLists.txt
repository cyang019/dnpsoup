target_sources(dnpsoup_core
  PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}//common.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/spin.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/Property.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/rho_eq.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/interactions/ScalarInteraction.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/interactions/OffsetInteraction.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/rotation/Quaternion.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/rotation/Coordinate.cpp"

  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/relaxation.cpp"
  
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinId.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/Observable.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/HamiltonianPacket.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/RelaxationPacket.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinSys.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinEntity.cpp"

  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/EMRadiation.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/Pulse.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/Delay.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/PulseSequence.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/Section.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/SubSequenceInterface.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/pulseseq/seq_common.cpp"

  "${CMAKE_CURRENT_LIST_DIR}/experiment/DnpRunner.cpp"
  )

