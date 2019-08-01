target_sources(dnpsoup_core
  PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}//common.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/spin.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/Property.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/hamiltonian/interactions/ScalarInteraction.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/rotation/Quaternion.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spin_physics_components/rotation/Coordinate.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinId.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/Observable.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SimulationPacket.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinSys.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/spinsys/SpinEntity.cpp"
  )

