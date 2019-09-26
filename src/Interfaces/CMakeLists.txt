SET(INTERFACEBASE_SRCS
   ESInterfaceBase.cpp ESHDF5/ESHDF5Interface.cpp  InterfaceBuilder.cpp PWSCF/ESPWSCFInterface.cpp
  )
INCLUDE_DIRECTORIES(${qmcpack_SOURCE_DIR}/src/Interfaces
${qmcpack_SOURCE_DIR}/src/Interfaces/ESHDF5 ${qmcpack_SOURCE_DIR}/src/
${qmcpack_SOURCE_DIR}/src/Interfaces/PWSCF)
ADD_LIBRARY(interface ${INTERFACEBASE_SRCS})
TARGET_LINK_LIBRARIES(interface ${qmcpack_SOURCE_DIR}/external_codes/quantum_espresso/q-e-qe-5.3/PW/src/libpwinterface.so mpi gfortran mpi_mpifh ) # MR: I'm not 100% sure this is actually required

