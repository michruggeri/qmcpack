SET(INTERFACEBASE_SRCS
   ESInterfaceBase.cpp ESHDF5/ESHDF5Interface.cpp  InterfaceBuilder.cpp
  )
INCLUDE_DIRECTORIES(${qmcpack_SOURCE_DIR}/src/Interfaces
${qmcpack_SOURCE_DIR}/src/Interfaces/ESHDF5 ${qmcpack_SOURCE_DIR}/src/
${qmcpack_SOURCE_DIR}/src/Interfaces/PWSCF)
ADD_LIBRARY(interface ${INTERFACEBASE_SRCS})

