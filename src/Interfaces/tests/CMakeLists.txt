#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////
    
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR interfaces)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/interfaces)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)

SET(FILES_TO_COPY pwscf.in O.BFD.upf)

FOREACH(fname ${FILES_TO_COPY})
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${fname}" ${UTEST_DIR})
ENDFOREACH()

IF ((NOT QMC_COMPLEX) AND (NOT QMC_MIXED_PRECISION) AND (NOT QMC_CUDA))
  SET(MO_SRCS test_MO.cpp test_multiquintic_spline.cpp)
  IF (ENABLE_SOA)
    SET(MO_SRCS ${MO_SRCS} test_soa_cusp_corr.cpp)
  ELSE()
    SET(MO_SRCS ${MO_SRCS} test_cusp_corr.cpp)
  ENDIF()
ENDIF()

ADD_EXECUTABLE(${UTEST_EXE}
               test_interface_HDF5.cpp test_interface_PWSCF.cpp
               ${MO_SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcwfs qmcbase qmcutil interface ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
