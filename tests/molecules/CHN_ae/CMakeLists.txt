 IF (NOT QMC_CUDA)
  IF (NOT QMC_COMPLEX)
#
# CHN molecular Cusp Correction test for all electron B3LYP using pyscf
# This test checks VMC energy. The use of Jastrow increases sensitivity to bad Cusp parameters    
# Energy from pyscf: E=-93.4096773055448
#
#   "kinetic" " 93.179 0.045") # kinetic energy
#   "totenergy" " -93.33934 0.00021 ") # total energy
#   "eeenergy" " 55.2876 0.0013") # e-e energy
#   "samples" "409600000 0") # samples


  LIST(APPEND CHN_SHORT "kinetic" "93.179 0.85") # kinetic energy
  LIST(APPEND CHN_SHORT "totenergy" "-93.33934 0.004") # total energy
  LIST(APPEND CHN_SHORT "eeenergy" "55.2876 0.023") # e-e energy
  LIST(APPEND CHN_SHORT "samples" "160000 0.0") # samples
  QMC_RUN_AND_CHECK(short-chn_ae_cuspCorrection-vmc
                    "${qmcpack_SOURCE_DIR}/tests/molecules/CHN_ae"
                    vmc_short
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 CHN_SHORT # VMC
                    )
 
  ELSE()
    MESSAGE_VERBOSE("Skipping CuspCorrection tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  ENDIF()
 ELSE()
    MESSAGE_VERBOSE("Skipping CuspCorrection tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
 ENDIF()
