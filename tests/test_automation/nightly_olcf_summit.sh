#!/bin/bash

echo ""
echo ""
echo "=================================================="
echo "job $LS_JOBID running at $(date)"
echo "=================================================="
echo ""
echo ""

#
# Setup for summit.olcf.ornl.gov
#
# Run the "short" nightlies
#

module load cmake
module load git
module load gcc/6.4.0
module load fftw
module load essl
module load netlib-lapack
module load hdf5
module load boost

export TEST_SITE_NAME=summit.olcf.ornl.gov
export N_PROCS_BUILD=32
export N_CONCURRENT_TESTS=1
export CC=mpicc
export CXX=mpicxx
export FFTW_HOME=${OLCF_FFTW_ROOT}

base_dir=/gpfs/alpine/proj-shared/mat151/periodic_qmc #Must be an absolute path
qmc_data_dir=${base_dir}/h5data

compiler=gcc6.4
branch=develop

source_dir=${base_dir}/qmcpack-${branch}

echo
echo --- Hostname --- $(hostname -f)
echo --- Checkout for branch ${branch} $(date)
echo

# get or refresh the source checkout
if [ ! -e ${source_dir} ]; then
  mkdir -p ${base_dir} && cd ${base_dir}
  echo --- Cloning QMCPACK git $(date)
  git clone --depth 1 https://github.com/QMCPACK/qmcpack.git ${source_dir}
else
  echo --- Updating local QMCPACK git $(date)
  cd ${source_dir}
  git pull
fi

if [ -e ${source_dir}/CMakeLists.txt ]; then

  cd ${source_dir}

  git checkout ${branch}

  for variant in Real-SoA #Real-Mixed-SoA Complex-SoA Complex-Mixed-SoA
  do

    echo --- Building for ${variant} $(date)

    # set up a build-dir for this variant
    build_dir=${source_dir}/build_${compiler}_${variant}
    rm -rf ${build_dir}
    mkdir ${build_dir} && cd ${build_dir}

    # create log file folder if not exist
    log_dir=${base_dir}/log/$(date +%y_%m_%d)
    mkdir -p ${log_dir}

    # set up ctest and run it
    CTEST_FLAGS="-DQMC_DATA=${qmc_data_dir} \
      -DENABLE_TIMERS=1 \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_C_COMPILER=mpicc \
      -DMPIEXEC_NUMPROC_FLAG=${source_dir}/tests/scripts/jsrunhelper.sh"

    [[ ${variant} == *"Complex"* ]] && CTEST_FLAGS="${CTEST_FLAGS} -D QMC_COMPLEX=1"
    [[ ${variant} == *"-SoA"* ]] && CTEST_FLAGS="${CTEST_FLAGS} -D ENABLE_SOA=1"
    [[ ${variant} == *"-Mixed"* ]] && CTEST_FLAGS="${CTEST_FLAGS} -D QMC_MIXED_PRECISION=1"

    export QMCPACK_TEST_SUBMIT_NAME=${compiler}-${variant}-Release

    ctest ${CTEST_FLAGS} \
      -S $(pwd)/../CMake/ctest_script.cmake,release \
      --stop-time $(date --date=now+117mins +%H:%M:%S) \
      -VV -E 'long' --timeout 800 &> \
      ${log_dir}/${QMCPACK_TEST_SUBMIT_NAME}.log

    echo --- Finished ${variant} $(date)
  done

else
  echo  "ERROR: No CMakeLists. Bad git clone."
  exit 1
fi

#
# Submit the next job
#

hour=1
minute=05

echo ""
echo ""
echo "=================================================="
echo "submitting job for ${hour}:${minute}"
echo "=================================================="
echo ""
echo ""

bsub -b ${hour}:${minute} \
  -P mat151 \
  -W 2:00 \
  -nnodes 1 \
  -alloc_flags smt1 \
  -N \
  -J periodic_qmc \
  -o periodic_qmc.%J \
  -e periodic_qmc.%J \
  < ${source_dir}/tests/test_automation/nightly_olcf_summit.sh
