#!/bin/bash -x

BUILD_DIR=$(pwd)
echo $BUILD_DIR

cat > $BUILD_TAG.pbs << EOF
#PBS -A MAT151ci
#PBS -N $BUILD_TAG
#PBS -j oe
#PBS -l walltime=1:00:00,nodes=1
#PBS -d $BUILD_DIR
#PBS -l partition=gpu

source /sw/rhea/lmod/7.8.2/rhel7.5_4.8.5/lmod/7.8.2/init/bash

module unload intel
module load gcc/6.2.0
module load openblas/0.3.5
module load cuda/10.0.130
module load fftw/3.3.8
export FFTW_HOME=\$OLCF_FFTW_ROOT
module load hdf5/1.10.3
module load git/2.18.0
module load cmake/3.13.4
module load boost/1.67.0

#// this is a workaround for picking up the exact libstdc++
#export LD_PRELOAD=/ccs/compilers/gcc/rhel7-x86_64/6.2.0/lib64/libstdc++.so.6

env
module list

echo ""
echo ""
echo "cache source in tmpfs"
echo "at $(date)"
echo ""
echo ""

rm -rf /tmp/${BUILD_TAG}-src
cp -R ${BUILD_DIR} /tmp/${BUILD_TAG}-src


echo ""
echo ""
echo "starting new test for real full precision"
echo "at $(date)"
echo ""
echo ""

mkdir -p /tmp/${BUILD_TAG}-build
cd /tmp/${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DENABLE_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

# hacky way to check on cmake. works for now
if ! ( grep -- '-- The C compiler identification is GNU 6.2.0' cmake.out && \
       grep -- '-- The CXX compiler identification is GNU 6.2.0' cmake.out ) ;
then
  echo "compiler version mismatch. exiting."
  exit 1
fi

# because Andreas tells me (and I observe) that GPU builds are unstable with Cmake
time make -j 24
time make -j 24

time ctest -L unit --output-on-failure



echo ""
echo ""
echo "starting new test for real mixed precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DENABLE_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure



echo ""
echo ""
echo "starting new test for complex full precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DENABLE_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure


echo ""
echo ""
echo "starting new test for complex mixed precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DENABLE_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure

# final cleanup of tmpfs
cd ../
rm -rf ./${BUILD_TAG}-build
rm -rf ./${BUILD_TAG}-src

EOF

/home/mat151ci_auser/blocking_qsub $BUILD_DIR $BUILD_TAG.pbs

cp $BUILD_DIR/$BUILD_TAG.o* ../

## this end of job logic could probably be more elegant
## hacks to get us going

cp $BUILD_DIR/$BUILD_TAG.o* ../

# check for correct test output from all builds
if [ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -ne 4 ]
then
   echo; echo
   echo One or more build variants failed. Check the build log for details.
   echo; echo
fi

# set the return code for the script
[ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -eq 4 ]
