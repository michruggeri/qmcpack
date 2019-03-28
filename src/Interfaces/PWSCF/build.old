export QE_HOME=/home/yyang173/soft/espresso-5.1

ifort -fPIC -openmp -mkl=parallel -O2 -assume byterecl -g -traceback -par-report0 -vec-report0 -nomodule -fpp -D__OPENMP -D__INTEL -D__FFTW3 -DEXX -DH5_USE_16_API  -I${QE_HOME}/include -I/usr/gapps/qmc/libs/INTEL/hdf5-1.8.8/include/  -I/opt/intel-14.0/mkl/include -I/usr/include/mpi -I${QE_HOME}/iotk/src -I${QE_HOME}/Modules -I${QE_HOME}/PW/src/ -c pwinterfacemod.f90 

ifort -fPIC -openmp -mkl=parallel -O2 -assume byterecl -g -traceback -par-report0 -vec-report0 -nomodule -fpp -D__OPENMP -D__INTEL -D__FFTW3 -DEXX -DH5_USE_16_API  -I${QE_HOME}/include -I/usr/gapps/qmc/libs/INTEL/hdf5-1.8.8/include/  -I/opt/intel-14.0/mkl/include -I/usr/include/mpi -I${QE_HOME}/iotk/src -I${QE_HOME}/Modules -I${QE_HOME}/PW/src/ -c pwlib_getinfo.f90 

ifort -fPIC -openmp -mkl=parallel -O2 -assume byterecl -g -traceback -par-report0 -vec-report0 -nomodule -fpp -D__OPENMP -D__INTEL -D__FFTW3 -DEXX -DH5_USE_16_API  -I${QE_HOME}/include -I/usr/gapps/qmc/libs/INTEL/hdf5-1.8.8/include/  -I/opt/intel-14.0/mkl/include -I/usr/include/mpi -I${QE_HOME}/iotk/src -I${QE_HOME}/Modules -I${QE_HOME}/PW/src/ -c pwlib_init.f90 

ifort -fPIC -openmp -mkl=parallel -O2 -assume byterecl -g -traceback -par-report0 -vec-report0 -nomodule -fpp -D__OPENMP -D__INTEL -D__FFTW3 -DEXX -DH5_USE_16_API  -I${QE_HOME}/include -I/usr/gapps/qmc/libs/INTEL/hdf5-1.8.8/include/  -I/opt/intel-14.0/mkl/include -I/usr/include/mpi -I${QE_HOME}/iotk/src -I${QE_HOME}/Modules -I${QE_HOME}/PW/src/ -c pwlib_scf.f90 

#mpic++ -c mpitest.cpp

mpiicpc -openmp -shared -o libpwinterface.so pwlib_getinfo.o pwlib_scf.o pwlib_init.o pwinterfacemod.o ${QE_HOME}/PW/src/libpw.a ${QE_HOME}/Modules/libqemod.a ${QE_HOME}/flib/ptools.a ${QE_HOME}/flib/flib.a ${QE_HOME}/clib/clib.a ${QE_HOME}/iotk/src/libiotk.a  ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl /opt/intel15/apps/qmcpack_rev6964/build/lib/libeinspline.a -L/usr/local/lib -lhdf5 -lhdf5_hl -lifcore -lifport -fPIC
 
#mpic++ -openmp -o wee mpitest.o pwlib_getinfo.o pwlib_scf.o pwlib_init.o pwinterfacemod.o ${QE_HOME}/PW/src/libpw.a ${QE_HOME}/Modules/libqemod.a ${QE_HOME}/flib/ptools.a ${QE_HOME}/flib/flib.a ${QE_HOME}/clib/clib.a ${QE_HOME}/iotk/src/libiotk.a -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64  -lfftw3  -L/usr/local/tools/mkl-10.3.8/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm -Wl,-rpath,/usr/local/tools/mkl-10.3.8/mkl/lib/intel64    /g/g14/clay8/lscratchd/qmcpack/build/lib/libeinspline.a -L/usr/local/lib -lhdf5 -lhdf5_hl -lifcore -lifport
