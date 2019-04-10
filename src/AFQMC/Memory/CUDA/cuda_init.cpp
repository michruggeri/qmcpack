//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////


#include "AFQMC/Memory/CUDA/cuda_init.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#include "Utilities/OutputManager.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

namespace qmc_cuda {

  // work around for problem with csrmm 
  extern boost::multi::array<std::complex<double>,1,qmc_cuda::cuda_gpu_allocator<std::complex<double>>>
                            *cusparse_buffer;

  extern cublasHandle_t afqmc_cublas_handle;
//extern  cublasXtHandle_t afqmc_cublasXt_handle;
  extern cusparseHandle_t afqmc_cusparse_handle;
  extern cusolverDnHandle_t afqmc_cusolverDn_handle;
  extern curandGenerator_t afqmc_curand_generator;
  extern bool afqmc_cuda_handles_init;
  extern cusparseMatDescr_t afqmc_cusparse_matrix_descr;

  extern std::vector<cudaStream_t> afqmc_cuda_streams;

  // need a cleanup routine
  void CUDA_INIT(boost::mpi3::shared_communicator& node, unsigned long long int iseed)
  {

    if(afqmc_cuda_handles_init) return;
    afqmc_cuda_handles_init=true;

    int num_devices=0;
    cudaGetDeviceCount(&num_devices);
    qmcplusplus::app_log()<<" Running in node with " <<num_devices <<" GPUs. \n";
    if(num_devices < node.size()) {
      qmcplusplus::app_error()<<"Error: # GPU < # tasks in node. " <<std::endl;
      qmcplusplus::app_error()<<"# GPU: " <<num_devices <<std::endl;
      qmcplusplus::app_error()<<"# tasks: " <<node.size() <<std::endl;
      APP_ABORT("");
    } else if(num_devices > node.size()) {
      qmcplusplus::app_log()<<"WARNING: Unused devices !!!!!!!!!!!!!! \n"
                                <<"         # tasks: " <<node.size() <<"\n"
                                <<"         num_devices: " <<num_devices <<std::endl;
    }

    cuda_check(cudaSetDevice(node.rank()),"cudaSetDevice()");

    cublas_check(cublasCreate (& afqmc_cublas_handle ), "cublasCreate");
//    cublas_check(cublasXtCreate (& afqmc_cublasXt_handle ), "cublasXtCreate");
    int devID[8] {0,1,2,3,4,5,6,7};
//    cublas_check(cublasXtDeviceSelect(afqmc_cublasXt_handle, 1, devID), "cublasXtDeviceSelect");
//    cublas_check(cublasXtSetPinningMemMode(afqmc_cublasXt_handle, CUBLASXT_PINNING_ENABLED), 
//                                            "cublasXtSetPinningMemMode");
    cusolver_check(cusolverDnCreate (& afqmc_cusolverDn_handle ), "cusolverDnCreate");
    //curand_check(curandCreateGenerator(&afqmc_curand_generator, CURAND_RNG_PSEUDO_DEFAULT),
    curand_check(curandCreateGenerator(&afqmc_curand_generator, CURAND_RNG_PSEUDO_MT19937),
                                            "curandCreateGenerator");
    curand_check(curandSetPseudoRandomGeneratorSeed(afqmc_curand_generator,iseed),
                                            "curandSetPseudoRandomGeneratorSeed");

    cusparse_check(cusparseCreate (& afqmc_cusparse_handle ), "cusparseCreate");
    cusparse_check(cusparseCreateMatDescr(&afqmc_cusparse_matrix_descr), 
            "cusparseCreateMatDescr: Matrix descriptor initialization failed"); 
    cusparseSetMatType(afqmc_cusparse_matrix_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(afqmc_cusparse_matrix_descr,CUSPARSE_INDEX_BASE_ZERO); 

    cusparse_buffer = new boost::multi::array<std::complex<double>,1,
                                 qmc_cuda::cuda_gpu_allocator<std::complex<double>>>(
                                 (typename boost::multi::layout_t<1u>::extensions_type{1},
                                 qmc_cuda::cuda_gpu_allocator<std::complex<double>>{}));

  }

}

