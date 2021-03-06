/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_STRUCTS_STD_H
#define MULTI_BSPLINE_STRUCTS_STD_H

#include <stdlib.h>
#include <stdint.h>

///////////////////////////
// Single precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride;
  Ugrid x_grid;
  BCtype_s xBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
  int num_splines;
} multi_UBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_3d_s;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride;
  Ugrid x_grid;
  BCtype_d xBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_1d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride, y_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
  int num_splines;
} multi_UBspline_2d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_3d_d;



//////////////////////////////
// Single precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride;
  Ugrid x_grid;
  BCtype_c xBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_1d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride, y_stride;
  Ugrid x_grid, y_grid;
  BCtype_c xBC, yBC;
  int num_splines;
  // temporary storage for laplacian components
  complex_float* restrict lapl2;
} multi_UBspline_2d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_c xBC, yBC, zBC;
  int num_splines;
  size_t coefs_size;
  // temporary storage for laplacian components
  complex_float* restrict lapl3;
} multi_UBspline_3d_c;


//////////////////////////////
// Double precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride;
  Ugrid x_grid;
  BCtype_z xBC;
  int num_splines;
  size_t coefs_size;
} multi_UBspline_1d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride, y_stride;
  Ugrid x_grid, y_grid;
  BCtype_z xBC, yBC;
  int num_splines;
  // temporary storage for laplacian components
  complex_double* restrict lapl2;
} multi_UBspline_2d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
  int num_splines;
  size_t coefs_size;
  // temporary storage for laplacian components
  complex_double* restrict lapl3;
} multi_UBspline_3d_z;


#endif
