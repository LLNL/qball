////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
//  BLAS Header file
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BLAS_H
#define BLAS_H

#include <config.h>
#include <complex>
using namespace std;

// default value for most compilers
#define FTN_LINK extern "C"

#ifdef USE_JAGGEMM
#define zgemm jag_zgemm
#define dgemm jag_dgemm
#endif

#define dcopy  FC_FUNC(dcopy, DCOPY) 
#define zcopy  FC_FUNC(zcopy, ZCOPY) 
#define daxpy  FC_FUNC(daxpy, DAXPY) 
#define zaxpy  FC_FUNC(zaxpy, ZAXPY) 
#define ddot   FC_FUNC(ddot,  DDOT)  
#define zdotu  FC_FUNC(zdotu, ZDOTU)  
#define zdotc  FC_FUNC(zdotc, ZDOTC)  
#define drot   FC_FUNC(drot,  DROT)  
#define dasum  FC_FUNC(dasum, DASUM) 
#define dsbmv  FC_FUNC(dsbmv, DSBMV) 
#define dgemm  FC_FUNC(dgemm, DGEMM) 
#define zgemm  FC_FUNC(zgemm, ZGEMM) 
#define dgesv  FC_FUNC(dgesv, DGESV) 
#define dgemv  FC_FUNC(dgemv, DGEMV) 
#define dscal  FC_FUNC(dscal, DSCAL) 
#define dsysv  FC_FUNC(dsysv, DSYSV)
#define dsyev  FC_FUNC(dsyev, DSYEV) 
#define zdscal FC_FUNC(zdscal, ZDSCAL)
#define idamax FC_FUNC(idamax, IDAMAX)
#define dvea   FC_FUNC(dvea, DVEA)  
#define dyax   FC_FUNC(dyax, DYAX)  
#define dnaxpy FC_FUNC(dnaxpy, DNAXPY)
#define dger   FC_FUNC(dger, DGER)  
#define zgeru  FC_FUNC(zgeru, ZGERU)  
#define zgerc  FC_FUNC(zgerc, ZGERC)  
#define dgetrf FC_FUNC(dgetrf, DGETRF)
#define dgetri FC_FUNC(dgetri, DGETRI)

#ifdef USE_JAGGEMM
#define zgemm jag_zgemm
#define dgemm jag_dgemm
#endif


#ifdef __cplusplus
FTN_LINK {
#endif

void dcopy(int *n, double *x, int *incx, 
double *y, int *incy );
void zcopy(int *n, complex<double> *x, int *incx, 
complex<double> *y, int *incy );
void daxpy(int *n, double *alpha, double *x, int *incx,
double *y, int *incy );
void zaxpy(int *n, complex<double> *alpha, complex<double> *x, int *incx,
complex<double> *y, int *incy );
double ddot(const int *n, const double *a, const int *inca, 
const double *b, const int *incb);
complex<double> zdotu(const int *n, const complex<double> *a, const int *inca, 
  const complex<double> *b, const int *incb);
complex<double> zdotc(int* n, complex<double> *zx, int* incx, 
  complex<double> *zy, int* incy);
void drot(int*, double*, int*, double*, int*, double*, double*);
void dgemm(char *ta, char *tb, int *m, int *n, int *k,
  double *alpha, double *a, int *lda, double *b, int *ldb,
  double *beta, double *c, int *ldc);
void zgemm(char *ta, char *tb, int *m, int *n, int *k,
  complex<double> *alpha, complex<double> *a, int *lda, complex<double> *b, int *ldb,
  complex<double> *beta, complex<double> *c, int *ldc);
void dgemv( char *ta, int *m, int *n,
                   double *alpha,  double *a, int *tda,
                   double *x,    int *incx,
                   double *beta,   double *y, int *incy );
 
void dger(int *,int *, double *, double *, int *,
          double *, int *, double *, int *);
void zgeru(int *,int *, complex<double> *, complex<double> *, int *,
          complex<double> *, int *, complex<double> *, int *);
void zgerc(int *,int *, complex<double> *, complex<double> *, int *,
          complex<double> *, int *, complex<double> *, int *);
 
void dscal(int *len, double *alpha, double *x, int *incx);
double dasum(int *len, double *x, int *incx);
int idamax(int *len, double *x, int *incx);
void dsyev(char *c1,char *c2,int *n, 
double *a,int *lda, double *wr,
double *wrk,int *lwrk, int *ierr);
void dsysv(char*,int*,int*,double*,int*,int*,double*,int*,double*,int*,int*);
void zdscal_(int *n,double *alpha,complex<double> *x,int *incx);
void dgbmv(char *trans, int *m, int *n,
int *kl, int *ku, double *alpha, double *a, 
int *lda, double *x, int *incx, double *beta,
double *y, int *incy);
void dsbmv(char *uplo, int *n, int *k,
double *alpha, double *a, int *lda, double *x, int *incx,
double *beta, double *y, int *incy);
void sspev(char *vec,char *uplo,int *size,double *ap,
double *wr,double *z,int *n,double *wrk,int *ierr);
void dgesv(int *n, int *nrhs, double *a, int *lda, int *ipiv, 
double *b, int *ldb, int *info);
void dgetrf(const int*, const int*, double*, const int*, int*, int*);
void dgetri(const int*, double*, const int*, int*, double*, int*, int*);
 
void dvea(int*,double*,int*,double*,int*,double*,int*);
void dyax(int*,double*,double*,int*,double*,int*);
void dnaxpy(int*,int*,double*,int*,double*,int*,int*,double*,int*,int*);
 
#ifdef __cplusplus
}
#endif
 
#endif
