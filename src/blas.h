////////////////////////////////////////////////////////////////////////////////
//
//  BLAS Header file
//
////////////////////////////////////////////////////////////////////////////////
// $Id: blas.h,v 1.5 2010/08/26 17:44:16 draeger1 Exp $

#ifndef BLAS_H
#define BLAS_H

#include <complex>
using namespace std;

// default value for most compilers
#define FTN_LINK extern "C"

#ifdef ADD_
#define dcopy  dcopy_  
#define zcopy  zcopy_  
#define daxpy  daxpy_  
#define zaxpy  zaxpy_  
#define ddot   ddot_   
#define zdotu  zdotu_   
#define zdotc  zdotc_   
#define drot   drot_   
#define dasum  dasum_  
#define dsbmv  dsbmv_  
#define dgemm  dgemm_  
#define zgemm  zgemm_  
#define dgesv  dgesv_  
#define dgemv  dgemv_  
#define dscal  dscal_  
#define dsysv  dsysv_
#define dsyev  dsyev_  
#define zdscal zdscal_ 
#define idamax idamax_ 
#define dvea   dvea_   
#define dyax   dyax_   
#define dnaxpy dnaxpy_ 
#define dger   dger_   
#define zgeru   zgeru_   
#define zgerc   zgerc_   
#define dgetrf     dgetrf_
#define dgetri     dgetri_
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
