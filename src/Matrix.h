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
// Matrix.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef MATRIX_H
#define MATRIX_H

#define MATRIX_DEF_BLOCK_SIZE 64

class Context;

#include <valarray>
#include <complex>
#include <cstring>
#ifdef USE_CTF
#include "cyclopstf.hpp"
#endif
using namespace std;

class ComplexMatrix;

class DoubleMatrix
{
  private:
  
    Context ctxt_;
    int ictxt_;
    int lld_;  // leading dimension of local matrix
    int m_, n_;   // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_;   // size of local array
    int size_; // size of local array;
    int nprow_, npcol_; // number of process rows and cols in context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block
    bool reference_; // object was created using the copy constructor
    double* val;
#ifdef USE_CTF
    CTF* myctf_;
#endif
    
  public:
  
    double* valptr(int i=0) { return &val[i]; }
    const double* cvalptr(int i=0) const { return &val[i]; }
    double& operator[] (int i) { return val[i]; }
    const double& operator[] (int i) const { return val[i]; }
    const Context& context(void) const { return ctxt_; }
    void set_context(Context newctxt);
    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int nloc(void) const { return nloc_; } // size of local array
    int size(void) const { return size_; } // size of local array
    int localsize(void) const { return mloc_*nloc_; } // local size of val
    double memsize(void) const { return (double)m_ * (double)n_ * sizeof(double); }
    double localmemsize(void) const
    { return (double) mloc_ * (double) nloc_ * sizeof(double); }
    const int* desc(void) const { return &desc_[0]; }
    
    // local block size of block (l,m)
    int mbs(int l) const
    {
      // if l is the last block and this process holds an incomplete
      // block, then the size is m_%mb_. Otherwise, it is mb_.
      return ( (l==mblocks_-1) && ( m_incomplete_ ) ) ? m_%mb_ : mb_;
    }
    int nbs(int m) const
    {
      // if m is the last block and this process holds an incomplete
      // block, then the size is n_%nb_. Otherwise, it is nb_.
      return ( (m==nblocks_-1) && ( n_incomplete_ ) ) ? n_%nb_ : nb_;
    }
    
    // number of local blocks
    int mblocks(void) const { return mblocks_; }     
    int nblocks(void) const { return nblocks_; }
     
    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)    
    int l(int i) const { return i/(nprow_*mb_); }
    int x(int i) const { return i % mb_; }
    int pr(int i) const { return (i/mb_) % nprow_; }
    
    int m(int j) const { return j/(npcol_*nb_); }
    int y(int j) const { return j % nb_; }
    int pc(int j) const { return (j/nb_) % npcol_; }
    
    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(int l, int x) const { return (l * nprow_ + myrow_) * mb_ + x; }     
    int j(int m, int y) const { return (m * npcol_ + mycol_) * nb_ + y; }
    
    // store element a(ii,jj) (where ii,jj are global indices)
    // in array val:
    //        int iii = l(ii) * mb_ + x(ii);
    //        int jjj = m(jj) * nb_ + y(jj);
    //        val[iii+mloc_*jjj] = a_ij;
    
    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }
    
    void init_size(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE, 
      int nb = MATRIX_DEF_BLOCK_SIZE);
          
    void resize(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE, 
      int nb = MATRIX_DEF_BLOCK_SIZE)
    {
      const int old_size = size_;
      init_size(m,n,mb,nb);
      if ( size_ == old_size ) return;
      init_size(m,n,mb,nb);
      delete[] val;
      val = new double[size_];
      clear();
    }    
    
    void print(ostream& os) const;
    
    explicit DoubleMatrix(const Context& ctxt) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0) {}
        
    // Construct a DoubleMatrix of dimensions m,n
    explicit DoubleMatrix(const Context& ctxt, int m, int n,
        int mb=MATRIX_DEF_BLOCK_SIZE, 
        int nb=MATRIX_DEF_BLOCK_SIZE) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0)
    {
      resize(m,n,mb,nb);
    }
    
    // copy constructor: create a separate copy of rhs
    explicit DoubleMatrix(const DoubleMatrix& rhs) : ctxt_(rhs.context()),
      reference_(false)
    {
      init_size(rhs.m(),rhs.n(),rhs.mb(),rhs.nb());
      val = new double[size_];
      memcpy(val, rhs.val, size_*sizeof(double));
    }
    
    // reference constructor create a proxy for a ComplexMatrix rhs
    explicit DoubleMatrix(ComplexMatrix& rhs);
    explicit DoubleMatrix(const ComplexMatrix& rhs);
        
    ~DoubleMatrix(void) 
    { 
      if ( !reference_ ) delete[] val;
    }
                  
    DoubleMatrix& operator=(const DoubleMatrix& a);
    
    DoubleMatrix& operator+=(const DoubleMatrix& a);
    DoubleMatrix& operator-=(const DoubleMatrix& a);
    
    DoubleMatrix& operator*=(double a);
    
    void matgather(double *a, int lda) const;

    void initdiag(const double* const a); // initialize diagonal from a[]
    void init(const double* const a, int lda);
    void set(char uplo, double x);
    void identity(void);
    void clear(void);
    
    double dot(const DoubleMatrix &a) const;
    void axpy(double alpha, const DoubleMatrix &a);
    void scal(double alpha);
    
    double nrm2(void) const;
    double asum(void) const;
    double amax(void) const;
    double trace(void) const;
    void symmetrize(char uplo);
    
    // rank-1 update: *this += alpha * x(kx) * (y(ky))^T
    // where x(kx) is row kx of x, and y(ky) is row ky of y
    void ger(double alpha, const DoubleMatrix& x, int kx,
                           const DoubleMatrix& y, int ky);

    // symmetric rank-1 update
    void syr(char uplo, double alpha, const DoubleMatrix& x,
             int ix, char rowcol);
    
    // copy data in place between overlapping contexts
    void copyInPlace(DoubleMatrix& a);

    // get submatrix A(ia:ia+m,ja:ja+n) of A
    void getsub(const DoubleMatrix& a,int m,int n,int ia,int ja);

    // get submatrix A(ia:ia+m,ja:ja+n) of A and store in
    // this(idest:idest+m,jdest:jdest+n)
    void getsub(const DoubleMatrix& a,int m,int n,int isrc,int jsrc,
                int idest, int jdest);
    
    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(char transa, char transb,
         double alpha, const DoubleMatrix& a,
         const DoubleMatrix& b, double beta);
    
    // symmetric_matrix * matrix
    // *this = alpha * A * B + beta * this
    void symm(char side, char uplo,
         double alpha, const DoubleMatrix& a,
         const DoubleMatrix& b, double beta);
    
    // symmetric rank k update
    // this = beta * this + alpha * A * A^T  (trans=='n')
    // this = beta * this + alpha * A^T * A  (trans=='t')
    void syrk(char uplo, char trans,
         double alpha, const DoubleMatrix& a, double beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(double alpha, const DoubleMatrix& a, double beta);
    void transpose(const DoubleMatrix& a);

    void trmm(char side, char uplo, char trans, char diag, 
              double alpha, const DoubleMatrix& a);
    
    // solve triangular system
    void trsm(char side, char uplo, char trans, char diag, 
              double alpha, const DoubleMatrix& a);
    void trtrs(char uplo, char trans, char diag, DoubleMatrix& b) const;
    
    // Cholesky decomposition of a symmetric matrix
    void potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    void potri(char uplo);
    
    // LU decomposition
    void lu(valarray<int>& ipiv);
    
    // compute inverse of a square matrix
    void inverse(void);
    
    // Inverse of triangular matrix
    void trtri(char uplo,char diag);
    
    double pocon(char) const;
    void sygst(int, char, const DoubleMatrix&);
    
    // compute eigenvalues and eigenvectors of symmetric matrix *this
    void syev(char uplo, valarray<double>& w, DoubleMatrix& z);
    
    // compute eigenvalues and eigenvectors of symmetric matrix *this
    // using the divide and conquer method of Tisseur and Dongarra
    void syevd(char uplo, valarray<double>& w, DoubleMatrix& z);
    
    // compute eigenvalues (only) of symmetric matrix *this
    void syev(char uplo, valarray<double>& w);

    // compute eigenvalues (only) of symmetric matrix *this
    // using the divide and conquer method of Tisseur and Dongarra
    void syevd(char uplo, valarray<double>& w);
};
ostream& operator << ( ostream& os, const DoubleMatrix& a );

class ComplexMatrix
{
  private:
  
    Context ctxt_;
    int ictxt_;
    int lld_;  // leading dimension of local matrix
    int m_, n_;   // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_;   // size of local array
    int size_; // size of local array;
    int nprow_, npcol_; // number of process rows and cols in  context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block
    bool reference_; // object was created using the copy constructor
    complex<double>* val;
#ifdef USE_CTF
    tCTF< std::complex<double> > * myctf_;
#endif

  public:
  
    complex<double>* valptr(int i=0) { return &val[i]; }
    const complex<double>* cvalptr(int i=0) const { return &val[i]; }
    complex<double>& operator[] (int i) { return val[i]; }
    const complex<double>& operator[] (int i) const { return val[i]; }
    const Context& context(void) const { return ctxt_; }
    void set_context(Context newctxt);
    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int nloc(void) const { return nloc_; } // size of local array
    int size(void) const { return size_; } // size of local array
    int localsize(void) const { return mloc_*nloc_; } // local size of val
    double memsize(void) const 
    { return (double) m_ * (double) n_ * sizeof(complex<double>); }
    double localmemsize(void) const
    { return (double) mloc_ * (double) nloc_ * sizeof(complex<double>); }
    const int* desc(void) const { return &desc_[0]; }
    
    // local block size of block (l,m)
    int mbs(int l) const
    {
      // if l is the last block and this process holds an incomplete
      // block, then the size is m_%mb_. Otherwise, it is mb_.
      return ( (l==mblocks_-1) && ( m_incomplete_ ) ) ? m_%mb_ : mb_;
    }
    int nbs(int m) const
    {
      // if m is the last block and this process holds an incomplete
      // block, then the size is n_%nb_. Otherwise, it is nb_.
      return ( (m==nblocks_-1) && ( n_incomplete_ ) ) ? n_%nb_ : nb_;
    }
    
    // number of local blocks
    int mblocks(void) const { return mblocks_; }     
    int nblocks(void) const { return nblocks_; }
     
    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)    
    int l(int i) const { return i/(nprow_*mb_); }
    int x(int i) const { return i % mb_; }
    int pr(int i) const { return (i/mb_) % nprow_; }
    
    int m(int j) const { return j/(npcol_*nb_); }
    int y(int j) const { return j % nb_; }
    int pc(int j) const { return (j/nb_) % npcol_; }
    
    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(int l, int x) const { return (l * nprow_ + myrow_) * mb_ + x; }     
    int j(int m, int y) const { return (m * npcol_ + mycol_) * nb_ + y; }
    
    // store element a(ii,jj) (where ii,jj are global indices)
    // in array val:
    //        int iii = l(ii) * mb_ + x(ii);
    //        int jjj = m(jj) * nb_ + y(jj);
    //        val[iii+mloc_*jjj] = a_ij;
    
    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }
    
    void init_size(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE, 
      int nb = MATRIX_DEF_BLOCK_SIZE);
          
    void resize(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE, 
      int nb = MATRIX_DEF_BLOCK_SIZE)
    {
      const int old_size = size_;
      init_size(m,n,mb,nb);
      if ( size_ == old_size ) return;
      delete[] val;
      val = new complex<double>[size_];
      clear();
    }    
    
    void print(ostream& os) const;
    void print_norm(std::ostream& os) const;                  
    
    explicit ComplexMatrix(const Context& ctxt) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0) {}
        
    // Construct a ComplexMatrix of dimensions m,n
    explicit ComplexMatrix(const Context& ctxt, int m, int n,
        int mb=MATRIX_DEF_BLOCK_SIZE, 
        int nb=MATRIX_DEF_BLOCK_SIZE) : ctxt_(ctxt), 
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0)
    {
      resize(m,n,mb,nb);
    }
    
    // copy constructor: create a separate copy of rhs
    explicit ComplexMatrix(const ComplexMatrix& rhs) : ctxt_(rhs.context()),
      reference_(false)
    {
      init_size(rhs.m(),rhs.n(),rhs.mb(),rhs.nb());
      val = new complex<double>[size_];
      memcpy(val, rhs.val, size_*sizeof(complex<double>));
    }
    
    // reference constructor: create a proxy for a DoubleMatrix rhs
    explicit ComplexMatrix(DoubleMatrix& rhs);
    explicit ComplexMatrix(const DoubleMatrix& rhs);
    
    ~ComplexMatrix(void)
    { 
      if ( !reference_ ) delete[] val;
    }
    ComplexMatrix& operator=(const ComplexMatrix& a);
    
    ComplexMatrix& operator+=(const ComplexMatrix& a);
    ComplexMatrix& operator-=(const ComplexMatrix& a);
    
    ComplexMatrix& operator*=(double a);
    ComplexMatrix& operator*=(complex<double> a);
    
    void matgather(complex<double> *a, int lda) const;

    void initdiag(const complex<double>* const a); // initialize diagonal 
    void init(const complex<double>* const a, int lda);
    void set(char uplo, complex<double> x);
    void identity(void);
    void clear(void);
    
    complex<double> dot(const ComplexMatrix &a) const;  // tr A^H * A
    complex<double> dotu(const ComplexMatrix &a) const; // tr A^T * A
    void axpy(complex<double> alpha, const ComplexMatrix &a);
    void axpy(double alpha, const ComplexMatrix &a);
    void scal(complex<double> alpha);
    void scal(double alpha);
    
    double nrm2(void) const;
    double asum(void) const;
    double amax(void) const;
    complex<double> trace(void) const;
    void symmetrize(char uplo);
    
    // rank-1 update: *this += alpha * x(kx) * (y(ky))^T
    // where x(kx) is row kx of x, and y(ky) is row ky of y
    void ger(complex<double> alpha, const ComplexMatrix& x,int kx,
                                    const ComplexMatrix& y,int ky);
    void geru(complex<double> alpha, const ComplexMatrix& x,int kx,
                                     const ComplexMatrix& y,int ky);
    void gerc(complex<double> alpha, const ComplexMatrix& x,int kx,
                                     const ComplexMatrix& y,int ky);

    // symmetric rank-1 update
    void her(char uplo, complex<double> alpha, 
             const ComplexMatrix& x, int ix, char rowcol);
    
    // copy data in place between overlapping contexts
    void copyInPlace(ComplexMatrix& a);

    // get submatrix A(ia:ia+m,ja:ja+n) of A
    void getsub(const ComplexMatrix& a, int m, int n, int ia, int ja);

    // get submatrix A(ia:ia+m,ja:ja+n) of A and store in
    // this(idest:idest+m,jdest:jdest+n)
    void getsub(const ComplexMatrix& a,int m,int n,int isrc,int jsrc,
                int idest, int jdest);
    
    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(char transa, char transb,
         complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, complex<double> beta);
    
    // hermitian_matrix * matrix
    // *this = alpha * A * B + beta * this
    void hemm(char side, char uplo,
         complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, complex<double> beta);
    
    // complex_symmetric_matrix * matrix
    // *this = alpha * A * B + beta * this
    void symm(char side, char uplo,
         complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, complex<double> beta);
    
    // hermitian rank k update
    // this = beta * this + alpha * A * A^T  (trans=='n')
    // this = beta * this + alpha * A^T * A  (trans=='t')
    void herk(char uplo, char trans,
         complex<double> alpha, const ComplexMatrix& a, 
         complex<double> beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(complex<double> alpha, const ComplexMatrix& a, 
                   complex<double> beta);
    void transpose(const ComplexMatrix& a);

    void trmm(char side, char uplo, char trans, char diag, 
              complex<double> alpha, const ComplexMatrix& a);
    void trsm(char side, char uplo, char trans, char diag, 
              complex<double> alpha, const ComplexMatrix& a);
    void trtrs(char uplo, char trans, char diag, ComplexMatrix& b) const;
    
    // Cholesky decomposition of a hermitian matrix
    void potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    void potri(char uplo);
    
    // compute inverse of a square matrix
    void inverse(void);
    
    // Inverse of triangular matrix
    void trtri(char uplo,char diag);
    
    double pocon(char) const;
    
    // compute eigenvalues and eigenvectors of hermitian matrix *this
    void heev(char uplo, valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues and eigenvectors of generalized (non-symmetric) matrix *this
    void hegv(char uplo, valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues and eigenvectors of hermitian matrix *this, using divide-and-conquer
    void heevd(char uplo, valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues and eigenvectors of hermitian matrix *this, using MRRR
    void heevr(char uplo, valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues and eigenvectors of hermitian matrix *this, using pzheevx
    void heevx(char uplo, valarray<double>& w, ComplexMatrix& z);
    void heevx(char uplo, valarray<double>& w, ComplexMatrix& z, double tolerance);
    // compute eigenvalues (only) of hermitian matrix *this
    void heev(char uplo, valarray<double>& w);
    // compute eigenvalues (only) of generalized (non-symmetric) matrix *this
    void hegv(char uplo, valarray<double>& w);
    // compute eigenvalues (only) of hermitian matrix *this, using divide-and-conquer
    void heevd(char uplo, valarray<double>& w);
};
ostream& operator << ( ostream& os, const ComplexMatrix& a );
#endif
