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
// FourierTransform.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "FourierTransform.h"
#include "Basis.h"
#include "blas.h"

#include <complex>
#include <algorithm>
#include <map>
#include <cassert>

#if _OPENMP
#include <omp.h>
#else
//#if defined(HAVE_FFTW3_THREADS)
//#error "Need OpenMP to use FFTW3 threads"
//#endif
#endif

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#if defined(HAVE_FFTW2) || defined(HAVE_FFTW3)
#ifdef ADD_
#define zdscal zdscal_
#define zcopy zcopy_
#endif
#endif

#if defined(HAVE_FFTW2)
#if defined(FFTWMEASURE)
#define FFTW_ALGO FFTW_MEASURE
#else
#define FFTW_ALGO FFTW_ESTIMATE
#endif
#endif

#if defined(HAVE_FFTW3)
#if defined(FFTWMEASURE)
#define FFTW_ALGO ( FFTW_MEASURE | FFTW_UNALIGNED )
#else
#define FFTW_ALGO ( FFTW_ESTIMATE | FFTW_UNALIGNED )
#endif
#endif

#if defined(HAVE_FFTW2) || defined(HAVE_FFTW3)
extern "C" void zdscal(int *n,double *alpha,std::complex<double> *x,int *incx);
#elif USE_ESSL_FFT
extern "C" {
  void dcft_(int *initflag, std::complex<double> *x, int *inc2x, int *inc3x,
             std::complex<double> *y, int *inc2y, int *inc3y,
             int *length, int *ntrans, int *isign,
             double *scale, double *aux1, int *naux1,
             double *aux2, int *naux2);
  void dcft2_(int *initflag, std::complex<double> *x, int *inc1x, int *inc2x,
             std::complex<double> *y, int *inc1y, int *inc2y,
             int *n1, int *n2, int *isign,
             double *scale, double *aux1, int *naux1,
             double *aux2, int *naux2);
#define USE_GATHER_SCATTER 1
}
#elif defined(FFT_NOLIB)
void cfftm ( std::complex<double> *ain, std::complex<double> *aout,
  double scale, int ntrans, int length, int ainc, int ajmp, int idir );
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if USE_GATHER_SCATTER
extern "C" {
  // zgthr: x(i) = y(indx(i))
  void zgthr_(int* n, std::complex<double>* y,
              std::complex<double>* x, int*indx);
  // zsctr: y(indx(i)) = x(i)
  void zsctr_(int* n, std::complex<double>* x, int* indx,
              std::complex<double>* y);
}
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////////////
FourierTransform::~FourierTransform()
{
#if HAVE_FFTW2
  fftw_destroy_plan(fwplan0);
  fftw_destroy_plan(fwplan1);
  fftw_destroy_plan(fwplan2);
  fftw_destroy_plan(bwplan0);
  fftw_destroy_plan(bwplan1);
  fftw_destroy_plan(bwplan2);
#endif

#if HAVE_FFTW3
#if HAVE_FFTW3_THREADS
  fftw_cleanup_threads();
#endif
#if defined(HAVE_FFTW3_2D) || defined(HAVE_FFTW3_THREADS)
  fftw_destroy_plan(fwplan2d);
  fftw_destroy_plan(bwplan2d);
#else
  fftw_destroy_plan(fwplanx);
  fftw_destroy_plan(bwplanx);
  fftw_destroy_plan(fwplany);
  fftw_destroy_plan(bwplany);
#endif
  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);
#endif
}

////////////////////////////////////////////////////////////////////////////////
FourierTransform::FourierTransform (const Basis &basis,
  int np0, int np1, int np2) : comm_(basis.context().comm()), basis_(basis),
  np0_(np0), np1_(np1), np2_(np2)
{
  MPI_Comm_size(comm_,&nprocs_);
  MPI_Comm_rank(comm_,&myproc_);

  np2_loc_.resize(nprocs_);
  np2_first_.resize(nprocs_);

  // Block-cyclic distribution for np2
  // Partition np2 into nprocs_ intervals and
  // store local sizes in np2_loc_[iproc]
  // Use same block distribution as in ScaLAPACK
  // Blocks 0,...,nprocs_-2 have size np2_block_size
  // Block nprocs_-1 may have a smaller size
  if ( np2_ % nprocs_ == 0 )
  {
    // all blocks have equal size
    const int np2_block_size = np2_ / nprocs_;
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
      np2_loc_[iproc] = np2_block_size;
  }
  else
  {
    // first k-1 blocks have same size, k_th block is smaller, others zero
    const int np2_block_size = np2_ / nprocs_ + 1;
    const int k = np2_ / np2_block_size;
    for ( int iproc = 0; iproc < k; iproc++ )
      np2_loc_[iproc] = np2_block_size;
    np2_loc_[k] = np2_ - k * np2_block_size;
    for ( int iproc = k+1; iproc < nprocs_; iproc++ )
      np2_loc_[iproc] = 0;
  }

  np2_first_[0] = 0;
  for ( int iproc = 1; iproc < nprocs_; iproc++ )
  {
    np2_first_[iproc] = np2_first_[iproc-1] + np2_loc_[iproc-1];
  }

  // number of local z vectors
  if ( basis_.real() )
  {
    if ( myproc_ == 0 )
      // rod(0,0) is mapped to only one z vector
      nvec_ = 2 * basis_.nrod_loc() - 1;
    else
      nvec_ = 2 * basis_.nrod_loc();
  }
  else
  {
    nvec_ = basis_.nrod_loc();
  }

  // compute number of transforms along the x,y and z directions
  // ntrans0_ is the number of transforms along x in one of the two blocks
  // of vectors corresponding to positive and negative y indices
  ntrans0_ = max(abs(basis_.idxmax(1)),abs(basis_.idxmin(1)))+1;
  ntrans1_ = np0_;
  ntrans2_ = nvec_;

  // resize array zvec holding columns
  zvec_.resize(nvec_ * np2_);

#if TIMING
  tm_init.start();
#endif
  // Initialize FT library auxiliary arrays
  init_lib();
#if TIMING
  tm_init.stop();
#endif

  // allocate send buffer
  sbuf.resize(nvec_ * np2_);

  // allocate receive buffer
  if ( basis_.real() )
    rbuf.resize((2 * basis_.nrods() - 1) * np2_loc_[myproc_]);
  else
    rbuf.resize(basis_.nrods() * np2_loc_[myproc_]);

  // compute send/receive counts and displacements in units of sizeof(double)

  scounts.resize(nprocs_);
  sdispl.resize(nprocs_);
  rcounts.resize(nprocs_);
  rdispl.resize(nprocs_);

  if ( basis_.real() )
  {
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
    {
      scounts[iproc] = 2 * nvec_ * np2_loc_[iproc];
      int nvec_iproc = iproc == 0 ? 2*basis_.nrod_loc(iproc)-1 :
                                2 * basis_.nrod_loc(iproc);
      rcounts[iproc] = 2 * nvec_iproc * np2_loc_[myproc_];
    }
  }
  else
  {
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
    {
      scounts[iproc] = 2 * nvec_ * np2_loc_[iproc];
      int nvec_iproc = basis_.nrod_loc(iproc);
      rcounts[iproc] = 2 * nvec_iproc * np2_loc_[myproc_];
    }
  }

  sdispl[0] = 0;
  rdispl[0] = 0;
  for ( int iproc = 1; iproc < nprocs_; iproc++ )
  {
    sdispl[iproc] = sdispl[iproc-1] + scounts[iproc-1];
    rdispl[iproc] = rdispl[iproc-1] + rcounts[iproc-1];
  }

  if ( basis_.real() )
  {
    // compute index arrays ifftp_ and ifftm_ for mapping vector->zvec
    ifftp_.resize(basis_.localsize());
    ifftm_.resize(basis_.localsize());

    if ( myproc_ == 0 )
    {
      // this process holds rod(0,0)
      // nvec_ == 2 * nrod_loc - 1

      // map rod(0,0)
      // the positive segment of rod(0,0) maps onto the first half of
      // the first column of zvec_, and the negative segment maps onto
      // the second half
      int ig = 0;
      ifftp_[0] = 0;
      ifftm_[0] = 0;
      ig++;
      for ( int l = 1; l < basis_.rod_size(0); l++ )
      {
        ifftp_[ig] = l;
        ifftm_[ig] = np2_ - l;
        ig++;
      }

      // map other rods(h,k) on pe 0, h!=0, k!=0
      // rod(h,k) maps onto column (2*irod-1)*np2_ of zvec_, irod=1,..,nrods-1
      // rod(-h,-k) maps onto column (2*irod)*np2_ of zvec_, irod=1,..,nrods-1
      for ( int irod = 1; irod < basis_.nrod_loc(); irod++ )
      {
        const int rodsize = basis_.rod_size(irod);
        for ( int i = 0; i < rodsize; i++ )
        {
          const int l = i + basis_.rod_lmin(irod);
          int izp =  l;
          int izm = -l;
          if ( izp < 0 ) izp += np2_;
          if ( izm < 0 ) izm += np2_;
          ifftp_[ig] = ( 2 * irod - 1 ) * np2_ + izp;
          ifftm_[ig] = ( 2 * irod ) * np2_ + izm;
          ig++;
        }
      }
      assert(ig == basis_.localsize());
    }
    else
    {
      // this process does not hold rod(0,0)
      // map rods(h,k)
      // rod(h,k)   maps onto column (2*irod)*np2_ of zvec_, irod=0,..,nrods-1
      // rod(-h,-k) maps onto column (2*irod+1)*np2_ of zvec_, irod=0,..,nrods-1
      int ig = 0;
      for ( int irod = 0; irod < basis_.nrod_loc(); irod++ )
      {
        const int rodsize = basis_.rod_size(irod);
        for ( int i = 0; i < rodsize; i++ )
        {
          const int l = i + basis_.rod_lmin(irod);
          int izp =  l;
          int izm = -l;
          if ( izp < 0 ) izp += np2_;
          if ( izm < 0 ) izm += np2_;
          ifftp_[ig] = ( 2 * irod ) * np2_ + izp;
          ifftm_[ig] = ( 2 * irod + 1 ) * np2_ + izm;
          ig++;
        }
      }
      assert(ig == basis_.localsize());
    }

    // compute ipack index array
    // used in packing zvec_ into sbuf
    // sbuf[ipack_[i]] = zvec_[i]
    ipack_.resize(nvec_*np2_);
    int idest = 0;
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
    {
      int isource = np2_first_[iproc];
      int sz = np2_loc_[iproc];
      for ( int ivec = 0; ivec < nvec_; ivec++ )
      {
        for ( int i = 0; i < sz; i++ )
        {
          ipack_[isource+i] = idest + i;
        }
        idest += sz;
        isource += np2_;
      }
    }

    // compute array iunpack
    // used in unpacking rbuf into val
    // val[iunpack[i]] = rbuf[i]

    // rbuf contains 2*_nrods-1 segments of size np2_loc[myproc]
    // the position of vector ivec in local rbuf[_nrods*np2_loc_] is
    // obtained from rod_h[iproc][irod], rod_k[irod][iproc]
    // compute iunpack[i], i = 0, .. , _nrods * np2_loc_
    iunpack_.resize((2*basis_.nrods()-1)*np2_loc_[myproc_]);

    // map rod(0,0)
    for ( int l = 0; l < np2_loc_[myproc_]; l++ )
    {
      iunpack_[l] = l * np0_ * np1_;
    }
    int isource_p = np2_loc_[myproc_];
    int isource_m = 2 * np2_loc_[myproc_];

    // all rods of pe 0
    for ( int irod = 1; irod < basis_.nrod_loc(0); irod++ )
    {
      // map rod(h,k) and rod(-h,-k) columns of zvec_

      // map rod(h,k)
      // find position of rod(h,k) in the slab
      int hp = basis_.rod_h(0,irod);
      int kp = basis_.rod_k(0,irod);
      if ( hp < 0 ) hp += np0_;
      if ( kp < 0 ) kp += np1_;

      int hm = -hp;
      int km = -kp;
      if ( hm < 0 ) hm += np0_;
      if ( km < 0 ) km += np1_;

      for ( int l = 0; l < np2_loc_[myproc_]; l++ )
      {
        int idest_p = hp + np0_ * ( kp + np1_ * l );
        iunpack_[isource_p+l] = idest_p;

        int idest_m = hm + np0_ * ( km + np1_ * l );
        iunpack_[isource_m+l] = idest_m;
      }
      isource_p += 2 * np2_loc_[myproc_];
      isource_m += 2 * np2_loc_[myproc_];
    }

    // pe's above pe0
    for ( int iproc = 1; iproc < nprocs_; iproc++ )
    {
      for ( int irod = 0; irod < basis_.nrod_loc(iproc); irod++ )
      {
        // map rod(h,k) and rod(-h,-k) columns of zvec_

        // map rod(h,k)
        // find position of rod(h,k) in the slab
        int hp = basis_.rod_h(iproc,irod);
        int kp = basis_.rod_k(iproc,irod);
        if ( hp < 0 ) hp += np0_;
        if ( kp < 0 ) kp += np1_;

        int hm = -hp;
        int km = -kp;
        if ( hm < 0 ) hm += np0_;
        if ( km < 0 ) km += np1_;

        for ( int l = 0; l < np2_loc_[myproc_]; l++ )
        {
          int idest_p = hp + np0_ * ( kp + np1_ * l );
          iunpack_[isource_p+l] = idest_p;

          int idest_m = hm + np0_ * ( km + np1_ * l );
          iunpack_[isource_m+l] = idest_m;
        }
        isource_p += 2 * np2_loc_[myproc_];
        isource_m += 2 * np2_loc_[myproc_];
      }
    }
  }
  else
  {
    // basis is complex
    // compute index array ifftp_ for mapping vector->zvec
    // Note: ifftm_ is not used
    ifftp_.resize(basis_.localsize());

    // map rods(h,k)
    // rod(h,k)   maps onto column irod*np2_ of zvec_, irod=0,..,nrods-1
    int ig = 0;
    for ( int irod = 0; irod < basis_.nrod_loc(); irod++ )
    {
      const int rodsize = basis_.rod_size(irod);
      for ( int i = 0; i < rodsize; i++ )
      {
        const int l = i + basis_.rod_lmin(irod);
        int iz =  l;
        if ( iz < 0 ) iz += np2_;
        ifftp_[ig] = irod * np2_ + iz;
        ig++;
      }
    }
    assert(ig == basis_.localsize());

    // compute ipack index array
    // used in packing zvec_ into sbuf
    // sbuf[ipack_[i]] = zvec_[i]
    ipack_.resize(nvec_*np2_);
    int idest = 0;
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
    {
      int isource = np2_first_[iproc];
      int sz = np2_loc_[iproc];
      for ( int ivec = 0; ivec < nvec_; ivec++ )
      {
        for ( int i = 0; i < sz; i++ )
        {
          ipack_[isource+i] = idest + i;
        }
        idest += sz;
        isource += np2_;
      }
    }

    // compute array iunpack
    // used in unpacking rbuf into val
    // val[iunpack[i]] = rbuf[i]

    // rbuf contains _nrods segments of size np2_loc[mype]
    // the position of vector ivec in local rbuf[_nrods*np2_loc_] is
    // obtained from rod_h[iproc][irod], rod_k[irod][iproc]
    // compute iunpack[i], i = 0, .. , _nrods * np2_loc_
    iunpack_.resize(basis_.nrods()*np2_loc_[myproc_]);

    int isource = 0;
    for ( int iproc = 0; iproc < nprocs_; iproc++ )
    {
      for ( int irod = 0; irod < basis_.nrod_loc(iproc); irod++ )
      {
        // map rod(h,k)
        // find position of rod(h,k) in the slab
        int h = basis_.rod_h(iproc,irod);
        int k = basis_.rod_k(iproc,irod);
        if ( h < 0 ) h += np0_;
        if ( k < 0 ) k += np1_;

        for ( int l = 0; l < np2_loc_[myproc_]; l++ )
        {
          int idest = h + np0_ * ( k + np1_ * l );
          iunpack_[isource+l] = idest;

        }
        isource += np2_loc_[myproc_];
      }
    }
  }

  // for ( int ig = 0; ig < basis_.localsize(); ig++ )
  // {
  //   assert(ifftp_[ig] >= 0 && ifftp_[ig] < zvec_.size());
  //   assert(ifftm_[ig] >= 0 && ifftm_[ig] < zvec_.size());
  // }

#if USE_GATHER_SCATTER
  // shift index array by one for fortran ZGTHR and ZSCTR calls
  for ( int i = 0; i < iunpack_.size(); i++ )
  {
    iunpack_[i]++;
  }
  for ( int i = 0; i < ipack_.size(); i++ )
  {
    ipack_[i]++;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::backward(const complex<double>* c, complex<double>* f)
{
#if TIMING
  tm_b_map.start();
#endif
  vector_to_zvec(c);
#if TIMING
  tm_b_map.stop();
#endif
  bwd(f);
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::forward(complex<double>* f, complex<double>* c)
{
  fwd(f);
#if TIMING
  tm_f_map.start();
#endif
  zvec_to_vector(c);
#if TIMING
  tm_f_map.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::backward(const complex<double>* c1,
                               const complex<double>* c2,
                               complex<double>* f)
{
#if TIMING
  tm_b_map.start();
#endif
  doublevector_to_zvec(c1,c2);
#if TIMING
  tm_b_map.stop();
#endif
  bwd(f);
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::forward(complex<double>* f,
  complex<double>* c1, complex<double>* c2)
{
  fwd(f);
#if TIMING
  tm_f_map.start();
#endif
  zvec_to_doublevector(c1,c2);
#if TIMING
  tm_f_map.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::bwd(complex<double>* val)
{
  // transform zvec along z, transpose and transform along x,y, store
  // result in val
  // The columns of zvec[nvec_ * np2_] contain the full vectors
  // to be transformed
  //
  // If the basis is real: Column (h,k) is followed by column (-h,-k),
  // except for (0,0)

#if TIMING
  tm_b_fft.start();
  tm_b_z.start();
#endif

#if USE_ESSL_FFT
  int inc1 = 1, inc2 = np2_, ntrans = nvec_, isign = -1, initflag = 0;
  double scale = 1.0;
  if (ntrans > 0)
     dcft_(&initflag,&zvec_[0],&inc1,&inc2,&zvec_[0],&inc1,&inc2,&np2_,&ntrans,
           &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
#elif HAVE_FFTW2
   /*
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
#if _OPENMP
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_one(bwplan2,(FFTW_COMPLEX*)&zvec_[i*np2_],(FFTW_COMPLEX*)0);
  }
#else
  int ntrans = nvec_, inc1 = 1, inc2 = np2_;
  fftw(bwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP

#elif HAVE_FFTW3 // HAVE_FFTW2

#if HAVE_FFTW3_THREADS
  fftw_execute_dft ( bwplan, (fftw_complex*)&zvec_[0],
                     (fftw_complex*)&zvec_[0]);
#else
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_execute_dft ( bwplan, (fftw_complex*)&zvec_[i*np2_],
                       (fftw_complex*)&zvec_[i*np2_]);
  }
#endif // HAVE_FFTW3_THREADS

#elif defined(FFT_NOLIB) // HAVE_FFTW3
  // No library
  /* Transform along z */
  int ntrans = nvec_;
  int length = np2_;
  int ainc   = 1;
  int ajmp   = np2_;
  double scale = 1.0;
  int idir = -1;
  cfftm ( &zvec_[0], &zvec_[0], scale, ntrans, length, ainc, ajmp, idir );
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif // HAVE_FFTW3

#if TIMING
  tm_b_z.stop();
  tm_b_com.start();
  tm_b_fft.stop();
  tm_b_pack.start();
#endif

  // scatter zvec_ to sbuf for transpose
#if USE_GATHER_SCATTER
  // zsctr: y(indx(i)) = x(i)
  // void zsctr_(int* n, complex<double>* x, int* indx, complex<double>* y);
  {
    complex<double>* y = &sbuf[0];
    complex<double>* x = &zvec_[0];
    int n = zvec_.size();
    zsctr_(&n,x,&ipack_[0],y);
  }
#else
  const int zvec_size = zvec_.size();
  double* const ps = (double*) &sbuf[0];
  const double* const pz = (double*) &zvec_[0];
  #pragma omp parallel for
  for ( int i = 0; i < zvec_size; i++ )
  {
    // sbuf[ipack_[i]] = zvec_[i];
    const int ip = ipack_[i];
    const double a = pz[2*i];
    const double b = pz[2*i+1];
    ps[2*ip]   = a;
    ps[2*ip+1] = b;
  }
#endif

  // segments of z-vectors are now in sbuf

#if TIMING
  tm_b_pack.stop();
  tm_b_mpi.start();
#endif

  // transpose
#if USE_MPI
  int status = MPI_Alltoallv((double*)&sbuf[0],&scounts[0],&sdispl[0],
      MPI_DOUBLE,(double*)&rbuf[0],&rcounts[0],&rdispl[0],MPI_DOUBLE,
      comm_);
  if ( status != 0 )
  {
    cout << " FourierTransform: status = " << status << endl;
    MPI_Abort(MPI_COMM_WORLD,2);
  }
#else
  assert(sbuf.size()==rbuf.size());
  rbuf = sbuf;
#endif

#if TIMING
  tm_b_mpi.stop();
  tm_b_zero.start();
#endif

  // copy from rbuf to val
  // scatter index array iunpack
  {
    const int len = np012loc() * 2;
    double* const pv = (double*) &val[0];
    #pragma omp parallel for
    for ( int i = 0; i < len; i++ )
    {
      pv[i]   = 0.0;
    }
  }

#if TIMING
  tm_b_zero.stop();
  tm_b_unpack.start();
#endif

#if USE_GATHER_SCATTER
  // zsctr(n,x,indx,y): y(indx(i)) = x(i)
  {
    complex<double>* y = &val[0];
    complex<double>* x = &rbuf[0];
    int n = rbuf.size();
    zsctr_(&n,x,&iunpack_[0],y);
  }
#else
  {
    const int rbuf_size = rbuf.size();
    const double* const pr = (double*) &rbuf[0];
    double* const pv = (double*) &val[0];
    #pragma omp parallel for
    for ( int i = 0; i < rbuf_size; i++ )
    {
      // val[iunpack_[i]] = rbuf[i];
      const int iu = iunpack_[i];
      const double a = pr[2*i];
      const double b = pr[2*i+1];
      pv[2*iu]   = a;
      pv[2*iu+1] = b;
    }
  }
#endif

#if TIMING
  tm_b_unpack.stop();
  tm_b_fft.start();
  tm_b_com.stop();
  tm_b_xy.start();
#endif

#if HAVE_FFTW3
#if HAVE_FFTW3_THREADS
  fftw_execute_dft ( bwplan2d, (fftw_complex*)&val[0],
                     (fftw_complex*)&val[0] );
#elif HAVE_FFTW3_2D
  #pragma omp parallel for
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
    fftw_execute_dft ( bwplan2d, (fftw_complex*)&val[k*np0_*np1_],
                       (fftw_complex*)&val[k*np0_*np1_] );
#else // FFTW3_2D
  // fftw3 1d
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    int ibase = k * np0_ * np1_;
#if TIMING
    tm_b_x.start();
#endif
    #pragma omp parallel for
    for ( int i = 0; i < ntrans0_; i++ )
    {
      // Transform first block along x: positive y indices
      fftw_execute_dft ( bwplanx, (fftw_complex*)&val[ibase+i*np0_],
                         (fftw_complex*)&val[ibase+i*np0_]);
      // Transform second block along x: negative y indices
      fftw_execute_dft ( bwplanx,
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_]);
    }
#if TIMING
    tm_b_x.stop();
    tm_b_y.start();
#endif
#if FFTW_TRANSPOSE
    #pragma omp parallel
    {
      vector<complex<double> >t_trans(np1_);
      #pragma omp for
      for ( int i = 0; i < np0_; i++ )
      {
        int length = t_trans.size();
        int inc1 = 1, inc2 = np0_;
        zcopy(&length, &val[ibase+i], &inc2, &t_trans[0], &inc1);
        fftw_execute_dft ( bwplany, (fftw_complex*)&t_trans[0],
                           (fftw_complex*)&t_trans[0]);
        zcopy(&length, &t_trans[0], &inc1, &val[ibase+i], &inc2);
      }
    }
#else // FFTW_TRANSPOSE
    #pragma omp parallel for
    for ( int i = 0; i < np0_; i++ )
    {
      fftw_execute_dft ( bwplany, (fftw_complex*)&val[ibase+i],
                         (fftw_complex*)&val[ibase+i]);
    }
#endif // FFTW_TRANSPOSE
#if TIMING
    tm_b_y.stop();
#endif
  }
#endif // HAVE_FFTW3_2D

#elif USE_ESSL_FFT
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if USE_ESSL_2DFFT

    // use 2D FFT for x and y transforms
    int inc1, inc2, istart, isign = -1, initflag = 0;
    double scale = 1.0;

    // xy transform
    istart = k * np0_ * np1_;
    inc1 = 1; inc2 = np0_;
    dcft2_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
          &np0_,&np1_,&isign,&scale,&aux1xyb[0],&naux1xy,&aux2[0],&naux2);

#else

    // use multiple 1-d FFTs for x and y transforms

    int inc1, inc2, ntrans, istart, length, isign = -1, initflag = 0;
    double scale = 1.0;
    // transform only non-zero vectors along x
    // First block: positive y indices: [0,ntrans0_]
    ntrans = ntrans0_;
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    length = np0_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);

    // Second block: negative y indices: [np1-ntrans0_,np1-1]
    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    length = np0_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);

    // transform along y for all values of x
    ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    length = np1_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1yb[0],&naux1y,&aux2[0],&naux2);
#endif // USE_ESSL_2DFFT
  } // k

#elif HAVE_FFTW2
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if _OPENMP
  int ibase = k * np0_ * np1_;
  #pragma omp parallel for
  for ( int i = 0; i < ntrans0_; i++ )
  {
    //#pragma omp task
    {
      // Transform first block along x: positive y indices
      fftw_one(bwplan0,(FFTW_COMPLEX*)&val[ibase+i*np0_],(FFTW_COMPLEX*)0);
      //fftw(bwplan0,1,(FFTW_COMPLEX*)&val[ibase+i*np0_],1,np0_,
      //               (FFTW_COMPLEX*)0,0,0);
      // Transform second block along x: negative y indices
      fftw_one(bwplan0,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                       (FFTW_COMPLEX*)0);
      //fftw(bwplan0,1,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],1,np0_,
      //               (FFTW_COMPLEX*)0,0,0);
    }
  }

  //complex<double> *tmp1 = new complex<double>[np1_];
  #pragma omp parallel for
  for ( int i = 0; i < np0_; i++ )
  {
    {
      // transform along y for all values of x
      // copy data to local array
      int one=1;
      #if 0
      zcopy_(&np1_,&val[ibase+i],&np0_,tmp1,&one);
      fftw_one(bwplan1,(FFTW_COMPLEX*)tmp1,(FFTW_COMPLEX*)0);
      zcopy_(&np1_,tmp1,&one,&val[ibase+i],&np0_);
      #else
      fftw(bwplan1,1,(FFTW_COMPLEX*)&val[ibase+i],np0_,one,
                     (FFTW_COMPLEX*)0,0,0);
      #endif
    }
  }
  //delete [] tmp1;
#else // _OPENMP
    int inc1, inc2, istart;

    int ntrans = ntrans0_;
    // Transform first block along x: positive y indices
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    fftw(bwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
    // Transform second block along x: negative y indices
    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    fftw(bwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);

    // transform along y for all values of x
    ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    fftw(bwplan1,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP
  } // k
#elif defined(FFT_NOLIB) // HAVE_FFTW2
  // No library
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
    // transform along x for non-zero elements
    // Transform first block along x: positive y indices
    int ntrans = ntrans0_;
    int istart = k * np0_ * np1_;
    int length = np0_;
    int ainc   = 1;
    int ajmp   = np0_;
    double scale = 1.0;
    int idir = -1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // Transform second block along x: negative y indices
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // transform along y for all values of x
    ntrans = np0_;
    istart = k * np0_ * np1_;
    length = np1_;
    ainc = np0_;
    ajmp = 1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );
  } // for k
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_b_xy.stop();
  tm_b_fft.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fwd(complex<double>* val)
{
#if TIMING
  tm_f_fft.start();
  tm_f_xy.start();
#endif

//fftw_execute_dft is thread safe
#if HAVE_FFTW3
#if HAVE_FFTW3_THREADS
  fftw_execute_dft ( fwplan2d, (fftw_complex*)&val[0],
                     (fftw_complex*)&val[0] );
#elif HAVE_FFTW3_2D // HAVE_FFTW3_2D
  #pragma omp parallel for
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
    fftw_execute_dft ( fwplan2d, (fftw_complex*)&val[k*np0_*np1_],
                       (fftw_complex*)&val[k*np0_*np1_] );
#else // HAVE_FFTW3_2D
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    const int ibase = k * np0_ * np1_;
#if TIMING
    tm_f_y.start();
#endif
#if FFTW_TRANSPOSE
    #pragma omp parallel
    {
      vector<complex<double> >t_trans(np1_);
      #pragma omp for
      for ( int i = 0; i < np0_; i++ )
      {
        int length = t_trans.size();
        int inc1 = 1, inc2 = np0_;
        zcopy(&length, &val[ibase+i], &inc2, &t_trans[0], &inc1);
        fftw_execute_dft ( fwplany, (fftw_complex*)&t_trans[0],
                         (fftw_complex*)&t_trans[0]);
        zcopy(&length, &t_trans[0], &inc1, &val[ibase+i], &inc2);
      }
    }
#else // FFTW_TRANSPOSE
    #pragma omp parallel for
    for ( int i = 0; i < np0_; i++ )
    {
      fftw_execute_dft ( fwplany, (fftw_complex*)&val[ibase+i],
                         (fftw_complex*)&val[ibase+i]);
    }
#endif // FFTW_TRANSPOSE
#if TIMING
    tm_f_y.stop();
    tm_f_x.start();
#endif
    #pragma omp parallel for
    for ( int i = 0; i < ntrans0_; i++ )
    {
      // Transform first block along x: positive y indices
      fftw_execute_dft ( fwplanx,(fftw_complex*)&val[ibase+i*np0_],
                         (fftw_complex*)&val[ibase+i*np0_]);

      // Transform second block along x: negative y indices
      fftw_execute_dft ( fwplanx,
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_]);
    }
#if TIMING
    tm_f_x.stop();
#endif
  }
#endif // HAVE_FFTW3_2D
#elif USE_ESSL_FFT
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if USE_ESSL_2DFFT

    // use 2D FFT for x and y transforms
    int inc1, inc2, istart, isign = 1, initflag = 0;
    double scale = 1.0;

    // xy transform
    istart = k * np0_ * np1_;
    inc1 = 1; inc2 = np0_;
    dcft2_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
          &np0_,&np1_,&isign,&scale,&aux1xyf[0],&naux1xy,&aux2[0],&naux2);

#else

    // use multiple 1-d FFTs for x and y transforms

    int inc1, inc2, ntrans, istart, length, isign = 1, initflag = 0;
    double scale = 1.0;
    // transform along y for all values of x
    ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    length = np1_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1yf[0],&naux1y,&aux2[0],&naux2);

    // transform only non-zero vectors along x
    ntrans = ntrans0_;
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    length = np0_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);

    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    length = np0_;
    dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
         &length,&ntrans,&isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);
#endif // USE_ESSL_2DFFT
  } // k
#elif HAVE_FFTW2
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if _OPENMP
  int ibase = k * np0_ * np1_;
  //complex<double> *tmp1 = new complex<double>[np1_];
  #pragma omp parallel for
  for ( int i = 0; i < np0_; i++ )
  {
    //#pragma omp task
    {
      // transform along y for all values of x
      // copy data to local array
      int one=1;
      #if 0
      zcopy_(&np1_,&val[ibase+i],&np0_,tmp1,&one);
      fftw_one(fwplan1,(FFTW_COMPLEX*)tmp1,(FFTW_COMPLEX*)0);
      zcopy_(&np1_,tmp1,&one,&val[ibase+i],&np0_);
      #else
      fftw(fwplan1,1,(FFTW_COMPLEX*)&val[ibase+i],np0_,one,
                     (FFTW_COMPLEX*)0,0,0);
      #endif
    }
  }
  //delete [] tmp1;

  #pragma omp parallel for
  for ( int i = 0; i < ntrans0_; i++ )
  {
    //#pragma omp task
    {
      // Transform first block along x: positive y indices
      fftw_one(fwplan0,(FFTW_COMPLEX*)&val[ibase+i*np0_],(FFTW_COMPLEX*)0);
      // Transform second block along x: negative y indices
      fftw_one(fwplan0,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                       (FFTW_COMPLEX*)0);
    }
  }
#else // _OPENMP
    int inc1, inc2, istart;

    // transform along y for all values of x
    int ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    fftw(fwplan1,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);

    ntrans = ntrans0_;
    // Transform first block along x: positive y indices
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    fftw(fwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
    // Transform second block along x: negative y indices
    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    fftw(fwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP
  } // k
#elif defined(FFT_NOLIB)
  // No library
  for ( int k = 0; k < np2_loc_[myproc_]; k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
    // transform along y for all values of x
    int ntrans = np0_;
    int istart = k * np0_ * np1_;
    int length = np1_;
    int ainc = np0_;
    int ajmp = 1;
    double scale = 1.0;
    int idir = 1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // transform along x for non-zero elements
    ntrans = ntrans0_;
    istart = k * np0_ * np1_;
    length = np0_;
    ainc   = 1;
    ajmp   = np0_;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );
  } // for k
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_f_xy.stop();
  tm_f_com.start();
  tm_f_fft.stop();
  tm_f_pack.start();
#endif

  // gather val into rbuf
#if USE_GATHER_SCATTER
  // zgthr: x(i) = y(indx(i))
  // void zgthr_(int* n, complex<double>* y, complex<double>* x, int*indx);
  {
    complex<double>* y = &val[0];
    complex<double>* x = &rbuf[0];
    int n = rbuf.size();
    zgthr_(&n,y,x,&iunpack_[0]);
  }
#else
  const int rbuf_size = rbuf.size();
  double* const pr = (double*) &rbuf[0];
  const double* const pv = (double*) &val[0];
  #pragma omp parallel for
  for ( int i = 0; i < rbuf_size; i++ )
  {
    // rbuf[i] = val[iunpack_[i]];
    const int iu = iunpack_[i];
    const double a = pv[2*iu];
    const double b = pv[2*iu+1];
    pr[2*i]   = a;
    pr[2*i+1] = b;
  }
#endif

  // transpose
#if TIMING
  tm_f_pack.stop();
  tm_f_mpi.start();
#endif
#if USE_MPI
  int status = MPI_Alltoallv((double*)&rbuf[0],&rcounts[0],&rdispl[0],
      MPI_DOUBLE,(double*)&sbuf[0],&scounts[0],&sdispl[0],MPI_DOUBLE,
      comm_);
  assert ( status == 0 );
#else
  assert(sbuf.size()==rbuf.size());
  sbuf = rbuf;
#endif

  // segments of z-vectors are now in sbuf
  // gather sbuf into zvec_
#if TIMING
  tm_f_mpi.stop();
  tm_f_unpack.start();
#endif

#if USE_GATHER_SCATTER
  // zgthr: x(i) = y(indx(i))
  // void zgthr_(int* n, complex<double>* y, complex<double>* x, int*indx);
  {
    complex<double>* y = &sbuf[0];
    complex<double>* x = &zvec_[0];
    int n = zvec_.size();
    zgthr_(&n,y,x,&ipack_[0]);
  }
#else // no gather scatter
  const int zvec_size = zvec_.size();
  const double* const ps = (double*) &sbuf[0];
  double* const pz = (double*) &zvec_[0];
  #pragma omp parallel for
  for ( int i = 0; i < zvec_size; i++ )
  {
    // zvec_[i] = sbuf[ipack_[i]];
    const int ip = ipack_[i];
    const double a = ps[2*ip];
    const double b = ps[2*ip+1];
    pz[2*i]   = a;
    pz[2*i+1] = b;
  }
#endif

  // transform along z
#if TIMING
  tm_f_unpack.stop();
  tm_f_fft.start();
  tm_f_com.stop();
  tm_f_z.start();
#endif

#if USE_ESSL_FFT
  int inc1 = 1, inc2 = np2_, ntrans = nvec_, isign = 1, initflag = 0;
  double scale = 1.0 / (np0_ * np1_ * np2_);

  if (ntrans > 0)
     dcft_(&initflag,&zvec_[0],&inc1,&inc2,&zvec_[0],&inc1,&inc2,&np2_,&ntrans,
           &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);

#elif HAVE_FFTW2
#if _OPENMP
  const double fac = 1.0 / ( np0_ * np1_ * np2_ );
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    //#pragma omp task
    fftw_one(fwplan2,(FFTW_COMPLEX*)&zvec_[i*np2_],(FFTW_COMPLEX*)0);
    for ( int j = 0; j < np2_; j++ )
      zvec_[i*np2_+j] *= fac;
  }
  // int inc1=1;
  // zdscal(&len,&fac,&zvec_[0],&inc1);
#else

 /*
  * void fftw(fftw_plan plan, int howmany,
  *    FFTW_COMPLEX *in, int istride, int idist,
  *    FFTW_COMPLEX *out, int ostride, int odist);
  */
  int ntrans, inc1, inc2;

  ntrans = nvec_;
  inc1 = 1;
  inc2 = np2_;
  fftw(fwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  int len = zvec_.size();
  double fac = 1.0 / ( np0_ * np1_ * np2_ );
  zdscal(&len,&fac,&zvec_[0],&inc1);
#endif
#elif HAVE_FFTW3

#if HAVE_FFTW3_THREADS
  fftw_execute_dft ( fwplan, (fftw_complex*)&zvec_[0],
                    (fftw_complex*)&zvec_[0]);
#else
  // do np2_ same for D_USE_1D or not
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_execute_dft ( fwplan, (fftw_complex*)&zvec_[i*np2_],
                      (fftw_complex*)&zvec_[i*np2_]);
  }
#endif
  // scale
  double fac = 1.0 / ( np0_ * np1_ * np2_ );
  int len = zvec_.size();
  int inc1 = 1;
  zdscal(&len,&fac,&zvec_[0],&inc1);
#elif defined(FFT_NOLIB)
  // No library
  /* Transform along z */
  int ntrans = nvec_;
  int length = np2_;
  int ainc   = 1;
  int ajmp   = np2_;
  double scale = 1.0 / ( np0_ * np1_ * np2_ );
  int idir = 1;
  cfftm ( &zvec_[0], &zvec_[0], scale, ntrans, length, ainc, ajmp, idir );
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_f_z.stop();
  tm_f_fft.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::init_lib(void)
{
  // initialization of FFT libs

   int nthreads = 1;
#if _OPENMP
   nthreads = omp_get_max_threads();
#endif  

#if USE_ESSL_FFT
  complex<double> *p = 0;
#if USE_ESSL_2DFFT
  // use 2D FFT for x and y transforms and 1D FFT for z transforms
  naux1xy = 40000 + 2.28 * (np0_+np1_);
  aux1xyf.resize(naux1xy);
  aux1xyb.resize(naux1xy);
  int r = max(np0_,np1_);
  int s = min(64,min(np0_,np1_));
  naux2 = 20000 + (2*r+256)*(s+2.28);

  naux1z = 20000 + 2.28 * np2_;
  aux1zf.resize(naux1z);
  aux1zb.resize(naux1z);

  int ntrans2 = nvec_;
  int naux2z = 20000 + 2.28 * np2_ + (256 + 2*np2_)*min(64,ntrans2);
  naux2 = max( naux2, naux2z );
  aux2.resize(naux2);

  double scale = 1.0;

  // initialize xy transforms
  int initflag = 1, inc1, inc2, isign = -1;
  inc1 = 1; inc2 = np0_;
  dcft2_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&np1_,
         &isign,&scale,&aux1xyb[0],&naux1xy,&aux2[0],&naux2);
  isign = 1;
  dcft2_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&np1_,
         &isign,&scale,&aux1xyf[0],&naux1xy,&aux2[0],&naux2);

  // initialize z transforms
  int ntrans = nvec_;
  inc1 = 1; inc2 = np2_;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
       &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
  isign = 1; scale = 1.0 / ( np0_ * np1_ * np2_ );
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
       &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);

  // 1d z transforms
  naux1z1d = (int) (20000 + 2.28 * np2_);
  aux1zf1d.resize(naux1z1d);
  aux1zb1d.resize(naux1z1d);
  scale = 1.0;
  inc1 = 1; inc2 = np2_; ntrans = 1;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zb1d[0],&naux1z1d,&aux2[0],&naux2);
  isign = 1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zf1d[0],&naux1z1d,&aux2[0],&naux2);

#else // USE_ESSL_2DFFT

  naux1x = (int) (20000 + 2.28 * np0_);
  naux1y = (int) (20000 + 2.28 * np1_);
  naux1z = (int) (20000 + 2.28 * np2_);
  naux1z1d = (int) (20000 + 2.28 * np2_);
  aux1xf.resize(naux1x);
  aux1yf.resize(naux1y);
  aux1zf.resize(naux1z);
  aux1zf1d.resize(naux1z1d);
  aux1xb.resize(naux1x);
  aux1yb.resize(naux1y);
  aux1zb.resize(naux1z);
  aux1zb1d.resize(naux1z1d);

  int naux2x = (int) (20000 + 2.28 * np0_ + (256 + 2*np0_)*min(64,ntrans0_));
  naux2 = naux2x;
  int naux2y = (int) (20000 + 2.28 * np1_ + (256 + 2*np1_)*min(64,ntrans1_));
  naux2 = max( naux2, naux2y );
  int naux2z = (int) (20000 + 2.28 * np2_ + (256 + 2*np2_)*min(64,ntrans2_));
  naux2 = max( naux2, naux2z );
  aux2.resize(naux2);

  // initialize x, y and z transforms

  int initflag = 1, inc1, inc2, ntrans, isign;
  double scale = 1.0;

  // x transforms
  inc1 = 1; inc2 = np0_; ntrans = ntrans0_;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&ntrans,
        &isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);
  isign = 1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&ntrans,
        &isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);

  // y transforms
  inc1 = np0_; inc2 = 1; ntrans = ntrans1_;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np1_,&ntrans,
        &isign,&scale,&aux1yb[0],&naux1y,&aux2[0],&naux2);
  isign = 1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np1_,&ntrans,
        &isign,&scale,&aux1yf[0],&naux1y,&aux2[0],&naux2);

  // z transforms
  inc1 = 1; inc2 = np2_; ntrans = ntrans2_;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
  isign = 1; scale = 1.0 / ( np0_ * np1_ * np2_ );
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);

  // 1d z transforms
  scale = 1.0;
  inc1 = 1; inc2 = np2_; ntrans = 1;
  isign = -1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zb1d[0],&naux1z1d,&aux2[0],&naux2);
  isign = 1;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
        &isign,&scale,&aux1zf1d[0],&naux1z1d,&aux2[0],&naux2);

#endif // USE_ESSL_2DFFT

#elif HAVE_FFTW2

  fwplan0 = fftw_create_plan(np0_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  fwplan1 = fftw_create_plan(np1_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  fwplan2 = fftw_create_plan(np2_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan0 = fftw_create_plan(np0_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan1 = fftw_create_plan(np1_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan2 = fftw_create_plan(np2_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);

#elif HAVE_FFTW3
  vector<complex<double> > aux(np0_*np1_);
  #if defined(HAVE_FFTW3MKL) && !defined(HAVE_FFTW3_THREADS) && _OPENMP
  fftw3_mkl.number_of_user_threads = nthreads;
  #endif


  // 1d fft plan for fft_1z
  zin_.resize(np2_);
  zout_.resize(np2_);
  fwplan1z = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zin_[0]),
                                (fftw_complex*)(&zout_[0]), -1,
                                FFTW_ALGO );
  bwplan1z = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zin_[0]),
                                (fftw_complex*)(&zout_[0]), 1,
                                FFTW_ALGO );
  
#if HAVE_FFTW3_THREADS
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
  vector<complex<double> > aux1(np0_*np1_*np2_loc_[myproc_]);

  // xy
  int rank = 2;
  int n[] = {np1_,np0_};
  int howmany = np2_loc_[myproc_];
  //int howmany = 1;
  int idist = np0_*np1_, odist = np0_*np1_;
  int istride = 1, ostride = 1; /* array is contiguous in memory */
  int *inembed = n, *onembed = n;

  fwplan2d = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux1[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux1[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplan2d = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux1[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux1[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);

  // z
  rank = 1;
  int nz[] = {np2_};
  howmany = nvec_;
  idist = np2_, odist = np2_;
  istride = 1, ostride = 1; /* array is contiguous in memory */
  inembed = nz, onembed = nz;

  fwplan = fftw_plan_many_dft(rank, nz, howmany, (fftw_complex*)&zvec_[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&zvec_[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplan = fftw_plan_many_dft(rank, nz, howmany, (fftw_complex*)&zvec_[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&zvec_[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);


#else // HAVE_FFTW3_THREADS
#if HAVE_FFTW3_2D
  // row major in FFTW3 2d plans
  fwplan2d = fftw_plan_dft_2d ( np1_, np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplan2d = fftw_plan_dft_2d ( np1_, np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );
#else // HAVE_FFTW3_2D
  // FFTW3 1D
  fwplanx = fftw_plan_dft_1d ( np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplanx = fftw_plan_dft_1d ( np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );

#if FFTW_TRANSPOSE
  fwplany = fftw_plan_dft_1d ( np1_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplany = fftw_plan_dft_1d ( np1_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );

#else // FFTW_TRANSPOSE
  // strided FFT
  int rank = 1;
  int n[] = {np1_};
  int howmany = 1;
  int idist = 1, odist = 1;
  int istride = np0_, ostride = np0_; /* array is contiguous in memory */
  int *inembed = n, *onembed = n;

  fwplany = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplany = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);
#endif // FFTW_TRANSPOSE
#endif // HAVE_FFTW3_2D
  // do z using 1d plans
  fwplan = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zvec_[0]),
                                (fftw_complex*)(&zvec_[0]), -1,
                                FFTW_ALGO );
  bwplan = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zvec_[0]),
                                (fftw_complex*)(&zvec_[0]), 1,
                                FFTW_ALGO );
#endif //HAVE_FFTW3_THREADS

#elif FFT_NOLIB // HAVE_FFTW3
  /* no library */
#else
#error "Must define HAVE_FFTW2, HAVE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::vector_to_zvec(const complex<double> *c)
{
  const int ng = basis_.localsize();
  const int zvec_size = zvec_.size();
  double* const pz = (double*) &zvec_[0];
  #pragma omp parallel for
  for ( int i = 0; i < zvec_size; i++ )
  {
    pz[2*i]   = 0.0;
    pz[2*i+1] = 0.0;
  }
  const double* const pc = (double*) &c[0];
  if ( basis_.real() )
  {
    #pragma omp parallel for
    for ( int ig = 0; ig < ng; ig++ )
    {
      // zvec_[ifftp_[ig]] = c[ig];
      // zvec_[ifftm_[ig]] = conj(c[ig]);
      const double a = pc[2*ig];
      const double b = pc[2*ig+1];
      const int ip = ifftp_[ig];
      const int im = ifftm_[ig];
      pz[2*ip] = a;
      pz[2*ip+1] = b;
      pz[2*im] = a;
      pz[2*im+1] = -b;
    }
  }
  else
    #pragma omp parallel for
    for ( int ig = 0; ig < ng; ig++ )
    {
      // zvec_[ifftp_[ig]] = c[ig];
      const double a = pc[2*ig];
      const double b = pc[2*ig+1];
      const int ip = ifftp_[ig];
      pz[2*ip] = a;
      pz[2*ip+1] = b;
    }
}
////////////////////////////////////////////////////////////////////////////////
void FourierTransform::zvec_to_vector(complex<double> *c)
{
  const int ng = basis_.localsize();
  const double* const pz = (double*) &zvec_[0];
  double* const pc = (double*) &c[0];
  #pragma omp parallel for
  for ( int ig = 0; ig < ng; ig++ )
  {
    // c[ig] = zvec_[ifftp_[ig]];
    const int ip = ifftp_[ig];
    const double pz0 = pz[2*ip];
    const double pz1 = pz[2*ip+1];
    pc[2*ig]   = pz0;
    pc[2*ig+1] = pz1;
  }
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::doublevector_to_zvec(const complex<double> *c1,
  const complex<double> *c2)
{
  // Mapping of two real functions onto zvec
  assert(basis_.real());
  const int zvec_size = zvec_.size();
  double* const pz = (double*) &zvec_[0];
  #pragma omp parallel for
  for ( int i = 0; i < zvec_size; i++ )
  {
    pz[2*i] = 0.0;
    pz[2*i+1] = 0.0;
  }
  const int ng = basis_.localsize();
  const double* const pc1 = (double*) &c1[0];
  const double* const pc2 = (double*) &c2[0];
  #pragma omp parallel for
  for ( int ig = 0; ig < ng; ig++ )
  {
    // const double a = c1[ig].real();
    // const double b = c1[ig].imag();
    // const double c = c2[ig].real();
    // const double d = c2[ig].imag();
    // zvec_[ip] = complex<double>(a-d, b+c);
    // zvec_[im] = complex<double>(a+d, c-b);
    const double a = pc1[2*ig];
    const double b = pc1[2*ig+1];
    const double c = pc2[2*ig];
    const double d = pc2[2*ig+1];
    const int ip = ifftp_[ig];
    const int im = ifftm_[ig];
    pz[2*ip]   = a - d;
    pz[2*ip+1] = b + c;
    pz[2*im]   = a + d;
    pz[2*im+1] = c - b;
  }
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::zvec_to_doublevector(complex<double> *c1,
  complex<double> *c2 )
{
  // Mapping of zvec onto two real functions
  assert(basis_.real());
  const int ng = basis_.localsize();
  const double* const pz = (double*) &zvec_[0];
  double* const pc1 = (double*) &c1[0];
  double* const pc2 = (double*) &c2[0];
  #pragma omp parallel for
  for ( int ig = 0; ig < ng; ig++ )
  {
    // const double a = 0.5*zvec_[ip].real();
    // const double b = 0.5*zvec_[ip].imag();
    // const double c = 0.5*zvec_[im].real();
    // const double d = 0.5*zvec_[im].imag();
    // c1[ig] = complex<double>(a+c, b-d);
    // c2[ig] = complex<double>(b+d, c-a);
    const int ip = ifftp_[ig];
    const int im = ifftm_[ig];
    const double a = pz[2*ip];
    const double b = pz[2*ip+1];
    const double c = pz[2*im];
    const double d = pz[2*im+1];
    pc1[2*ig]   = 0.5 * ( a + c );
    pc1[2*ig+1] = 0.5 * ( b - d );
    pc2[2*ig]   = 0.5 * ( b + d );
    pc2[2*ig+1] = 0.5 * ( c - a );
  }
}


#if defined(FFT_NOLIB)

////////////////////////////////////////////////////////////////////////////////
//
//     /* no library: use cfftm function */
//
//
//     /* Transform along x */
//     int i2;
//     int ntrans = np1*np2;
//     int length = np0;
//     int ainc   = 1;
//     int ajmp   = np0i;
//     int idir;
//     double scale;
//     if ( dir == R_TO_K )
//     {
//       idir = -1;
//       scale = 1.0 / ((double) np0*np1*np2);
//     }
//     else
//     {
//       idir = 1;
//       scale = 1.0;
//     }
//
//     cfftm ( &val[0], &val[0], scale, ntrans, length, ainc, ajmp, idir );
//
//     /* Transform along y */
//     for ( i2 = 0; i2 < np2; i2++ )
//     {
//       int ist = i2 * np0i * np1;
//       ntrans = np0;
//       length = np1;
//       ainc   = np0i;
//       ajmp   = 1;
//       scale = 1.0;
//       cfftm ( &val[ist], &val[ist], scale, ntrans, length, ainc, ajmp, idir );
//     }
//
//     /* Transform along z */
//     ntrans = np0i*np1;
//     length = np2;
//     ainc   = np0i*np1;
//     ajmp   = 1;
//     scale = 1.0;
//     cfftm ( &val[0], &val[0], scale, ntrans, length, ainc, ajmp, idir );
//
////////////////////////////////////////////////////////////////////////////////

/* define multiple FFT function here */

void cfftm ( complex<double> *ain, complex<double> *aout, double scale,
  int ntrans, int length,
  int ainc, int ajmp, int idir )
/*
 *  cfftm: multiple one-dimensional complex FFT
 *
 *  ain     complex array (input)
 *  aout    complex array (output)
 *  scale   global scaling factor
 *  ntrans  number of transforms
 *  length  length of each transform (in complex numbers)
 *  ainc    distance between elements within a transform (in complex numbers)
 *  ajmp    distance between first elements of transforms (in complex numbers)
 *  idir    direction of transform
 */

{
  void cfft ( int idir, complex<double> *z1, complex<double> *z2, int n,
    int *inzee );
  vector<complex<double> > tmpa(length), tmpb(length);
  for ( int it = 0; it < ntrans; it++ )
  {
    int ibase = it * ajmp;
    for ( int i = 0; i < length; i++ )
    {
      tmpa[i] = ain[ibase+i*ainc];
    }
    int inzee = 1;
    cfft ( idir, &tmpa[0], &tmpb[0], length, &inzee );
    if ( inzee == 1 )
      for ( int i = 0; i < length; i++ )
      {
        aout[ibase+i*ainc] = tmpa[i];
      }
    else
      for ( int i = 0; i < length; i++ )
      {
        aout[ibase+i*ainc] = tmpb[i];
      }
    for ( int i = 0; i < length; i++ )
    {
      aout[ibase+i*ainc] *= scale;
    }
  }
}

/*******************************************************************************
 *
 *  Complex FFT
 *  C version
 *
 *  From: C++ Language System Release 3.0 Library Manual
 *  Transcription from 'FFT as Nested Multiplication, with a twist'
 *  C. de Boor, SIAM Sci. Stat. Comput., Vol 1, No 1, March 1980
 *
 *  Adapted to C by F.Gygi, 17 Feb 1993, 9 Dec 1993
 *
 ******************************************************************************/

#include <math.h>

#define NEXTMX 12

void cfft ( int idir, complex<double> *z1, complex<double> *z2, int n,
  int *inzee )
{
  // Compute the discrete Fourier transform of z1 (or z2) in
  // the Cooley-Tukey way, but with a twist.
  // z1[before], z2[before]
  // *inzee == 1 means input in z1; *inzee == 2 means input in z2

  void fftstp ( int idir, complex<double> *zin, int after,
                int now, int before, complex<double> *zout );
  static int prime[NEXTMX] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
  int before = n;
  int after = 1;
  int next = 0;
  int now;


  do
  {
    int np = prime[next];
    if ( (before/np)*np < before )
    {
      if ( ++next < NEXTMX ) continue;
      now = before;
      before = 1;
    }
    else
    {
      now = np;
      before /= np;
    }
    if ( *inzee == 1 )
      fftstp ( idir, z1, after, now, before, z2 );
    else
      fftstp ( idir, z2, after, now, before, z1 );
    *inzee = 3 - *inzee;
    after *= now;
  } while ( before > 1 );

}

void fftstp ( int idir, complex<double> *zin, int after,
              int now, int before, complex<double> *zout )
{

  static const double twopi = 2 * 3.141592653589793;
  double angle;
  complex<double> omega;
  complex<double> arg,value;
  int ia,ib,j,k,l,in;

  angle = twopi/(now*after);
  omega =  complex<double>(cos ( angle ),-idir * sin ( angle ));
  arg = 1.0;
  for ( j = 0; j < now; j++ )
  {
    for ( ia = 0; ia < after; ia++ )
    {
      for ( ib = 0; ib < before; ib++ )
      {
        /* value = zin(ia,ib,now) */
        k = ia + ib*after + (now-1)*before*after;
        value = zin[k];
        for ( in = now-2; in >= 0; in-- )
        {
          /* value = value*arg + zin(ia,ib,in) */
          /* zin(ia,ib,in) = zin[ia + ib*after + in*before*after]; */
          l = ia + ib*after + in*before*after;
          value = value * arg + zin[l];
        }
        /* zout(ia,j,ib) = value */
        /* zout[ia + j*after + ib*now*after] = value; */
        l = ia + j*after + ib*now*after;
        zout[l] = value;
      }

      /* arg *= omega; */
      arg *= omega;

    }
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::reset_timers(void)
{
#if TIMING
  tm_f_map.reset();
  tm_f_fft.reset();
  tm_f_pack.reset();
  tm_f_mpi.reset();
  tm_f_zero.reset();
  tm_f_unpack.reset();

  tm_b_map.reset();
  tm_b_fft.reset();
  tm_b_pack.reset();
  tm_b_mpi.reset();
  tm_b_zero.reset();
  tm_b_unpack.reset();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::forward_1z( vector<complex<double> >& val_in,
                                   vector<complex<double> >& val_out )
{
  fft_1z( val_in, val_out, 1 );
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::backward_1z( vector<complex<double> >& val_in,
                                    vector<complex<double> >& val_out )
{
  fft_1z( val_in, val_out, -1 );
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fft_1z( vector<complex<double> >& val_in,
                               vector<complex<double> >& val_out, int isign )
{
  val_out.resize(np2_);

#if USE_ESSL_FFT
  int inc1 = 1, inc2 = np2_, ntrans = 1, initflag = 0;
  double scale = 1.0;

  if (isign > 0)
     dcft_(&initflag,&val_in[0],&inc1,&inc2,&val_out[0],&inc1,&inc2,&np2_,&ntrans,
           &isign,&scale,&aux1zf1d[0],&naux1z1d,&aux2[0],&naux2);
  else
     dcft_(&initflag,&val_in[0],&inc1,&inc2,&val_out[0],&inc1,&inc2,&np2_,&ntrans,
           &isign,&scale,&aux1zb1d[0],&naux1z1d,&aux2[0],&naux2);

#elif HAVE_FFTW2
   /*
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
  val_out = val_in;
  if( isign > 0 )
  {
    fftw_one(fwplan2,(FFTW_COMPLEX*)&val_out[0],(FFTW_COMPLEX*)0);
  }
  else {
    fftw_one(bwplan2,(FFTW_COMPLEX*)&val_out[0],(FFTW_COMPLEX*)0);
  }
#elif HAVE_FFTW3
  for (int ii=0; ii<np2_; ii++)
     zin_[ii] = val_in[ii];
  if( isign > 0 )
     fftw_execute_dft ( fwplan1z, (fftw_complex*)&zin_[0],
                        (fftw_complex*)&zout_[0]);
  else
     fftw_execute_dft ( bwplan1z, (fftw_complex*)&zin_[0],
                        (fftw_complex*)&zout_[0]);
  for (int ii=0; ii<np2_; ii++)
    val_out[ii] = zout_[ii];
  
#else
  // No library
  /* Transform along z */
  int ntrans = 1;
  int length = np2_;
  int ainc   = 1;
  int ajmp   = np2_;
  double scale = 1.0;
  cfftm ( &val_in[0], &val_out[0], scale, ntrans, length, ainc, ajmp, isign );
#endif

  if( isign > 0 )
  {
    const double fac = 1.0 / np2_;
    for( int i = 0; i < np2_; ++i  )
    {
      val_out[i] *= fac;
    }
  }
}

