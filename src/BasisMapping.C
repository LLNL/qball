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
//  BasisMapping.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Basis.h"
#include "Context.h"
#include "BasisMapping.h"

#include <iostream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BasisMapping::BasisMapping (const Basis &basis) : ctxt_(basis.context()),
 basis_(basis)

{
  assert(ctxt_.npcol() == 1);
  nprocs_ = ctxt_.size();
  myproc_ = ctxt_.myproc();

  np0_ = basis.np(0);
  np1_ = basis.np(1);
  np2_ = basis.np(2);
  np012_ = np0_ * np1_ * np2_;

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
  np012loc_ = np0_ * np1_ * np2_loc_[myproc_];

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
    // compute index arrays ip_ and im_ for mapping vector->zvec
    ip_.resize(basis_.localsize());
    im_.resize(basis_.localsize());

    if ( myproc_ == 0 )
    {
      // this process holds rod(0,0)
      // nvec_ == 2 * nrod_loc - 1

      // map rod(0,0)
      // the positive segment of rod(0,0) maps onto the first half of
      // the first column of zvec_, and the negative segment maps onto
      // the second half
      int ig = 0;
      ip_[0] = 0;
      im_[0] = 0;
      ig++;
      for ( int l = 1; l < basis_.rod_size(0); l++ )
      {
        ip_[ig] = l;
        im_[ig] = np2_ - l;
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
          ip_[ig] = ( 2 * irod - 1 ) * np2_ + izp;
          im_[ig] = ( 2 * irod ) * np2_ + izm;
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
      // rod(-h,-k) maps onto column (2*irod+1)*np2_ of zvec_,irod=0,..,nrods-1
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
          ip_[ig] = ( 2 * irod ) * np2_ + izp;
          im_[ig] = ( 2 * irod + 1 ) * np2_ + izm;
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
    // compute index array ip_ for mapping vector->zvec
    // Note: im_ is not used
    ip_.resize(basis_.localsize());

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
        ip_[ig] = irod * np2_ + iz;
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
void BasisMapping::transpose_fwd(const complex<double> *zvec,
  complex<double> *ct)
{
  // Transpose zvec to ct
  // scatter zvec to sbuf for transpose
#if USE_GATHER_SCATTER
  // zsctr: y(indx(i)) = x(i)
  // void zsctr_(int* n, complex<double>* x, int* indx, complex<double>* y);
  {
    complex<double>* y = &sbuf[0];
    complex<double>* x = const_cast<complex<double>*>(zvec);
    int n = zvec_.size();
    zsctr_(&n,x,&ipack_[0],y);
  }
#else
  const int len = zvec_size();
  double* const ps = (double*) &sbuf[0];
  const double* const pz = (const double*) zvec;
  for ( int i = 0; i < len; i++ )
  {
    // sbuf[ipack_[i]] = zvec[i];
    const int ip = ipack_[i];
    const double a = pz[2*i];
    const double b = pz[2*i+1];
    ps[2*ip]   = a;
    ps[2*ip+1] = b;
  }
#endif

  // segments of z-vectors are now in sbuf

  // transpose
#if USE_MPI
  int status = MPI_Alltoallv((double*)&sbuf[0],&scounts[0],&sdispl[0],
      MPI_DOUBLE,(double*)&rbuf[0],&rcounts[0],&rdispl[0],MPI_DOUBLE,
      ctxt_.comm());
  if ( status != 0 )
  {
    cout << " BasisMapping: status = " << status << endl;
    ctxt_.abort(2);
  }
#else
  assert(sbuf.size()==rbuf.size());
  rbuf = sbuf;
#endif

  // copy from rbuf to ct
  // scatter index array iunpack
  {
    const int len = np012loc_;
    double* const pv = (double*) ct;
    for ( int i = 0; i < len; i++ )
    {
      pv[2*i]   = 0.0;
      pv[2*i+1] = 0.0;
    }
  }

#if USE_GATHER_SCATTER
  // zsctr(n,x,indx,y): y(indx(i)) = x(i)
  {
    complex<double>* y = ct;
    complex<double>* x = &rbuf[0];
    int n = rbuf.size();
    zsctr_(&n,x,&iunpack_[0],y);
  }
#else
  {
    const int rbuf_size = rbuf.size();
    const double* const pr = (double*) &rbuf[0];
    double* const pv = (double*) ct;
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

  // coefficients are now in ct
}

////////////////////////////////////////////////////////////////////////////////
void BasisMapping::transpose_bwd(const complex<double> *ct,
  complex<double> *zvec)
{
  // transpose back distributed array ct into zvec
  // gather ct into rbuf
#if USE_GATHER_SCATTER
  // zgthr: x(i) = y(indx(i))
  // void zgthr_(int* n, complex<double>* y, complex<double>* x, int*indx);
  {
    complex<double>* y = const_cast<complex<double>*>(ct);
    complex<double>* x = &rbuf[0];
    int n = rbuf.size();
    zgthr_(&n,y,x,&iunpack_[0]);
  }
#else
  const int rbuf_size = rbuf.size();
  double* const pr = (double*) &rbuf[0];
  const double* const pv = (const double*) ct;
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
#if USE_MPI
  int status = MPI_Alltoallv((double*)&rbuf[0],&rcounts[0],&rdispl[0],
      MPI_DOUBLE,(double*)&sbuf[0],&scounts[0],&sdispl[0],MPI_DOUBLE,
      ctxt_.comm());
  assert ( status == 0 );
#else
  assert(sbuf.size()==rbuf.size());
  rbuf = sbuf;
#endif

  // segments of z-vectors are now in sbuf
  // gather sbuf into zvec_

#if USE_GATHER_SCATTER
  // zgthr: x(i) = y(indx(i))
  // void zgthr_(int* n, complex<double>* y, complex<double>* x, int*indx);
  {
    complex<double>* y = &sbuf[0];
    complex<double>* x = zvec;
    int n = zvec_.size();
    zgthr_(&n,y,x,&ipack_[0]);
  }
#else
  const int len = zvec_size();
  const double* const ps = (double*) &sbuf[0];
  double* const pz = (double*) zvec;
  for ( int i = 0; i < len; i++ )
  {
    // zvec[i] = sbuf[ipack_[i]];
    const int ip = ipack_[i];
    const double a = ps[2*ip];
    const double b = ps[2*ip+1];
    pz[2*i]   = a;
    pz[2*i+1] = b;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void BasisMapping::vector_to_zvec(const complex<double> *c,
  complex<double> *zvec)
{
  // map coefficients from the basis order to a zvec
  const int ng = basis_.localsize();
  const int len = zvec_size();
  double* const pz = (double*) zvec;
  for ( int i = 0; i < len; i++ )
  {
    pz[2*i]   = 0.0;
    pz[2*i+1] = 0.0;
  }
  const double* const pc = (const double*) c;
  if ( basis_.real() )
  {
    for ( int ig = 0; ig < ng; ig++ )
    {
      // zvec[ip_[ig]] = c[ig];
      // zvec[im_[ig]] = conj(c[ig]);
      const double a = pc[2*ig];
      const double b = pc[2*ig+1];
      const int ip = ip_[ig];
      const int im = im_[ig];
      pz[2*ip] = a;
      pz[2*ip+1] = b;
      pz[2*im] = a;
      pz[2*im+1] = -b;
    }
  }
  else
    for ( int ig = 0; ig < ng; ig++ )
    {
      // zvec[ip_[ig]] = c[ig];
      const double a = pc[2*ig];
      const double b = pc[2*ig+1];
      const int ip = ip_[ig];
      pz[2*ip] = a;
      pz[2*ip+1] = b;
    }
}
////////////////////////////////////////////////////////////////////////////////
void BasisMapping::zvec_to_vector(const complex<double> *zvec,
  complex<double> *c)
{
  const int ng = basis_.localsize();
  const double* const pz = (const double*) zvec;
  double* const pc = (double*) c;
  for ( int ig = 0; ig < ng; ig++ )
  {
    // c[ig] = zvec[ip_[ig]];
    const int ip = ip_[ig];
    const double pz0 = pz[2*ip];
    const double pz1 = pz[2*ip+1];
    pc[2*ig]   = pz0;
    pc[2*ig+1] = pz1;
  }
}
