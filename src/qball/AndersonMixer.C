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
// AndersonMixer.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "AndersonMixer.h"
#include <math/blas.h>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::restart(void)
{
  n_ = -1;
  k_ = -1;
}

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::update(double* x, double* f, double* xbar, double* fbar)
{
  // update:
  // input: x, f
  // output: xbar, fbar
  //
  // Computes the pair (xbar,fbar) using pairs (x,f) used
  // in previous updates, according to the Anderson algorithm.

  // increment index of current vector
  k_ = ( k_ + 1 ) % ( nmax_ + 1 );
  // increment current number of vectors
  if ( n_ < nmax_ ) n_++;

  // save vectors

  for ( int i = 0; i < m_; i++ )
  {
    x_[k_][i] = x[i];
    f_[k_][i] = f[i];
  }

  valarray<double> a;
  valarray<double> b;
  valarray<double> theta;
  if ( n_ > 0 )
  {
    // compute matrix A = F^T F and rhs b = F^T f
    // compute the lower part of A only (i>=j)
    a.resize(n_*n_);
    b.resize(n_);
    theta.resize(n_);
    for ( int i = 0; i < n_; i++ )
    {
      const int kmi = ( k_ - i + nmax_ ) % ( nmax_ + 1 );
      assert(kmi>=0);
      assert(kmi<nmax_+1);
      // cout << "k=" << k_ << " i=" << i << " kmi=" << kmi << endl;
      for ( int j = 0; j <= i; j++ )
      {
        const int kmj = ( k_ - j + nmax_ ) % ( nmax_ + 1 );
        assert(kmj>=0);
        assert(kmj<nmax_+1);
        // cout << "k=" << k_ << " j=" << j << " kmj=" << kmj << endl;
        double sum = 0.0;
        for ( int l = 0; l < m_; l++ )
          sum += (f_[k_][l] - f_[kmi][l]) * (f_[k_][l] - f_[kmj][l]);
        a[i+j*n_] = sum;
      }
      double bsum = 0.0;
      for ( int l = 0; l < m_; l++ )
        bsum += ( f_[k_][l] - f_[kmi][l] ) * f_[k_][l];
      b[i] = bsum;
    }

    if ( pctxt_ != 0 )
    {
      pctxt_->dsum(n_*n_,1,&a[0],n_*n_);
      pctxt_->dsum(n_,1,&b[0],n_);
    }

#if 0
    // print matrix a and rhs b
    if ( pctxt_ != 0 )
    {
      // print matrix a and rhs b
      double anrm = 0.0;
      for ( int i = 0; i < n_; i++ )
        for ( int j = 0; j <=i; j++ )
        {
          cout << "AMIX4a, mype = " << pctxt_->mype() << ", "
               << "a("<<i<<","<<j<<")=" << a[i+j*n_] << endl;
          anrm += a[i+j*n_]*a[i+j*n_];
        }
      for ( int i = 0; i < n_; i++ )
        cout << "AMIX4b, mype = " << pctxt_->mype() << ", "
             <<  "b("<<i<<")=" << b[i] << endl;
      cout << " AndersonMixer: n=" << n_ << " anorm = " << anrm << endl;
    }
#endif

    // solve the linear system a * theta = b
    // solve on task 0 and bcast result

    //ewd: try doing this on every task, w. no communication
    //if ( pctxt_ == 0 || pctxt_->oncoutpe() )
    {
      const bool diag = false;
      if ( diag )
      {
        // solve the linear system using eigenvalues and eigenvectors
        // compute eigenvalues of a
        char jobz = 'V';
        char uplo = 'L';
        valarray<double> w(n_);
        int lwork = 3*n_;
        valarray<double> work(lwork);
        int info;
        dsyev(&jobz,&uplo,&n_,&a[0],&n_,&w[0],&work[0],&lwork,&info);

        cout << "AndersonMixer: eigenvalues: ";
        for ( int i = 0; i < n_; i++ )
          cout << w[i] << "  ";
        cout << endl;
        if ( info != 0 )
        {
          cerr << " AndersonMixer: Error in dsyev" << endl;
          cerr << " info = " << info << endl;
          exit(1);
        }

        // solve for theta
        // theta_i = sum_j
        for ( int k = 0; k < n_; k++ )
        {
          // correct only if eigenvalue w[k] is large enough compared to the
          // largest eigenvalue
          const double eig_ratio = 1.e-14;
          if ( w[k] > eig_ratio * w[n_-1] )
          {
            const double fac = 1.0 / w[k];
            for ( int i = 0; i < n_; i++ )
              for ( int j = 0; j < n_; j++ )
                theta[i] += fac * a[i+k*n_] * a[j+k*n_] * b[j];
          }
        }
      }
      else
      {
        // solve the linear system directly

        // Tikhonov regularization parameter
        // adjust the parameter until the norm of theta is < 1.0
        double tikhonov_parameter = 1.e-12;
        bool norm_ok = false;
        valarray<double> asave(a);
        valarray<double> bsave(b);
        int iter = 0;
        const int maxiter = 100;
        while ( !norm_ok && iter < maxiter )
        {
          a = asave;
          b = bsave;
          for ( int i = 0; i < n_; i++ )
            a[i+i*n_] += tikhonov_parameter;

          char uplo = 'L';
          int nrhs = 1;
          valarray<int> ipiv(n_);
          valarray<double> work(n_);
          int info;
          dsysv(&uplo,&n_,&nrhs,&a[0],&n_,&ipiv[0],
                &b[0],&n_,&work[0],&n_,&info);
          if ( info != 0 )
          {
            cerr << " AndersonMixer: Error in dsysv" << endl;
            cerr << " info = " << info << endl;
            exit(1);
          }
          // the vector b now contains the solution
          theta = b;

          // check condition on the norm of theta
          norm_ok = true;
#if 0
          // unit simplex criterion
          double theta_sum = 0.0;
          for ( int i = 0; i < theta.size(); i++ )
          {
            theta_sum += theta[i];
            norm_ok &= theta[i] >= 0.0;
          }
          norm_ok &= fabs(theta_sum) <= 1.0;
#endif
#if 0
          // infinity norm criterion
          for ( int i = 0; i < theta.size(); i++ )
            norm_ok &= fabs(theta[i]) <  3.0;
#endif
#if 1
          // 2-norm criterion
          double theta_sum = 0.0;
          for ( int i = 0; i < theta.size(); i++ )
          {
            theta_sum += theta[i] * theta[i];
          }
          norm_ok = theta_sum <= 1.0;
#endif
          // cout << " tp = " << tikhonov_parameter
          //      << " AndersonMixer: theta = ";
          // for ( int i = 0; i < theta.size(); i++ )
          //   cout << theta[i] << " ";
          // cout << endl;

          tikhonov_parameter *= 2.0;
          iter++;
        }
      }

      if ( pctxt_->oncoutpe() ) {
        cout << " AndersonMixer: theta = ";
        for ( int i = 0; i < theta.size(); i++ )
          cout << theta[i] << " ";
        cout << endl;
      }
    }

    //ewd
    /*
    // broadcast theta from task 0
    if ( pctxt_ != 0 )
    {
      if ( pctxt_->oncoutpe() )
        pctxt_->dbcast_send(n_,1,&theta[0],n_);
      else
        pctxt_->dbcast_recv(n_,1,&theta[0],n_,0,0);
    }
    */
#if 0
    for ( int ip = 0; ip < pctxt_->size(); ip++ )
    {
      pctxt_->barrier();
      if ( pctxt_->mype() == ip )
      {
        cout << pctxt_->mype() << ": ";
        cout << " AndersonMixer: theta = ";
        for ( int i = 0; i < theta.size(); i++ )
          cout << theta[i] << " ";
        cout << endl;
      }
    }
#endif

  } // if n_ > 0

  // fbar = f[k] + sum_i theta_i * ( f[k] - f[kmi] )
  // xbar = x[k] + sum_i theta_i * ( x[k] - x[kmi] )
  for ( int l = 0; l < m_; l++ )
  {
    fbar[l] = f_[k_][l];
    xbar[l] = x_[k_][l];
  }
  for ( int i = 0; i < n_; i++ )
  {
    const int kmi = ( k_ - i + nmax_ ) % ( nmax_ + 1 );
    assert(kmi>=0);
    assert(kmi<nmax_+1);
    for ( int l = 0; l < m_; l++ )
    {
      fbar[l] -= theta[i] * ( f_[k_][l] - f_[kmi][l] );
      xbar[l] -= theta[i] * ( x_[k_][l] - x_[kmi][l] );
    }
  }
}
