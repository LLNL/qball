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
// FourierTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FourierTransform.h,v 1.13 2008-09-08 15:56:18 fgygi Exp $

#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include <complex>
#include <vector>

#if USE_FFTW || USE_SPIRAL
#include "fftw.h"
#endif

#include "Timer.h"

class Basis;
class Context;

class FourierTransform
{
  private:

  const Context& ctxt_;
  const Basis& basis_;
  int nprocs_, myproc_;

  int np0_,np1_,np2_;
  int ntrans0_,ntrans1_,ntrans2_;

  int nvec_;

  std::vector<int> np2_loc_; // np2_loc_[iproc], iproc=0, nprocs_-1
  std::vector<int> np2_first_; // np2_first_[iproc], iproc=0, nprocs_-1
  std::vector<std::complex<double> > zvec_;
#ifdef USE_SPIRAL
  std::vector<std::complex<double> > zout_;
  std::vector<std::complex<double> > vin_;
  std::vector<std::complex<double> > vout_;
#endif
  std::vector<int> scounts, sdispl, rcounts, rdispl;
  std::vector<std::complex<double> > sbuf, rbuf;

  std::vector<int> ifftp_, ifftm_;
  std::vector<int> ipack_, iunpack_;

  void init_lib(void);

#if USE_ESSL
#if USE_ESSL_2DFFT
  std::vector<double> aux1xyf;
  std::vector<double> aux1xyb;
  int naux1xy;
#else
  std::vector<double> aux1xf, aux1yf, aux1zf;
  std::vector<double> aux1xb, aux1yb, aux1zb;
  std::vector<double> aux2;
  int naux1x,naux1y,naux1z,naux2;
#endif
#elif USE_FFTW || USE_FFTW3 || USE_SPIRAL
  fftw_plan fwplan0,fwplan1,fwplan2,bwplan0,bwplan1,bwplan2;
#else
  // no library
#endif

  void vector_to_zvec(const std::complex<double>* c);
  void zvec_to_vector(std::complex<double>* c);
  void doublevector_to_zvec(const std::complex<double>* c1,
       const std::complex<double> *c2);
  void zvec_to_doublevector(std::complex<double>* c1, std::complex<double>* c2);
  void fwd(std::complex<double>* val);
  void bwd(std::complex<double>* val);

  public:

  FourierTransform (const Basis &basis, int np0, int np1, int np2);
  ~FourierTransform ();
  const Context& context(void) const { return ctxt_; }

  // backward: Fourier synthesis, compute real-space function
  // forward:  Fourier analysis, compute Fourier coefficients
  // forward transform includes scaling by 1/np012
  // single transforms: c -> f, f -> c
  void backward (const std::complex<double>* c, std::complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f, std::complex<double>* c);

  // double transforms: c1 + i*c2 -> f, f -> c1 + i*c2
  void backward (const std::complex<double>* c1,
                 const std::complex<double>* c2, std::complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f,
               std::complex<double>* c1, std::complex<double>* c2);

  int np0() const { return np0_; }
  int np1() const { return np1_; }
  int np2() const { return np2_; }
  int np2_loc() const { return np2_loc_[myproc_]; }
  int np2_loc(int iproc) const { return np2_loc_[iproc]; }
  int np2_first() const { return np2_first_[myproc_]; }
  int np2_first(int iproc) const { return np2_first_[iproc]; }
  int np012() const { return np0_ * np1_ * np2_; }
  int np012loc(int iproc) const { return np0_ * np1_ * np2_loc_[iproc]; }
  int np012loc() const { return np0_ * np1_ * np2_loc_[myproc_]; }
  int index(int i, int j, int k) const
  { return i + np0_ * ( j +  np1_ * k ); }

  void reset_timers(void);
  Timer tm_f_map, tm_f_fft, tm_f_pack, tm_f_mpi, tm_f_zero, tm_f_unpack,
        tm_b_map, tm_b_fft, tm_b_pack, tm_b_mpi, tm_b_zero, tm_b_unpack;
};
#endif
