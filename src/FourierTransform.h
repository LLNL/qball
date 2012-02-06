////////////////////////////////////////////////////////////////////////////////
//
// FourierTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FourierTransform.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include <complex>
#include <vector>
using namespace std;

#if USE_FFTW
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
  
  vector<int> np2_loc_; // np2_loc_[iproc], iproc=0, nprocs_-1
  vector<int> np2_first_; // np2_first_[iproc], iproc=0, nprocs_-1
  vector<complex<double> > zvec_;
  
  vector<int> scounts, sdispl, rcounts, rdispl;
  vector<complex<double> > sbuf, rbuf;

  vector<int> ifftp_, ifftm_;
  vector<int> ipack_, iunpack_;
  
  void init_lib(void);
  
#if USE_ESSL
#if USE_ESSL_2DFFT
  vector<double> aux1xyf;
  vector<double> aux1xyb;
  int naux1xy;
#else
  vector<double> aux1xf, aux1yf, aux1zf;
  vector<double> aux1xb, aux1yb, aux1zb;
  vector<double> aux2;
  int naux1x,naux1y,naux1z,naux2;
#endif
#elif USE_FFTW || USE_FFTW3
  fftw_plan fwplan0,fwplan1,fwplan2,bwplan0,bwplan1,bwplan2;
#else
  // no library
#endif

  void vector_to_zvec(const complex<double>* c);
  void zvec_to_vector(complex<double>* c);
  void doublevector_to_zvec(const complex<double>* c1,
       const complex<double> *c2);
  void zvec_to_doublevector(complex<double>* c1, complex<double>* c2);
  void fwd(complex<double>* val);
  void bwd(complex<double>* val);
       
  public:

  FourierTransform (const Basis &basis, int np0, int np1, int np2);
  ~FourierTransform ();
  const Context& context(void) const { return ctxt_; }
  
  // backward: Fourier synthesis, compute real-space function
  // forward:  Fourier analysis, compute Fourier coefficients
  // forward transform includes scaling by 1/np012
  // single transforms: c -> f, f -> c
  void backward (const complex<double>* c, complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(complex<double>* f, complex<double>* c);
  
  // double transforms: c1 + i*c2 -> f, f -> c1 + i*c2
  void backward (const complex<double>* c1, const complex<double>* c2, 
                 complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(complex<double>* f, 
               complex<double>* c1, complex<double>* c2);
  
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
