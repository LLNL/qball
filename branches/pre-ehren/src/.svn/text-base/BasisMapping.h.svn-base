////////////////////////////////////////////////////////////////////////////////
//
//  BasisMapping.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BasisMapping.h,v 1.1 2008/05/13 18:28:54 draeger1 Exp $

#ifndef BASISMAPPING_H
#define BASISMAPPING_H

#include <complex>
#include <vector>

class Basis;
class Context;

class BasisMapping
{
  private:

  const Context& ctxt_;
  const Basis& basis_;
  int nprocs_, myproc_;

  int np0_, np1_, np2_, np012_, np012loc_;
  int nvec_;

  std::vector<int> np2_loc_;   // np2_loc_[iproc], iproc=0, nprocs_-1
  std::vector<int> np2_first_; // np2_first_[iproc], iproc=0, nprocs_-1

  std::vector<int> scounts, sdispl, rcounts, rdispl;
  std::vector<std::complex<double> > sbuf, rbuf;

  std::vector<int> ip_, im_;
  std::vector<int> ipack_, iunpack_;

  public:

  BasisMapping (const Basis &basis);
  int np0(void) const { return np0_; }
  int np1(void) const { return np1_; }
  int np2(void) const { return np2_; }
  int np2loc(void) const { return np2_loc_[myproc_]; }
  int np012(void) const { return np012_; }
  int np012loc(void) const { return np012loc_; }
  int nvec(void) const { return nvec_; }
  int zvec_size(void) const { return nvec_ * np2_; }

  const Context& context(void) const { return ctxt_; }

  // map a function c(G) to zvec_
  void vector_to_zvec(const std::complex<double> *c,
                      std::complex<double> *zvec);
  // map zvec_ to a function c(G)
  void zvec_to_vector(const std::complex<double> *zvec,
                      std::complex<double> *c);

  void transpose_fwd(const std::complex<double> *zvec,
                     std::complex<double> *ct);
  void transpose_bwd(const std::complex<double> *ct,
                     std::complex<double> *zvec);
};
#endif
