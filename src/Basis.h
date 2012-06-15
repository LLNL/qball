////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
//  Basis.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Basis.h,v 1.11 2008-09-08 15:56:18 fgygi Exp $

#ifndef BASIS_H
#define BASIS_H

#include "D3vector.h"
#include "UnitCell.h"
#include "Context.h"
#include <vector>

class Context;

class Basis
{
  private:

  Context ctxt_;
  int nprow_, myrow_;

  UnitCell cell_;         // cell dimensions
  UnitCell refcell_;      // reference cell dimensions
  D3vector kpoint_;       // k-point in units of b0,b1,b2
  double ecut_;           // energy cutoff of wavefunctions in Rydberg
  int idxmin_[3];          // minimum index in each direction
  int idxmax_[3];          // maximum index in each direction
  int size_;              // basis size
  int nrods_;             // total number of rods
  std::vector<int> localsize_; // localsize_[ipe]
  int maxlocalsize_, minlocalsize_;
  std::vector<int> nrod_loc_;
  std::vector<std::vector<int> > rod_h_;
  std::vector<std::vector<int> > rod_k_;
  std::vector<std::vector<int> > rod_lmin_;
  std::vector<std::vector<int> > rod_size_;
  std::vector<std::vector<int> > rod_first_;

  std::vector<int>    idx_;   // 3-d index of vectors idx[i*3+j]
  std::vector<double> g_;     // norm of g vectors g[localsize]
  std::vector<double> kpg_;   // norm of g vectors g[localsize]
  std::vector<double> gi_;    // inverse norm of g vectors gi[localsize]
  std::vector<double> kpgi_;  // inverse norm of k+g vectors kpgi[localsize]
  std::vector<double> g2_;    // 2-norm of g vectors g2[localsize]
  std::vector<double> kpg2_;  // 2-norm of g vectors g2[localsize]
  std::vector<double> g2i_;   // inverse square norm of g vec g2i[localsize]
  std::vector<double> kpg2i_; // inverse square norm of k+g vec kpg2i[localsize]
  int np_[3];            // cache for the function np
  std::vector<double> gx_;    // g vec components gx[j*localsize+i], j=0,1,2
  std::vector<double> kpgx_;  // k+g vec components kpgx[j*localsize+i], j=0,1,2
  std::vector<int> isort_loc; // index array to access locally sorted vectors
                         // kpg2_[isort_loc[i]] < kpg2_[isort_loc[j]] if i < j
  bool real_;            // true if k=0
  bool complex_forced_;   
  void update_g(void);

  public:

  const Context& context(void) const; // context on which Basis is defined

  const UnitCell& cell() const;   // cell dimensions
  const UnitCell& refcell() const;// reference cell dimensions
  const D3vector kpoint() const; // k-point in units of b0,b1,b2
  int np(int i) const;           // good size of FFT grid in direction i
  bool factorizable(int n) const;// check if n is factorizable with low factors
  int idxmin(int i) const;       // smallest index in direction i
  int idxmax(int i) const;       // largest index in direction i
  double ecut() const;           // energy cutoff in Hartree
  bool real() const;             // return true if kpoint == (0,0,0)
  void force_complex();          // force complex basis even if kpoint == (0,0,0)
  bool complex_forced() const;   // return complex_forced_

  int size() const;              // total number of g vectors
  int localsize() const;         // local number of g vectors on current process
  int localsize(int ipe) const;  // local number of g vectors on process ipe
  int maxlocalsize() const;      // largest local size
  int minlocalsize() const;      // smallest local size

  int nrods() const;             // total number of rods
  int nrod_loc() const;          // local number of rods on current process
  int nrod_loc(int ipe) const;   // local number of rods on process ipe

  int rod_h(int irod) const;     // h-position of rod irod on current process
  int rod_h(int ipe, int irod) const; // h-position of rod irod on process ipe

  int rod_k(int irod) const;     // k-position of rod irod on current process
  int rod_k(int ipe, int irod) const; // k-position of rod irod on process ipe

  int rod_lmin(int irod) const;  // lmin-position of rod irod on current process
  int rod_lmin(int ipe, int irod) const; // lmin-pos. of rod irod on process ipe

  // size of rod irod
  int rod_size(int irod) const;
  int rod_size(int ipe, int irod) const;

  // local position of first elem. of rod irod
  int rod_first(int irod) const;
  int rod_first(int ipe, int irod) const;

  int    idx(int i) const;   // integer indices of vectors idx[i*3+j]
  double g(int i) const;     // norm of g vectors g[i]
  double kpg(int i) const;   // norm of k+g vectors kpg[i]
  double gi(int i) const;    // inverse norm of g vectors gi[i]
  double kpgi(int i) const;  // inverse norm of k+g vectors kpgi[i]
  double g2(int i) const;    // 2-norm of g vectors g2[i]
  double kpg2(int i) const;  // 2-norm of k+g vectors kpg[i]
  double g2i(int i) const;   // inverse square norm of g g2i[i]
  double kpg2i(int i) const; // inverse square norm of k+g kpg2i[i]
  double gx(int i) const;    // g vectors gx[i+localsize*j],j=0,1,2
  double kpgx(int i) const;  // k+g vectors kpgx[i+localsize*j],j=0,1,2

  int isort(int i) const;    // index of vectors locally sorted by norm

  const int*    idx_ptr(void) const;
  const double* g_ptr(void) const;
  const double* kpg_ptr(void) const;
  const double* gi_ptr(void) const;
  const double* kpgi_ptr(void) const;
  const double* g2_ptr(void) const;
  const double* kpg2_ptr(void) const;
  const double* g2i_ptr(void) const;
  const double* kpg2i_ptr(void) const;
  const double* gx_ptr(int j) const;
  const double* kpgx_ptr(int j) const;

  double memsize(void) const;
  double localmemsize(void) const;

  Basis(const Context &ctxt, D3vector kpoint);
  Basis(const Context &ctxt, D3vector kpoint, bool force_complex);
  //Basis(const Basis &b);
  ~Basis(void);

  bool resize(const UnitCell& cell, const UnitCell& refcell, double ecut);
  void print(std::ostream& os);
  void print_casino(ostream& os);

};
std::ostream& operator << ( std::ostream& os, Basis& b );
#endif
