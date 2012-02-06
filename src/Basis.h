////////////////////////////////////////////////////////////////////////////////
//
//  Basis.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Basis.h,v 1.4 2008/09/25 18:45:24 draeger1 Exp $

#ifndef BASIS_H
#define BASIS_H

#include "D3vector.h"
#include "UnitCell.h"

class Context;

class Basis
{
  private:
  
  struct BasisImpl* pimpl_;
  
  public:
  
  const Context& context(void) const; // context on which Basis is defined

  const UnitCell& cell() const;   // cell dimensions
  const UnitCell& refcell() const;// reference cell dimensions
  const D3vector kpoint() const; // k-point in units of b0,b1,b2
  int np(int i) const;           // good size of FFT grid in direction i
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
  double kpgi(int i) const;    // inverse norm of k+g vectors kpgi[i]
  double g2(int i) const;    // 2-norm of g vectors g2[i]
  double kpg2(int i) const;  // 2-norm of k+g vectors kpg[i]
  double g2i(int i) const;   // inverse square norm of g g2i[i]
  double kpg2i(int i) const;   // inverse square norm of k+g kpg2i[i]
  double gx(int i) const;    // g vectors gx[i+localsize*j],j=0,1,2
  double kpgx(int i) const;    // k+g vectors kpgx[i+localsize*j],j=0,1,2
  double gx2(int i) const;   // g vectors components^2 gx2[i+localsize*j]
  
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
  const double* gx2_ptr(int j) const;
  
  double memsize(void) const;
  double localmemsize(void) const;

  Basis(const Context &ctxt, D3vector kpoint);
  Basis(const Context &ctxt, D3vector kpoint, bool force_complex);
  Basis(const Basis &b);
  ~Basis(void);
  
  bool resize(const UnitCell& cell, const UnitCell& refcell, double ecut);
  void print(ostream& os);
  void print_casino(ostream& os);
};
ostream& operator << ( ostream& os, Basis& b );
#endif
