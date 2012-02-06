////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.h,v 1.3 2010/03/01 18:05:45 draeger1 Exp $

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <vector>
#include <valarray>
#include <cassert>
#include "Context.h"

class AndersonMixer
{
  // nmax is the dimension of the subspace of previous search directions
  // nmax=0: use simple mixing (no acceleration)
  // nmax=1: use one previous direction
  int     m_;                    // dimension of vectors
  int     nmax_;                 // maximum number of vectors (without current)
  int     n_;                    // number of vectors
  int     k_;                    // index of current vector
  const   Context* const pctxt_; // pointer to relevant Context, null if local

  std::vector<std::valarray<double> > x_,f_;

  public:

  AndersonMixer(const int m, const int nmax, const Context* const pctxt) :
    m_(m), nmax_(nmax), pctxt_(pctxt)
  {
    assert( nmax >= 0 );
    x_.resize(nmax_+1);
    f_.resize(nmax_+1);
    for ( int n = 0; n < nmax_+1; n++ )
    {
      x_[n].resize(m_);
      f_[n].resize(m_);
    }
    restart();
  }

  void update(double* x, double* f, double* xbar, double* fbar);
  void restart(void);
};
#endif
