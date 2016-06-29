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
///////////////////////////////////////////////////////////////////////////////
//
// spline.C
//
///////////////////////////////////////////////////////////////////////////////
// $Id: spline.C,v 1.5 2008/09/08 15:56:20 fgygi Exp $
#include "spline.h"
#include <cassert>

void tridsolve(int n, double* d, double* e, double* f, double* x)
{
  // solve the tridiagonal system Ax=b
  // d[i] = a(i,i)
  // e[i] = a(i,i+1) (superdiagonal of A, e[n-1] not defined)
  // f[i] = a(i,i-1) (subdiagonal of A, f[0] not defined)
  // x[i] = right-hand side b as input
  // x[i] = solution as output

  for ( int i = 1; i < n; i++ )
  {
    f[i] /= d[i-1];
    d[i] -= f[i]*e[i-1];
  }

  for ( int i = 1; i < n; i++ )
    x[i] -= f[i]*x[i-1];

  x[n-1] /= d[n-1];

  for ( int i = n-2; i >= 0; i-- )
    x[i] = (x[i]-e[i]*x[i+1])/d[i];
}

void spline(int n, double *x, double *y, double yp_left, double yp_right,
  int bcnat_left, int bcnat_right, double *y2)
{
  const double third = 1.0/3.0;
  const double sixth = 1.0/6.0;
  double *d = new double[n];
  double *e = new double[n];
  double *f = new double[n];
  if ( bcnat_left == 0 )
  {
    // use derivative yp_left at x[0]
    const double h = x[1]-x[0];
    assert(h>0.0);
    d[0] = third*h;
    e[0] = sixth*h;
    f[0] = 0.0;
    y2[0] = (y[1]-y[0])/h - yp_left;
  }
  else
  {
    // use natural spline at x[0]
    d[0] = 1.0;
    e[0] = 0.0;
    f[0] = 0.0;
    y2[0] = 0.0;
  }
  if ( bcnat_right == 0 )
  {
    // use derivative yp_right at x[n-1]
    const double h = x[n-1]-x[n-2];
    assert(h>0.0);
    d[n-1] = third*h;
    e[n-1] = 0.0;
    f[n-1] = sixth*h;
    y2[n-1] = yp_right - (y[n-1]-y[n-2])/h;
  }
  else
  {
    // use natural spline at x[n-1]
    d[n-1] = 1.0;
    e[n-1] = 0.0;
    f[n-1] = 0.0;
    y2[n-1] = 0.0;
  }

  // tridiagonal matrix
  for ( int i = 1; i < n-1; i++ )
  {
    const double hp = x[i+1]-x[i];
    const double hm = x[i]-x[i-1];
    assert(hp>0.0);
    assert(hm>0.0);
    d[i] = third * (hp+hm);
    e[i] = sixth * hp;
    f[i] = sixth * hm;
    y2[i] = (y[i+1]-y[i])/hp - (y[i]-y[i-1])/hm;
  }

  tridsolve(n,d,e,f,y2);

  delete [] d;
  delete [] e;
  delete [] f;
}

void splint (int n, double *xa, double *ya, double *y2a, double x, double *y)
{
  int k;
  double a,b,h;

  int kl = 0;
  int kh = n-1;

  while ( kh - kl > 1 )
  {
    k = ( kh + kl ) / 2;
    if ( xa[k] > x )
      kh = k;
    else
      kl = k;
  }

  h = xa[kh] - xa[kl];
  assert ( h > 0.0 );

  a = ( xa[kh] - x ) / h;
  b = ( x - xa[kl] ) / h;

  *y = a * ya[kl] + b * ya[kh] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[kl] + (b*b*b-b) * y2a[kh] );

}

void splintd (int n, double *xa, double *ya, double *y2a,
              double x, double *y, double *dy)
{
  int k;
  double a,b,h;

  int kl = 0;
  int kh = n-1;

  while ( kh - kl > 1 )
  {
    k = ( kh + kl ) / 2;
    if ( xa[k] > x )
      kh = k;
    else
      kl = k;
  }

  h = xa[kh] - xa[kl];
  assert ( h > 0.0 );

  a = ( xa[kh] - x ) / h;
  b = ( x - xa[kl] ) / h;

  *y = a * ya[kl] + b * ya[kh] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[kl] + (b*b*b-b) * y2a[kh] );

  *dy = ( ya[kh] - ya[kl] ) / h +
        h * ( ( (1.0/6.0) - 0.5 * a * a ) * y2a[kl] +
              ( 0.5 * b * b - (1.0/6.0) ) * y2a[kh] );
}
