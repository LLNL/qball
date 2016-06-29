/*******************************************************************************
 *
 * spline.c
 *
 ******************************************************************************/

#include "spline.h"
#include <assert.h>
//ewd DEBUG
#include <iostream>

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{

  int i,k;
  double p,qn,sig,un,*u = new double[n];

  if ( yp1 >= 1.e30 )
  {
    y2[0] = 0.0;
    u[0] = 0.0;
  }
  else
  {
    y2[0] = -0.5;
    assert ( x[1] - x[0] > 0.0 );
    u[0] = ( 3.0 / (x[1]-x[0]) ) * ( (y[1]-y[0]) / (x[1]-x[0]) - yp1 );
  }

  for ( i = 1; i < n-1; i++ )
  {
    assert ( x[i+1] > x[i] );
    sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
    p = sig * y2[i-1] + 2.0;
    y2[i] = ( sig - 1.0 ) / p;
    u[i] = ( 6.0 * ( ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) -
                     ( y[i] - y[i-1] ) / ( x[i] - x[i-1] ) ) /
             ( x[i+1] - x[i-1] ) - sig * u[i-1] ) / p;
  }

  if ( ypn >= 1.e30 )
  {
    qn = 0.0;
    un = 0.0;
  }
  else
  {
    qn = 0.5;
    un = ( 3.0 / (x[n-1]-x[n-2]) ) * 
         ( ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]) );
  }

  y2[n-1] = ( un - qn * u[n-2] ) / ( qn * y2[n-2] + 1.0 );

  for ( k = n-2; k >= 0; k-- )
  {
    y2[k] = y2[k] * y2[k+1] + u[k];
  }

  delete [] u;
}

void splint (double *xa, double *ya, double *y2a, int n, double x, double *y)
{
  int k,khi,klo;
  double a,b,h;

  klo = 0;
  khi = n-1;

  while ( khi - klo > 1 )
  {
    k = ( khi + klo ) / 2;
    if ( xa[k] > x )
      khi = k;
    else
      klo = k;
  }

  //ewd DEBUG
  if (khi > n-1) {
    std::cout << "ERROR.SPLINT:  khi = " << khi << ", n = " << n << std::endl;
    return;
  }
  if (klo > n-1) {
    std::cout << "ERROR.SPLINT:  klo = " << klo << ", n = " << n << std::endl;
    return;
  }
  
  h = xa[khi] - xa[klo];
  assert ( h > 0.0 );

  a = ( xa[khi] - x ) / h;
  b = ( x - xa[klo] ) / h;

  *y = a * ya[klo] + b * ya[khi] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi] );

}

void splintd (double *xa, double *ya, double *y2a,
 int n, double x, double *y, double *dy)
{
  int k,khi,klo;
  double a,b,h;

  klo = 0;
  khi = n-1;

  while ( khi - klo > 1 )
  {
    k = ( khi + klo ) / 2;
    if ( xa[k] > x )
      khi = k;
    else
      klo = k;
  }

  h = xa[khi] - xa[klo];
  assert ( h > 0.0 );

  a = ( xa[khi] - x ) / h;
  b = ( x - xa[klo] ) / h;

  *y = a * ya[klo] + b * ya[khi] + h * h * (1.0/6.0) *
       ( (a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi] );

  *dy = ( ya[khi] - ya[klo] ) / h +
        h * ( ( (1.0/6.0) - 0.5 * a * a ) * y2a[klo] +
              ( 0.5 * b * b - (1.0/6.0) ) * y2a[khi] );
}
