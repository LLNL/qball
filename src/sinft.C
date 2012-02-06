/*******************************************************************************
 *
 * sinft.c
 *
 ******************************************************************************/

#include "sinft.h"
#include <math.h>
#include <assert.h>

void sinft ( double *y, int n )
{
  /*
   * Calculates the sine transform of a set of n-real-valued data points
   * stored in array y[0:n-1]. The number n must be a power of 2. On exit
   * y is replaced by its transform. This program, without change, also
   * calculates the inverse sine transform, but in this case the output
   * array should be multiplied by 2/n.
   * Numerical Recipes, 2nd Edition.
   */

  int j;
  double sum,y1,y2,wi,wpi,wr,wpr,wtemp;
  double theta = 3.141592653589793 /((double) n);

  wr = 1.0;
  wi = 0.0;
  wtemp = sin ( 0.5 * theta );
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin ( theta );

  y[0] = 0.0;
  for ( j = 0; j < n/2; j++ )
  {
    wtemp = wr;
    wr = wr * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
    y1 = wi * ( y[j+1] + y[n-j-1] );
    y2 = 0.5 * ( y[j+1] - y[n-j-1] );
    y[j+1] = y1 + y2;
    y[n-j-1] = y1 - y2;
  }

  realft(y,n,1);

  sum = 0.0;
  y[0] *= 0.5;
  y[1] = 0.0;
  for ( j = 0; j < n-1; j++,j++ )
  {
    sum += y[j];
    y[j] = y[j+1];
    y[j+1] = sum;
  }
}

void cosft1 ( double *y, int n )
{
  /* Note: the array  y contains n+1 elements */

  int j;
  double sum,y1,y2,wi,wpi,wr,wpr,wtemp;
  double theta = 3.141592653589793 / ( (double) n );

  wr = 1.0;
  wi = 0.0;
  wtemp = sin ( 0.5 * theta );
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin ( theta );

  sum = 0.5 * ( y[0] - y[n] );
  y[0] = 0.5 * ( y[0] + y[n] );
  for ( j = 0; j < n/2-1; j++ )
  {
    wtemp = wr;
    wr = wr * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
    y1 = 0.5 * ( y[j+1] + y[n-j-1] );
    y2 = y[j+1] - y[n-j-1];
    y[j+1] = y1 - wi * y2;
    y[n-j-1] = y1 + wi * y2;
    sum = sum + wr * y2;
  }

  realft(y,n,1);

  y[n] = y[1];
  y[1] = sum;
  for ( j = 3; j < n; j++,j++ )
  {
    sum += y[j];
    y[j] = sum;
  }
}

void realft ( double *data, int n, int isign )
{
  int i,i1,i2,i3,i4,np3;
  double c1,c2,h1i,h1r,h2i,h2r,wis,wrs,wr,wi,wpr,wpi,wtemp;
  double theta = 3.141592653589793 / ( (double) n/2 );

  c1 = 0.5;
  if ( isign == 1 )
  {
    c2 = -0.5;
    four1(data,n/2,1);
  }
  else
  {
    c2 = 0.5;
    theta = -theta;
  }
	  
  wtemp = sin ( 0.5 * theta );
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin ( theta );

  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n+3;

  /*
   Ftn pgm:
   i1 = 3 5 7 9
   i2 = 4 6 8 10
   i3 = n-1 n-3 n-5 ...
   i4 = n   n-2 n-4 ...

   C pgm:
   i1 = 2 4 6 8
   i2 = 3 5 7 9
   i3 = n-2 n-4 n-6
   i4 = n-1 n-3 n-5

  */

  for ( i = 2; i <= n/4; i++ )
  {
    i1 = i+i-1;
    i2 = i1+1;
    i3 = np3-i2;
    i4 = i3+1;
    assert ( i1 >= 1 ); assert ( i1 <= n );
    assert ( i2 >= 1 ); assert ( i2 <= n );
    assert ( i3 >= 1 ); assert ( i3 <= n );
    assert ( i4 >= 1 ); assert ( i4 <= n );
    wrs = wr;
    wis = wi;
    h1r =  c1 * ( data[i1-1] + data[i3-1] );
    h1i =  c1 * ( data[i2-1] - data[i4-1] );
    h2r = -c2 * ( data[i2-1] + data[i4-1] );
    h2i =  c2 * ( data[i1-1] - data[i3-1] );
    data[i1-1] =  h1r + wrs * h2r - wis * h2i;
    data[i2-1] =  h1i + wrs * h2i + wis * h2r;
    data[i3-1] =  h1r - wrs * h2r + wis * h2i;
    data[i4-1] = -h1i + wrs * h2i + wis * h2r;
    wtemp = wr;
    wr = wr * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
  }

  if ( isign == 1 )
  {
    h1r = data[0];
    data[0] = h1r + data[1];
    data[1] = h1r - data[1];
  }
  else
  {
    h1r = data[0];
    data[0] = c1 * ( h1r + data[1] );
    data[1] = c1 * ( h1r - data[1] );
    four1(data,n/2,-1);
  }
}

void four1 ( double *data, int nn, int isign )
{
  int i,istep,j,m,mmax,n;
  double tempr,tempi,wr,wi,wpr,wpi,wtemp,theta;
      
  n = 2 * nn;
  j = 1;

  for ( i = 1; i < n; i++,i++ )
  {
    if ( j > i )
    {
      tempr     = data[j-1];
      tempi     = data[j];
      data[j-1] = data[i-1];
      data[j]   = data[i];
      data[i-1] = tempr;
      data[i]   = tempi;
    }
    m = n/2;
    while ( ( m >= 2 ) && ( j > m ) )
    {
      j -= m;
      m /= 2;
    }
    j += m;
  }
  mmax = 2;

  while ( mmax < n )
  {
    istep = 2 * mmax;
    theta = isign * 2.0 * ( 3.141592653589793 / mmax );
    wtemp = sin ( 0.5 * theta );
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin ( theta );
    wr = 1.0;
    wi = 0.0;
    for ( m = 1; m < mmax; m++,m++ )
    {
      for ( i = m; i <= n; i += istep )
      {
        j = i + mmax;
        tempr = wr * data[j-1] - wi * data[j];
        tempi = wr * data[j] + wi * data[j-1];
        data[j-1] = data[i-1] - tempr;
        data[j] = data[i] - tempi;
        data[i-1] = data[i-1] + tempr;
        data[i] = data[i] + tempi;
      }  
      wtemp = wr;
      wr = wr * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}
