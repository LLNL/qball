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
// LDAFunctional.C
//
// LDA Exchange-correlation energy and potential
// Ceperley & Alder, parametrized by Perdew and Zunger
//
////////////////////////////////////////////////////////////////////////////////

#include "LDAFunctional.h"
#include <cmath>
#include <cassert>
#include <vector>

void LDAFunctional::setxc(void) 
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert(rho != 0);
    assert(exc != 0);
    assert(vxc1 != 0);
    for ( int ir = 0; ir < _np; ir++ )
    {
      xc_unpolarized(rho[ir],exc[ir],vxc1[ir]);
    }
  }
  else
  {
    // spin polarized
    assert(rho_up != 0);
    assert(rho_dn != 0);
    assert(exc != 0);
    assert(vxc1_up != 0);
    assert(vxc1_dn != 0);
    const double fz_prefac = 1.0 / ( cbrt(2.0)*2.0 - 2.0 );
    const double dfz_prefac = (4.0/3.0) * fz_prefac;
    for ( int ir = 0; ir < _np; ir++ )
    {
      double excir = 0.0;
      double v_up = 0.0;
      double v_dn = 0.0;

      double roe_up = rho_up[ir];
      double roe_dn = rho_dn[ir];
      double roe = roe_up + roe_dn;

      if ( roe > 0.0 )
      {
        double zeta = ( roe_up - roe_dn ) / roe;

        double zp1 = 1.0 + zeta;
        double zm1 = 1.0 - zeta;
        double zp1_13 = cbrt(zp1);
        double zm1_13 = cbrt(zm1);
        double fz = fz_prefac * ( zp1_13 * zp1 + zm1_13 * zm1 - 2.0 );
        double dfz = dfz_prefac * ( zp1_13 - zm1_13 );

        double xc_u, xc_p, v_u, v_p;
        xc_unpolarized(roe,xc_u,v_u);
        xc_polarized(roe,xc_p,v_p);

        double xc_pu = xc_p - xc_u;
        excir = xc_u + fz * xc_pu;

        double v = v_u + fz * ( v_p - v_u );
        v_up = v + xc_pu * (  1.0 - zeta ) * dfz;
        v_dn = v + xc_pu * ( -1.0 - zeta ) * dfz;
      }

      vxc1_up[ir] = v_up;
      vxc1_dn[ir] = v_dn;
      exc[ir] = excir;
    }
  }
}

void LDAFunctional::xc_unpolarized(const double rh, double &ee, double &vv)
{
  // compute LDA xc energy and potential, unpolarized 
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711
  const double c3 = -0.610887057711;

  const double A  =  0.0311;
  const double B  = -0.048;
  const double b1 =  1.0529;
  const double b2 =  0.3334;
  const double G  = -0.1423;
  
  // C from the PZ paper: const double C  =  0.0020;
  // D from the PZ paper: const double D  = -0.0116;
  // C and D by matching Ec and Vc at rs=1
  const double D = G / ( 1.0 + b1 + b2 ) - B;
  const double C = -A - D - G * ( (b1/2.0 + b2) / ((1.0+b1+b2)*(1.0+b1+b2)));
  
  ee = 0.0;
  vv = 0.0;

  if ( rh > 0.0 )
  {
    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    double ex=0.0,vx=0.0,ec=0.0,vc=0.0;

    // Next line : exchange part in Hartree units
    vx = c3 / rs;
    ex = 0.75 * vx;
 
    // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
    if ( rs < 1.0 )
    {
      double logrs = log(rs);
      ec = A * logrs + B + C * rs * logrs + D * rs;
      vc = A * logrs + ( B - A / 3.0 ) +
           (2.0/3.0) * C * rs * logrs +
           ( ( 2.0 * D - C ) / 3.0 ) * rs;
    }
    else
    {
      double sqrtrs = sqrt(rs);
      double den = 1.0 + b1 * sqrtrs + b2 * rs;
      ec = G / den;
      vc = ec * ( 1.0 + (7.0/6.0) * b1 * sqrtrs +
                  (4.0/3.0) * b2 * rs ) / den;
    }
    ee = ex + ec;
    vv = vx + vc;
  }
}

void LDAFunctional::xc_polarized(const double rh, double &ee, double &vv)
{
  // compute LDA polarized XC energy and potential

  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711
  // const double c3 = -0.610887057711;
  // c4 = 2**third * c3 
  const double c4 = -0.769669463118;

  const double A =  0.01555;
  const double B = -0.0269;
  const double b1  =  1.3981;
  const double b2  =  0.2611;
  const double G   = -0.0843;
  // C from PZ paper: const double C   =  0.0007;
  // D from PZ paper: const double D   = -0.0048;
  // C and D by matching Ec and Vc at rs=1
  const double D = G / ( 1.0 + b1 + b2 ) - B;
  const double C = -A - D - G * ( (b1/2.0 + b2) / ((1.0+b1+b2)*(1.0+b1+b2)));

  ee = 0.0;
  vv = 0.0;

  if ( rh > 0.0 )
  {
    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    double ex=0.0,vx=0.0,ec=0.0,vc=0.0;

    // Next line : exchange part in Hartree units
    vx = c4 / rs;
    ex = 0.75 * vx;
 
    // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
    if ( rs < 1.0 )
    {
      double logrs = log(rs);
      ec = A * logrs + B + C * rs * logrs + D * rs;
      vc = A * logrs + ( B - A / 3.0 ) +
           (2.0/3.0) * C * rs * logrs +
           ( ( 2.0 * D - C ) / 3.0 ) * rs;
    }
    else
    {
      double sqrtrs = sqrt(rs);
      double den = 1.0 + b1 * sqrtrs + b2 * rs;
      ec = G / den;
      vc = ec * ( 1.0 + (7.0/6.0) * b1 * sqrtrs +
                  (4.0/3.0) * b2 * rs ) / den;
    }
    ee = ex + ec;
    vv = vx + vc;
  }
}
