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
// PBEFunctional.C
//
////////////////////////////////////////////////////////////////////////////////

#include "PBEFunctional.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

PBEFunctional::PBEFunctional(const vector<vector<double> > &rhoe) {
  _nspin = rhoe.size();
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
  _np = rhoe[0].size();
  if ( _nspin == 1 ) {
    _exc.resize(_np);
    _vxc1.resize(_np);
    _vxc2.resize(_np);
    _grad_rho[0].resize(_np);
    _grad_rho[1].resize(_np);
    _grad_rho[2].resize(_np);
    rho = &rhoe[0][0];
    grad_rho[0] = &_grad_rho[0][0];
    grad_rho[1] = &_grad_rho[1][0];
    grad_rho[2] = &_grad_rho[2][0];
    exc = &_exc[0];
    vxc1 = &_vxc1[0];
    vxc2 = &_vxc2[0];
  }
  else {
    _exc_up.resize(_np);
    _exc_dn.resize(_np);
    _vxc1_up.resize(_np);
    _vxc1_dn.resize(_np);
    _vxc2_upup.resize(_np);
    _vxc2_updn.resize(_np);
    _vxc2_dnup.resize(_np);
    _vxc2_dndn.resize(_np);
    _grad_rho_up[0].resize(_np);
    _grad_rho_up[1].resize(_np);
    _grad_rho_up[2].resize(_np);
    _grad_rho_dn[0].resize(_np);
    _grad_rho_dn[1].resize(_np);
    _grad_rho_dn[2].resize(_np);
 
    rho_up = &rhoe[0][0];
    rho_dn = &rhoe[1][0];
    grad_rho_up[0] = &_grad_rho_up[0][0];
    grad_rho_up[1] = &_grad_rho_up[1][0];
    grad_rho_up[2] = &_grad_rho_up[2][0];
    grad_rho_dn[0] = &_grad_rho_dn[0][0];
    grad_rho_dn[1] = &_grad_rho_dn[1][0];
    grad_rho_dn[2] = &_grad_rho_dn[2][0];
    exc_up = &_exc_up[0];
    exc_dn = &_exc_dn[0];
    vxc1_up = &_vxc1_up[0];
    vxc1_dn = &_vxc1_dn[0];
    vxc2_upup = &_vxc2_upup[0];
    vxc2_updn = &_vxc2_updn[0];
    vxc2_dnup = &_vxc2_dnup[0];
    vxc2_dndn = &_vxc2_dndn[0];
  }
}

void PBEFunctional::setxc(void) {
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert( rho != 0 );
    assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
    assert( exc != 0 );
    assert( vxc1 != 0 );
    assert( vxc2 != 0 );
    for ( int i = 0; i < _np; i++ )
    {
      double grad = sqrt(grad_rho[0][i]*grad_rho[0][i] +
                         grad_rho[1][i]*grad_rho[1][i] +
                         grad_rho[2][i]*grad_rho[2][i] );
      excpbe(rho[i],grad,&exc[i],&vxc1[i],&vxc2[i]);
    }
  }
  else
  {
    assert( rho_up != 0 );
    assert( rho_dn != 0 );
    assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
    assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
    assert( exc_up != 0 );
    assert( exc_dn != 0 );
    assert( vxc1_up != 0 );
    assert( vxc1_dn != 0 );
    assert( vxc2_upup != 0 );
    assert( vxc2_updn != 0 );
    assert( vxc2_dnup != 0 );
    assert( vxc2_dndn != 0 );

    for ( int i = 0; i < _np; i++ )
    {
      double grx_up = grad_rho_up[0][i];
      double gry_up = grad_rho_up[1][i];
      double grz_up = grad_rho_up[2][i];
      double grx_dn = grad_rho_dn[0][i];
      double gry_dn = grad_rho_dn[1][i];
      double grz_dn = grad_rho_dn[2][i];
      double grx = grx_up + grx_dn;
      double gry = gry_up + gry_dn;
      double grz = grz_up + grz_dn;
      double grad_up = sqrt(grx_up*grx_up + gry_up*gry_up + grz_up*grz_up);
      double grad_dn = sqrt(grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn);
      double grad    = sqrt(grx*grx + gry*gry + grz*grz);
      excpbe_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&exc_up[i],&exc_dn[i],
                &vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],&vxc2_dndn[i],
                &vxc2_updn[i], &vxc2_dnup[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// 
//  excpbe: PBE exchange-correlation 
//  K.Burke's modification of PW91 codes, May 14, 1996.
//  Modified again by K.Burke, June 29, 1996, with simpler Fx(s)
//  Translated into C and modified by F.Gygi, Dec 9, 1996.
// 
//  input:
//    rho:  density
//    grad: abs(grad(rho))
//  output:
//    exc: exchange-correlation energy per electron
//    vxc1, vxc2 : quantities such that the total exchange potential is:
// 
//      vxc = vxc1 + div ( vxc2 * grad(n) )
// 
//  References:
//  [a] J.P.Perdew, K.Burke, and M.Ernzerhof,
//      "Generalized gradient approximation made simple, 
//      Phys.Rev.Lett. 77, 3865, (1996).
//  [b] J.P.Perdew and Y.Wang, Phys.Rev. B33, 8800 (1986),
//      Phys.Rev. B40, 3399 (1989) (E).
// 
////////////////////////////////////////////////////////////////////////////////
           
void PBEFunctional::excpbe(double rho, double grad, double *exc, double *vxc1, double *vxc2) {
  const double third  = 1.0 / 3.0;
  const double third4 = 4.0 / 3.0;
  const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double um = 0.2195149727645171;
  const double uk = 0.804;
  const double ul = um / uk;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
  const double alpha = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
  const double seven_sixth  =  7.0 / 6.0;
  const double four_over_pi = 1.27323954473516;
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double bet = 0.06672455060314922; /* see [a] (4) */
  const double delt = bet / gamma;

  double rtrs,fk,twoks,rs,t,h,
         ecrs,pon,b,b2,t2,t4,q4,q5,
         t6,rsthrd,fac,bec,q8,q9,hb,hrs,ht,vc;

  double rh13,exunif,s,s2,p0,fxpbe,fs;
  double ex,vx1,vx2,ec,vc1,vc2;

  *exc = 0.0;
  *vxc1 = 0.0;
  *vxc2 = 0.0;

  if ( rho < 1.e-18  ) {
    return;
  }

  /* exchange */

  rh13 = pow ( rho, third );

  /* LDA exchange energy density */
  exunif = ax * rh13;

  /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
  fk = pi32third * rh13;
  s  = grad / ( 2.0 * fk * rho );

  /* PBE enhancement factor */

  s2 = s * s;
  p0 = 1.0 + ul * s2;
  fxpbe = 1.0 + uk - uk / p0;

  ex = exunif * fxpbe;

  /* energy done, now the potential */
  /* find first derivative of Fx w.r.t the variable s. */
  /* fs = (1/s) * d Fx / d s */

  fs = 2.0 * uk * ul / ( p0 * p0 );

  vx1 = third4 * exunif * ( fxpbe - s2 * fs );
  vx2 = - exunif * fs / ( rho * 4.0 * fk * fk );

  /* correlation */

  /* Find LSD contributions, using [c] (10) and Table I of [c]. */
  /* ec = unpolarized LSD correlation energy */
  /* ecrs = d ec / d rs */
  /* construct ec, using [c] (8) */

  rs = alpha / fk;
  twoks = 2.0 * sqrt( four_over_pi * fk );
  t = grad / ( twoks * rho );

  rtrs = sqrt(rs);
  gcor2 ( 0.0310907, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294, 
          rtrs, &ec, &ecrs );

  /* LSD potential from [c] (A1) */
  /* ecrs = d ec / d rs [c] (A2) */

  vc = ec - rs * ecrs * third;

  /* PBE correlation energy */
  /* b = A of [a] (8) */

  pon = - ec / gamma;
  b = delt / ( exp ( pon ) - 1.0 );
  b2 = b * b;
  t2 = t * t;
  t4 = t2 * t2;
  q4 = 1.0 + b * t2;
  q5 = q4 + b2 * t4;
  h = gamma * log ( 1.0 + delt * q4 * t2 / q5 );

  // Energy done, now the potential, using appendix E of [b]

  t6 = t4 * t2;
  rsthrd = rs * third;
  fac = delt / b + 1.0;
  bec = b2 * fac / bet;
  q8 = q5 * q5 + delt * q4 * q5 * t2;
  q9 = 1.0 + 2.0 * b * t2;
  hb = - bet * b * t6 * ( 2.0 + b * t2 ) / q8;
  hrs = -rsthrd * hb * bec * ecrs;
  ht = 2.0 * bet * q9 / q8;

  vc1 = vc + h + hrs - t2 * ht * seven_sixth;
  vc2 = - ht / ( rho * twoks * twoks );

  *exc = ex + ec + h;
  *vxc1 = vx1 + vc1;
  *vxc2 = vx2 + vc2;
}

////////////////////////////////////////////////////////////////////////////////

void PBEFunctional::excpbe_sp(double rho_up, double rho_dn, 
  double grad_up, double grad_dn, double grad, double *exc_up, double *exc_dn,
  double *vxc1_up, double *vxc1_dn, double *vxc2_upup, double *vxc2_dndn,
  double *vxc2_updn, double *vxc2_dnup) {

  const double third  = 1.0 / 3.0;
  const double third2 =  2.0 / 3.0;
  const double third4 =  4.0 / 3.0;
  const double sixthm = -1.0 / 6.0;
  const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double um = 0.2195149727645171;
  const double uk = 0.804;
  const double ul = um / uk;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
  const double alpha = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
  const double seven_sixth  =  7.0 / 6.0;
  const double four_over_pi = 1.27323954473516;
  const double gam = 0.5198420997897463; /* gam = 2^(4/3) - 2 */
  const double fzz = 8.0 / ( 9.0 * gam );
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double bet = 0.06672455060314922; /* see [a] (4) */
  const double delt = bet / gamma;
  const double eta = 1.e-12; // small number to avoid blowup as |zeta|->1

  double eu,eurs,ep,eprs,alfm,alfrsm;
  double ex_up,ex_dn,vx1_up,vx1_dn,vx2_up,vx2_dn,ec,vc1_up,vc1_dn,vc2;

  *exc_up = 0.0;
  *exc_dn = 0.0;
  *vxc1_up = 0.0;
  *vxc1_dn = 0.0;
  *vxc2_upup = 0.0;
  *vxc2_updn = 0.0;
  *vxc2_dnup = 0.0;
  *vxc2_dndn = 0.0;

  if ( rho_up < 1.e-18 && rho_dn < 1.e-18  ) {
    return;
  }

  /* exchange up */

  ex_up = 0.0;
  vx1_up = 0.0;
  vx2_up = 0.0;
  if ( rho_up > 1.e-18 ) 
  {
    double tworho = 2.0 * rho_up;
    double gr = 2.0 * grad_up;
 
    double rh13 = pow ( tworho, third );
    /* LDA exchange energy density */
    double exunif = ax * rh13;
    /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
    double fk = pi32third * rh13;
    double s  = gr / ( 2.0 * fk * tworho );
    /* PBE enhancement factor */
    double s2 = s * s;
    double p0 = 1.0 + ul * s2;
    double fxpbe = 1.0 + uk - uk / p0;
    ex_up = exunif * fxpbe;
    /* energy done, now the potential */
    /* find first derivative of Fx w.r.t the variable s. */
    /* fs = (1/s) * d Fx / d s */
    double fs = 2.0 * uk * ul / ( p0 * p0 );
    vx1_up = third4 * exunif * ( fxpbe - s2 * fs );
    vx2_up = - exunif * fs / ( tworho * 4.0 * fk * fk );
  }

  /* exchange dn */

  ex_dn = 0.0;
  vx1_dn = 0.0;
  vx2_dn = 0.0;
  if ( rho_dn > 1.e-18 ) {
    double tworho = 2.0 * rho_dn;
    double gr = 2.0 * grad_dn;
 
    double rh13 = pow ( tworho, third );
    /* LDA exchange energy density */
    double exunif = ax * rh13;
    /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
    double fk = pi32third * rh13;
    double s  = gr / ( 2.0 * fk * tworho );
    /* PBE enhancement factor */
    double s2 = s * s;
    double p0 = 1.0 + ul * s2;
    double fxpbe = 1.0 + uk - uk / p0;
    ex_dn = exunif * fxpbe;
    /* energy done, now the potential */
    /* find first derivative of Fx w.r.t the variable s. */
    /* fs = (1/s) * d Fx / d s */
    double fs = 2.0 * uk * ul / ( p0 * p0 );
    vx1_dn = third4 * exunif * ( fxpbe - s2 * fs );
    vx2_dn = - exunif * fs / ( tworho * 4.0 * fk * fk );
  }

  /* correlation */

  // Find LSD contributions, using [c] (10) and Table I of [c].
  // eu = unpolarized LSD correlation energy
  // eurs = d eu / d rs
  // ep = fully polarized LSD correlation energy
  // eprs = d ep / d rs
  // alfm = - spin stiffness, [c] (3)
  // alfrsm = -d alpha / d rs
  // f = spin-scaling factor from [c] (9)
  // construct ec, using [c] (8)

  double rhotot = rho_up + rho_dn;
  double rh13 = pow ( rhotot, third );
  double zet = ( rho_up - rho_dn ) / rhotot;
  double g = 0.5 * ( pow(1.0+zet, third2) + pow(1.0-zet, third2) );       
  double fk = pi32third * rh13;
  double rs = alpha / fk;
  double twoksg = 2.0 * sqrt( four_over_pi * fk ) *g;
  double t = grad / ( twoksg * rhotot );

  //ewd: needed to avoid nans in spin-polarized calculations with low density regions
  if (zet > 1.0 || zet < -1.0)
    return;
  
  double rtrs = sqrt(rs);
  gcor2 ( 0.0310907, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294, 
          rtrs, &eu, &eurs );
  gcor2 ( 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, 
          rtrs, &ep, &eprs );
  gcor2 ( 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671,
          rtrs, &alfm, &alfrsm );
  double z4 = zet * zet * zet * zet;
  double f = (pow(1.0+zet,third4)+pow(1.0-zet,third4)-2.0)/gam;
  ec = eu * ( 1.0 - f * z4 ) + ep * f * z4 - alfm * f * (1.0-z4) / fzz;

  /* LSD potential from [c] (A1) */
  /* ecrs = d ec / d rs [c] (A2) */
  double ecrs = eurs * ( 1.0 - f * z4 ) + eprs * f * z4 - alfrsm * f * (1.0-z4)/fzz;
  double fz = third4 * ( pow(1.0+zet,third) - pow(1.0-zet,third))/gam;
  double eczet = 4.0 * (zet*zet*zet) * f * ( ep - eu + alfm/fzz ) +
          fz * ( z4 * ep - z4 * eu - (1.0-z4) * alfm/fzz );
  double comm = ec - rs * ecrs * third - zet * eczet;
  vc1_up = comm + eczet;
  vc1_dn = comm - eczet;

  /* PBE correlation energy */
  /* b = A of [a] (8) */

  double g3 = g * g * g;
  double pon = - ec / (g3 * gamma);
  double b = delt / ( exp ( pon ) - 1.0 );
  double b2 = b * b;
  double t2 = t * t;
  double t4 = t2 * t2;
  double q4 = 1.0 + b * t2;
  double q5 = q4 + b2 * t4;
  double h = g3 * gamma * log ( 1.0 + delt * q4 * t2 / q5 );

  /* Energy done, now the potential, using appendix E of [b] */

  double g4 = g3 * g;
  double t6 = t4 * t2;
  double rsthrd = rs * third;
  double gz = ( pow ( (1.0+zet)*(1.0+zet) + eta, sixthm ) -
         pow ( (1.0-zet)*(1.0-zet) + eta, sixthm ) ) * third;
  double fac = delt / b + 1.0;
  double bg = -3.0 * b2 * ec * fac / ( bet * g4 );
  double bec = b2 * fac / ( bet * g3 );
  double q8 = q5 * q5 + delt * q4 * q5 * t2;
  double q9 = 1.0 + 2.0 * b * t2;
  double hb = - bet * g3 * b * t6 * ( 2.0 + b * t2 ) / q8;
  double hrs = -rsthrd * hb * bec * ecrs;
  double hzed = 3.0 * gz * h / g + hb * ( bg * gz + bec * eczet );
  double ht = 2.0 * bet * g3 * q9 / q8;

  double ccomm = h + hrs - t2 * ht * seven_sixth;
  double pref = hzed - gz * t2 * ht / g;
  
  ccomm -= pref * zet;
  
  vc1_up += ccomm + pref;
  vc1_dn += ccomm - pref;
  vc2 = - ht / ( rhotot * twoksg * twoksg );

  *exc_up = ex_up + ec + h;
  *exc_dn = ex_dn + ec + h;
  *vxc1_up = vx1_up + vc1_up;
  *vxc1_dn = vx1_dn + vc1_dn;
  *vxc2_upup = 2 * vx2_up + vc2;
  *vxc2_dndn = 2 * vx2_dn + vc2;
  *vxc2_updn = vc2;
  *vxc2_dnup = vc2;
}

////////////////////////////////////////////////////////////////////////////////
// 
//  gcor2.c: Interpolate LSD correlation energy
//  as given by (10) of Perdew & Wang, Phys Rev B45 13244 (1992)
//  Translated into C by F.Gygi, Dec 9, 1996
// 
////////////////////////////////////////////////////////////////////////////////

void PBEFunctional::gcor2(double a, double a1, double b1, double b2, double b3,
  double b4, double rtrs, double *gg, double *ggrs) {

  double q0,q1,q2,q3;
  q0 = -2.0 * a * ( 1.0 + a1 * rtrs * rtrs );
  q1 = 2.0 * a * rtrs * ( b1 + rtrs * ( b2 + rtrs * ( b3 + rtrs * b4 ) ) );
  q2 = log ( 1.0 + 1.0 / q1 );
  *gg = q0 * q2;
  q3 = a * ( b1 / rtrs + 2.0 * b2 + rtrs * ( 3.0 * b3 + 4.0 * b4 * rtrs ));
  *ggrs = -2.0 * a * a1 * q2 - q0 * q3 / ( q1 * ( 1.0 + q1 ));
}
