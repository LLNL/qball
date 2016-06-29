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
#include "SphericalIntegration.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SphericalIntegration::SphericalIntegration(void) {
}
////////////////////////////////////////////////////////////////////////////////
SphericalIntegration::~SphericalIntegration(void) {
}
////////////////////////////////////////////////////////////////////////////////
double SphericalIntegration::clebsch_gordan(int j1, int m1, int j2, int m2, int j3, int m3) {
  // calculates the Clebsch-Gordan coefficients defined by 
  //   Y_l1^m1*Y_l2^m2 = SUM_LM c_l1l2m1m2^LM Y_L^M
  // where Y_l^m are spherical harmonics (NOT real form)

  if ((m1+m2) != m3) return 0.0;
  if ((j1-m1) < 0) return 0.0;
  if ((j2-m2) < 0) return 0.0;
  if ((j3-m3) < 0) return 0.0;
  if ((j1-j2+j3) < 0) return 0.0;
  if ((j2-j1+j3) < 0) return 0.0;
  if ((j1+j2-j3) < 0) return 0.0;
  
  int rmax = j1+j2-j3;
  if ((j1+m1) > rmax) rmax = j1+m1;
  if ((j1-m2-j3) > rmax) rmax = j1-m2-j3;
  if ((j2+m1-j3) > rmax) rmax = j2+m1-j3;
  if ((j2-m2) > rmax) rmax = j2-m2;

  double pref = (double)((2.*j3+1.)*dfactorial(j1-j2+j3)*dfactorial(j2-j1+j3)*dfactorial(j1+j2-j3))/(double)dfactorial(j1+j2+j3+1);

  double num = (double)(dfactorial(j3+m3)*dfactorial(j3-m3)*dfactorial(j2+m2)*dfactorial(j2-m2)*dfactorial(j1+m1)*dfactorial(j1-m1));

  double sum = 0.0;
  for (int r=0; r<=rmax; r++) {
    double denom = (double)(dfactorial(r)*dfactorial(j1+j2-j3-r)*dfactorial(j1+m1-r)*dfactorial(j3-j1+m2+r)*dfactorial(j3-m1-j2+r)*dfactorial(j2-m2-r));
    if (denom > 0.0) {
      double neg = pow(-1.0,r+j1+j2-j3);
      sum += (double)neg*sqrt(pref)*sqrt(num)/denom;
    }
  }
  return sum;
}
////////////////////////////////////////////////////////////////////////////////
double SphericalIntegration::clebsch_gordan_real(int j1, int m1, int j2, int m2, int j3, int m3) {
  // calculates the real form of the Clebsch-Gordan coefficients defined by 
  //   Y_l1m1*Y_l2m2 = SUM_LM c_l1l2m1m2LM Y_L_M
  // where Y_lm are REAL spherical harmonics

  const double pi = M_PI;
  vector<double> rgrid;
  vector<double> wt;

  init_ll74grid(rgrid,wt); // high-order Lebedev-Laikov quadrature grid of unit sphere

  // integrate product of three real ylms over unit sphere to get C-G coeff
  double ylmsum = 0.0;
  for (int i=0; i<wt.size(); i++) {
    double x = rgrid[3*i+0];
    double y = rgrid[3*i+1];
    double z = rgrid[3*i+2];
    double ylm1 = ylm_real(j1,m1,x,y,z);
    double ylm2 = ylm_real(j2,m2,x,y,z);
    double ylm3 = ylm_real(j3,m3,x,y,z);
    ylmsum += wt[i]*ylm1*ylm2*ylm3;
  }
  ylmsum *= 4.*pi;

  return ylmsum;
}
////////////////////////////////////////////////////////////////////////////////
int SphericalIntegration::ifactorial(int v) {
  if (v < 0)
    return 0;
  else if (v == 0)
    return 1;
  else {
    int f = 1;
    for (int i=1; i<=v; i++)
      f*=i;
    return f;
  }
}
////////////////////////////////////////////////////////////////////////////////
double SphericalIntegration::dfactorial(int v) {
  int f = ifactorial(v);
  return (double)f;
}
////////////////////////////////////////////////////////////////////////////////
double SphericalIntegration::ylm_real(int l, int m, double gx, double gy, double gz)
{
  // real spherical harmonics up to l = 6

  const double pi = M_PI;
  const double fpi = 4.0 * pi;

  double f;
  double gnorm = sqrt(gx*gx+gy*gy+gz*gz);
  double gi = 1./gnorm;
  if (gnorm == 0.0) gi = 0.0;
  double gxi = gx*gi;
  double gyi = gy*gi;
  double gzi = gz*gi;
  double gxi2 = gxi*gxi;
  double gyi2 = gyi*gyi;
  double gzi2 = gzi*gzi;
  
  // compute real spherical harmonics 
  if (l == 0) {
    const double s14pi = sqrt(1.0/fpi);
    f = s14pi;
  }
  else if (l == 1) {
    const double s34pi = sqrt(3.0/fpi);  
    if (m == 0) f = s34pi*gzi;
    else if (m == 1) f = s34pi*gxi;
    else if (m == -1) f = s34pi*gyi;
  }
  else if (l == 2) {
    const double s54pi = sqrt(5.0/fpi);
    const double s3 = sqrt(3.0);
    if (m == 0) f = 0.5*s54pi*(2.*gzi2-gyi2-gxi2);
    else if (m == 1) f = s3*s54pi*gxi*gzi;
    else if (m == -1) f = s3*s54pi*gyi*gzi;
    else if (m == 2) f = 0.5*s3*s54pi*(gxi2-gyi2);
    else if (m == -2) f = s3*s54pi*gxi*gyi;
  }
  else if (l == 3) {
    const double s74pi = sqrt(7.0/fpi);
    const double s2132pi = sqrt(21.0/(32.*pi));
    const double s3532pi = sqrt(35.0/(32.*pi));
    const double s1054pi = sqrt(105.0/fpi);
    if (m == 0) f = 0.5*s74pi*gzi*(2.*gzi2-3.*gxi2-3.*gyi2);
    else if (m == 1) f = s2132pi*gxi*(4.*gzi2-gxi2-gyi2);
    else if (m == -1) f = s2132pi*gyi*(4.*gzi2-gxi2-gyi2);
    else if (m == 2) f = 0.5*s1054pi*gzi*(gxi2-gyi2);
    else if (m == -2) f = s1054pi*gxi*gyi*gzi;
    else if (m == 3) f = s3532pi*gxi*(gxi2-3.*gyi2);
    else if (m == -3) f = s3532pi*gyi*(3.*gxi2-gyi2);
  }
  else if (l == 4) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s14pi = sqrt(1.0/fpi);
    const double s52pi = sqrt(10.0/fpi);
    const double s54pi = sqrt(5.0/fpi);
    const double s351pi = sqrt(35.0/pi);
    const double s3532pi = sqrt(35.0/(32.*pi));
    if (m == 0) f = 0.375*s14pi*(35.*gzi2*gzi2-30.*gzi2+3.);
    else if (m == 1) f = 0.75*s52pi*gxi*gzi*(7.*gzi2-3.);
    else if (m == -1) f = 0.75*s52pi*gyi*gzi*(7.*gzi2-3.);
    else if (m == 2) f = 0.1875*s54pi*(gxi2-gyi2)*(7.*gzi2-1.);
    else if (m == -2) f = 0.375*s54pi*gxi*gyi*(7.*gzi2-1.);
    else if (m == 3) f = 0.1875*s3532pi*gxi*gzi*(gxi2-3.*gyi2);
    else if (m == -3) f = 0.1875*s3532pi*gyi*gzi*(3.*gxi2-gyi2);
    else if (m == 4) f = 0.1875*s351pi*gxi2*(gxi2-3.*gyi2)-gyi2*(3.*gxi2-gyi2);
    else if (m == -4) f = 0.75*s351pi*gxi*gyi*(gxi2-gyi2);
  }
  else if (l == 5) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s11pi = sqrt(11.0/pi);
    const double s77pi = sqrt(77.0/pi);
    const double s1654pi = sqrt(165.0/fpi);
    const double s3852pi = sqrt(385.0/(2.*pi));
    const double s3854pi = sqrt(385.0/fpi);
    const double s11554pi = sqrt(1155.0/fpi);
    if (m == 0) f = 0.0625*s11pi*(63.*gzi3*gzi2-70.*gzi3+15.*gzi);
    else if (m == 1) f = -0.125*s1654pi*gxi*(21.*gzi2*gzi2-14.*gzi2+1.);
    else if (m == -1) f = -0.125*s1654pi*gyi*(21.*gzi2*gzi2-14.*gzi2+1.);
    else if (m == 2) f = 0.125*s11554pi*(gxi2-gyi2)*(3.*gzi3-gzi);
    else if (m == -2) f = 0.5*s11554pi*gxi*gyi*(3.*gzi3-gzi);
    else if (m == 3) f = -0.0625*s3852pi*(gxi3-3.*gxi*gyi2)*(9.*gzi2-1.);
    else if (m == -3) f = -0.0625*s3852pi*(3.*gxi2*gyi-gyi3)*(9.*gzi2-1.);
    else if (m == 4) f = 0.375*s3854pi*gzi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2);
    else if (m == -4) f = 1.5*s3854pi*gzi*(gxi3*gyi-gxi*gyi3);
    else if (m == 5) f = -0.1875*s77pi*(gxi3*gxi2-10.*gxi3*gyi2+gxi*gyi2*gyi2);
    else if (m == -5) f = -0.1875*s77pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2);
  }
  else if (l == 6) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s13pi = sqrt(13./pi);
    const double s2734pi = sqrt(273./fpi);
    const double s10012pi = sqrt(1001.*0.5/pi);
    const double s13652pi = sqrt(1365.*0.5/pi);
    const double s30032pi = sqrt(3003.*0.5/pi);
    const double s914pi = sqrt(91./fpi);
    if (m == 0) f = 0.03125*s13pi*(231.*gzi3*gzi3-315.*gzi2*gzi2+105.*gzi2-5.);
    else if (m == 1) f = -0.125*s2734pi*gxi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
    else if (m == -1) f = -0.125*s2734pi*gyi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
    else if (m == 2) f = 0.015625*s13652pi*(gxi2-gyi2)*(33.*gzi2*gzi2-18.*gzi2+1.);
    else if (m == -2) f = 0.0625*s13652pi*gxi*gyi*(33.*gzi2*gzi2-18.*gzi2+1.);
    else if (m == 3) f = -0.0625*s13652pi*(gxi3-3.*gxi*gyi2)*(11.*gzi3-3.*gzi);
    else if (m == -3) f = -0.0625*s13652pi*(3.*gxi2*gyi-gyi3)*(11.*gzi3-3.*gzi);
    else if (m == 4) f = 0.1875*s914pi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2)*(11.*gzi2-1.);
    else if (m == -4) f = 0.375*s914pi*(gxi3*gyi-gxi*gyi3)*(11.*gzi2-1.);
    else if (m == 5) f = -0.1875*s10012pi*(gxi3*gxi2-10.*gxi3*gyi2+5.*gxi*gyi2*gyi2)*gzi;
    else if (m == -5) f = -0.1875*s10012pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2)*gzi;
    else if (m == 6) f = 0.03125*s30032pi*(gxi3*gxi3-15.*gxi2*gxi2*gyi2+15.*gxi2*gyi2*gyi2-gyi3*gyi3);
    else if (m == -6) f = 0.03125*s30032pi*(6.*gxi3*gxi2*gyi-20.*gxi3*gyi3+6.*gxi*gyi3*gyi2);
  }
  return f;
}
////////////////////////////////////////////////////////////////////////////////
void SphericalIntegration::init_ll74grid(vector<double> &llsph_r, vector<double> &llsph_wt) {
  //  implementation of the 74-pt quadrature grid for the unit sphere
  // 
  //   V.I. Lebedev, and D.N. Laikov, "A quadrature formula for the sphere of the 131st
  //      algebraic order of accuracy", Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.

  llsph_r.clear();
  llsph_wt.clear();

  double a = 0.0;  // default values (overridden)
  double b = 0.0;  // default values (overridden)

  double v=0.5130671797338464E-03;
  llsph_gen_oh(1,a,b,v,llsph_r,llsph_wt);

  v=0.01660406956574204;
  llsph_gen_oh(2,a,b,v,llsph_r,llsph_wt);

  v=-0.02958603896103896;
  llsph_gen_oh(3,a,b,v,llsph_r,llsph_wt);

  a=0.4803844614152614;
  v=0.02657620708215946;
  llsph_gen_oh(4,a,b,v,llsph_r,llsph_wt);

  a=0.3207726489807764;
  v=0.01652217099371571;
  llsph_gen_oh(5,a,b,v,llsph_r,llsph_wt);

}
////////////////////////////////////////////////////////////////////////////////
void SphericalIntegration::llsph_gen_oh(int code, double a, double b, double v, 
                                        vector<double>& r, vector<double>& wt) {

  if (code == 1) { 
    a = 1.0;
    r.push_back(a);  r.push_back(0.0);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(-a);  r.push_back(0.0);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(a);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(-a);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(0.0);  r.push_back(a);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(0.0);  r.push_back(-a);  wt.push_back(v);
  }
  else if (code == 2) { 
    a=sqrt(0.5);
    r.push_back(0.0);  r.push_back(a);  r.push_back(a);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(-a);  r.push_back(a);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(a);  r.push_back(-a);  wt.push_back(v);
    r.push_back(0.0);  r.push_back(-a);  r.push_back(-a);  wt.push_back(v);
    r.push_back(a);  r.push_back(0.0);  r.push_back(a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(0.0);  r.push_back(a);  wt.push_back(v);
    r.push_back(a);  r.push_back(0.0);  r.push_back(-a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(0.0);  r.push_back(-a);  wt.push_back(v);
    r.push_back(a);  r.push_back(a);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(-a);  r.push_back(a);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(a);  r.push_back(-a);  r.push_back(0.0);  wt.push_back(v);
    r.push_back(-a);  r.push_back(-a);  r.push_back(0.0);  wt.push_back(v);
  }
  else if (code == 3) {
    a = sqrt(1./3.);
    r.push_back(a);  r.push_back(a);  r.push_back(a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(a);  r.push_back(a);  wt.push_back(v);
    r.push_back(a);  r.push_back(-a);  r.push_back(a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(-a);  r.push_back(a);  wt.push_back(v);
    r.push_back(a);  r.push_back(a);  r.push_back(-a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(a);  r.push_back(-a);  wt.push_back(v);
    r.push_back(a);  r.push_back(-a);  r.push_back(-a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(-a);  r.push_back(-a);  wt.push_back(v);
  }
  else if (code == 4) { 
    a = 0.4803844614152614;
    b = sqrt(1. - 2.*a*a);
    r.push_back(a);  r.push_back(a);  r.push_back(b);  wt.push_back(v);
    r.push_back(-a);  r.push_back(a);  r.push_back(b);   wt.push_back(v);
    r.push_back(a);  r.push_back(-a);  r.push_back(b);   wt.push_back(v);
    r.push_back(-a);  r.push_back(-a);  r.push_back(b);  wt.push_back(v);
    r.push_back(a);  r.push_back(a);  r.push_back(-b);  wt.push_back(v);
    r.push_back(-a);  r.push_back(a);  r.push_back(-b);   wt.push_back(v);
    r.push_back(a);  r.push_back(-a);  r.push_back(-b);  wt.push_back(v);
    r.push_back(-a);  r.push_back(-a);  r.push_back(-b);   wt.push_back(v);
    r.push_back(a);  r.push_back(b);  r.push_back(a);   wt.push_back(v);
    r.push_back(-a);  r.push_back(b);  r.push_back(a);   wt.push_back(v);
    r.push_back(a);  r.push_back(-b);  r.push_back(a);  wt.push_back(v);
    r.push_back(-a);  r.push_back(-b);  r.push_back(a);   wt.push_back(v);
    r.push_back(a);  r.push_back(b);  r.push_back(-a);   wt.push_back(v);
    r.push_back(-a);  r.push_back(b);  r.push_back(-a);   wt.push_back(v);
    r.push_back(a);  r.push_back(-b);  r.push_back(-a);   wt.push_back(v);
    r.push_back(-a);  r.push_back(-b);  r.push_back(-a);   wt.push_back(v);
    r.push_back(b);  r.push_back(a);  r.push_back(a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(a);  r.push_back(a);   wt.push_back(v);
    r.push_back(b);  r.push_back(-a);  r.push_back(a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(-a);  r.push_back(a);   wt.push_back(v);
    r.push_back(b);  r.push_back(a);  r.push_back(-a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(a);  r.push_back(-a);   wt.push_back(v);
    r.push_back(b);  r.push_back(-a);  r.push_back(-a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(-a);  r.push_back(-a);   wt.push_back(v);
  }
  else if (code == 5) {
    b=sqrt(1.-a*a);
    r.push_back(a);  r.push_back(b);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(-a);  r.push_back(b);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(a);  r.push_back(-b);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(-a);  r.push_back(-b);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(b);  r.push_back(a);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(-b);  r.push_back(a);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(b);  r.push_back(-a);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(-b);  r.push_back(-a);  r.push_back(0.0);   wt.push_back(v);
    r.push_back(a);  r.push_back(0.0);  r.push_back(b);   wt.push_back(v);
    r.push_back(-a);  r.push_back(0.0);  r.push_back(b);   wt.push_back(v);
    r.push_back(a);  r.push_back(0.0);  r.push_back(-b);   wt.push_back(v);
    r.push_back(-a);  r.push_back(0.0);  r.push_back(-b);   wt.push_back(v);
    r.push_back(b);  r.push_back(0.0);  r.push_back(a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(0.0);  r.push_back(a);   wt.push_back(v);
    r.push_back(b);  r.push_back(0.0);  r.push_back(-a);   wt.push_back(v);
    r.push_back(-b);  r.push_back(0.0);  r.push_back(-a);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(a);  r.push_back(b);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(-a);  r.push_back(b);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(a);  r.push_back(-b);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(-a);  r.push_back(-b);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(b);  r.push_back(a);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(-b);  r.push_back(a);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(b);  r.push_back(-a);   wt.push_back(v);
    r.push_back(0.0);  r.push_back(-b);  r.push_back(-a);   wt.push_back(v);
  }
  else if (code == 6) { 
    cout << "Code 6 not implemented!" << endl;
  }
}
