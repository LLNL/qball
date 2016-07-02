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

#include <config.h>
#include "ClebschGordan.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ClebschGordan::ClebschGordan(void) {
}
////////////////////////////////////////////////////////////////////////////////
ClebschGordan::~ClebschGordan(void) {
}
////////////////////////////////////////////////////////////////////////////////
double ClebschGordan::coefficient(int j1, int m1, int j2, int m2, int j3, int m3) {

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

  //cout << "rmax = " << rmax << endl;
  
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
int ClebschGordan::ifactorial(int v) {
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
double ClebschGordan::dfactorial(int v) {
  int f = ifactorial(v);
  return (double)f;
}
