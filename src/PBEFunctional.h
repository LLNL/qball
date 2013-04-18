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
// PBEFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PBEFUNCTIONAL_H
#define PBEFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>
using namespace std;

class PBEFunctional : public XCFunctional
{
  PBEFunctional();
  
  vector<double> _exc, _exc_up, _exc_dn;
  vector<double> _vxc1, _vxc1_up, _vxc1_dn, 
                 _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn;
  vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];
  
  void gcor2(double a, double a1, 
    double b1, double b2, double b3,
    double b4, double rtrs, double *gg, double *ggrs);
    
  void excpbe(double rho, double grad, 
    double *exc, double *vxc1, double *vxc2);
    
  void excpbe_sp(double rho_up, double rho_dn, 
    double grad_up, double grad_dn, double grad,  
    double *exc_up, double *exc_dn,
    double *vxc1_up, double *vxc1_dn, double *vxc2_upup, double *vxc2_dndn,
    double *vxc2_updn, double *vxc2_dnup);

  public:
  
  PBEFunctional(const vector<vector<double> > &rhoe);
  
  bool isGGA() { return true; };
  string name() { return "PBE"; };
  void setxc(void); 
};
#endif
