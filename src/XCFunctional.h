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
// XCFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

//
// Abstract base class for density functionals
// Input variables are: rho, rho_up, rho_dn, grad_rho, grad_rho_up, grad_rho_dn
//
// Output quantities:
// The exchange-correlation energy is expressed as
//
// Exc = int{ rho[i] * exc[i] + rho_up[i] * exc_up[i] + rho_dn[i] * exc_dn[i] }
//
// It is assumed that exchange correlation potentials can be 
// written as:
// 
// vxc_up = vxc1 + vxc1_up +
//          div ( vxc2_upup grad_rho_up ) + div ( vxc2_updn grad_rho_dn )
// 
// vxc_dn = vxc1 + vxc1_dn + 
//          div ( vxc2_dndn grad_rho_dn ) + div ( vxc2_dnup grad_rho_up )
//
// Not all input quantities are needed, and not all output quantities are
// computed by certain functionals. Example:
//
// LDAFunctional:
//   without spin: Input: rho
//                 Output: exc, vxc1
//   with spin:    Input: rho_up, rho_dn
//                 Output: vxc1_up, vxc1_dn
//
// PBEFunctional:
//   without spin: Input:  rho, grad_rho
//                 Output: exc, vxc1, vxc2
//   with spin:    Input:  rho_up, rho_dn, grad_rho_up, grad_rho_dn,
//                 Output: exc_up, exc_dn, vxc1_up, vxc1_dn,
//                         vxc2_upup, vxc2_dndn, vxc2_updn, vxc2_dnup

#ifndef XCFUNCTIONAL_H
#define XCFUNCTIONAL_H

#include <string>
using namespace std;

class XCFunctional {

  protected:

  int _np, _nspin;
 
  public:
  
  const double *rho, *rho_up, *rho_dn;
  double *grad_rho[3], *grad_rho_up[3], *grad_rho_dn[3];
  double *exc, *exc_up, *exc_dn;
  double *vxc1, *vxc1_up, *vxc1_dn;
  double *vxc2, *vxc2_upup, *vxc2_dndn, *vxc2_updn, *vxc2_dnup;

  virtual bool isGGA(void) = 0;
  virtual string name(void) = 0;
  int np(void) { return _np; };
  int nspin(void) { return _nspin; };
  
  XCFunctional() {
    rho = rho_up = rho_dn = 0;
    grad_rho[0] = grad_rho[1] = grad_rho[2] = 0;
    grad_rho_up[0] = grad_rho_up[1] = grad_rho_up[2] = 0;
    grad_rho_dn[0] = grad_rho_dn[1] = grad_rho_dn[2] = 0;
    exc = exc_up = exc_dn = 0;
    vxc1 = vxc1_up = vxc1_dn = 0;
    vxc2 = vxc2_upup = vxc2_dndn = vxc2_updn = vxc2_dnup = 0;
  }

  // virtual destructor needed to ensure proper deallocation
  virtual ~XCFunctional() {}
  
  virtual void setxc(void) = 0; 

  virtual void setxc(int start, int end) = 0; 
};
#endif
