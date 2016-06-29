////////////////////////////////////////////////////////////////////////////////
//
// XCFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XCFunctional.h,v 1.3 2005/11/29 18:59:47 draeger1 Exp $

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
};
#endif
