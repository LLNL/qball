////////////////////////////////////////////////////////////////////////////////
//
// LDAFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LDAFunctional.h,v 1.3 2005/11/29 18:59:47 draeger1 Exp $

#ifndef LDAFUNCTIONAL_H
#define LDAFUNCTIONAL_H

#include <vector>
#include <cassert>
using namespace std;
#include "XCFunctional.h"

class LDAFunctional : public XCFunctional {
  void xc_unpolarized(const double rh, double &ee, double &vv);
  void xc_polarized(const double rh, double &ee, double &vv);    
  vector<double> _exc;
  vector<vector<double> > _vxc;
  
  LDAFunctional();
  
  public:
  
  LDAFunctional(const vector<vector<double> > &rhoe) {
    _nspin = rhoe.size();
    if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
    _np = rhoe[0].size();
    _exc.resize(_np);
    _vxc.resize(_nspin);
    for ( int i = 0; i < _nspin; i++ )
    {
      _vxc[i].resize(_np);
    }
    
    if ( _nspin == 1 ) {
      rho = &rhoe[0][0];
      exc = &_exc[0];
      vxc1 = &_vxc[0][0];
    }
    else {
      rho_up = &rhoe[0][0];
      rho_dn = &rhoe[1][0];
      exc = &_exc[0];
      vxc1_up = &_vxc[0][0];
      vxc1_dn = &_vxc[1][0];
    }
  };
  
  bool isGGA() { return false; };
  string name() { return "LDA"; };
  void setxc(void);
};
#endif
