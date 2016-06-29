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
// LDAFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

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
