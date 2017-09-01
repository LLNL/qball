////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// LIBXCFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LIBXCFunctional.h,v Yi Yao $
// 
// usage: 
// LDA: set xc LIBXC LDA_X LDA_C_PZ_MOD
// PBE: set xc LIBXC GGA_X_PBE:1.0 GGA_C_PBE:1.0
// BLYP: set xc LIBXC GGA_X_B88:1.0 GGA_C_LYP:1.0
// PBESOL: set xc LIBXC GGA_X_PBE_SOL:1.0 GGA_C_PBE_SOL:1.0
// PBEREV: set xc LIBXC GGA_X_PBE_R:1.0 GGA_C_PBE:1.0

#include <config.h>

#ifndef LIBXCFUNCTIONAL_H
#define LIBXCFUNCTIONAL_H

#include <vector>
#include <cassert>
#include <sstream>
#include <iostream>
#include <cstdlib>
using namespace std;
#include "XCFunctional.h"

#ifdef HAVE_LIBXC
#include "xc.h"
#endif

class LIBXCFunctional : public XCFunctional
{
  int _xcfamily;
  // LDA
  vector<double> _exc;
  vector<vector<double> > _vxc;

  // GGA
  vector<double> _exc_up, _exc_dn;
  vector<double> _vxc1, _vxc1_up, _vxc1_dn,
            _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn;
  vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];


#ifdef HAVE_LIBXC
  xc_func_type func_;
#endif
  vector<int> func_ids_;
  vector<int> func_families_;
  vector<double> func_coeffs_;

  LIBXCFunctional();

  public:

  LIBXCFunctional(const vector<vector<double> > &rhoe, 
    const string functional_full_name)
  {
#ifdef HAVE_LIBXC
    string tempbuf;
    stringstream ss(functional_full_name);
    string functional_name_temp;
    double functional_coeff_temp;
    int n = 0;

    _xcfamily = XC_FAMILY_LDA;
    while ( ss >> tempbuf )
    {
      if (n > 0)
      {
        functional_name_temp =
          tempbuf.substr(0, tempbuf.find(":"));
        functional_coeff_temp =
          atof(tempbuf.substr( tempbuf.find(":")+1,
                                          tempbuf.length()).c_str());
        if (functional_name_temp != "HF") {
          func_ids_.push_back(xc_functional_get_number(functional_name_temp.c_str()));
          func_coeffs_.push_back(functional_coeff_temp);
        }
      }
      n++;
    }
    for (int i = 0; i < func_ids_.size();i ++) 
      func_families_.push_back(xc_family_from_id(func_ids_[i],NULL,NULL));

    for (int i = 0; i < func_ids_.size();i ++) {
      if  (func_families_[i] == XC_FAMILY_GGA 
           || func_families_[i] == XC_FAMILY_HYB_GGA
          ) 
          _xcfamily = XC_FAMILY_GGA;
    }
    for (int i = 0; i < func_ids_.size();i ++) {
      if  ( func_families_[i] == XC_FAMILY_MGGA 
           || func_families_[i] == XC_FAMILY_HYB_MGGA
          ) 
          _xcfamily = XC_FAMILY_MGGA;
    }

    switch (_xcfamily)
    {
    case XC_FAMILY_LDA:
    {
      _nspin = rhoe.size();
      if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
      _np = rhoe[0].size();
      _exc.resize(_np);
      _vxc.resize(_nspin);
      for ( int i = 0; i < _nspin; i++ )
      {
        _vxc[i].resize(_np);
      }

      if ( _nspin == 1 )
      {
        rho = &rhoe[0][0];
        exc = &_exc[0];
        vxc1 = &_vxc[0][0];
      }
      else
      {
        rho_up = &rhoe[0][0];
        rho_dn = &rhoe[1][0];
        exc = &_exc[0];
        vxc1_up = &_vxc[0][0];
        vxc1_dn = &_vxc[1][0];
      }
    }
    break;
    case XC_FAMILY_GGA:
      {
        {
          _nspin = rhoe.size();
          if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
          _np = rhoe[0].size();
        
          if ( _nspin == 1 )
          {
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
          else
          {
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

      }
      break;
    }
#endif
  };

#ifdef HAVE_LIBXC
  bool isGGA() { if (_xcfamily == XC_FAMILY_LDA) 
                         return false; 
                       else
                         return true;
                     };
#endif

  string name() { return "LIBXC"; };
  void setxc(void);
};
#endif

