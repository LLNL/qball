#ifndef PSEUDO_UPF_HPP
#define PSEUDO_UPF_HPP

/*
 Copyright (C) 2018 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <fstream>
#include <vector>
#include <cassert>
#include <sstream>
#include <iostream>
#include <cmath>

#include "anygrid.hpp"

namespace pseudopotential {

  class upf : public pseudopotential::anygrid {

  public:

    upf(bool uniform_grid):
      pseudopotential::anygrid(uniform_grid){
    }
    
    double d_ij(int l, int i, int j) const {
      assert(l >= 0 && l <= lmax_);
      assert(i >= 0 && i <= nchannels());
      assert(j >= 0 && j <= nchannels());

      return dij_[l*nchannels()*nchannels() + i*nchannels() + j];
    }

  protected:

    double & d_ij(int l, int i, int j) {
      assert(l >= 0 && l <= lmax_);
      assert(i >= 0 && i <= nchannels());
      assert(j >= 0 && j <= nchannels());

      return dij_[l*nchannels()*nchannels() + i*nchannels() + j];
    }
    
    void extrapolate_first_point(std::vector<double> & function_) const{

      assert(function_.size() >= 4);
      assert(grid_.size() >= 4);
      
      double x1 = grid_[1];
      double x2 = grid_[2];
      double x3 = grid_[3];
      double f1 = function_[1];
      double f2 = function_[2];
      double f3 = function_[3];


      // obtained from:
      // http://www.wolframalpha.com/input/?i=solve+%7Bb*x1%5E2+%2B+c*x1+%2B+d+%3D%3D+f1,++b*x2%5E2+%2B+c*x2+%2B+d+%3D%3D+f2,+b*x3%5E2+%2B+c*x3+%2B+d+%3D%3D+f3+%7D++for+b,+c,+d
      
      function_[0] = f1*x2*x3*(x2 - x3) + f2*x1*x3*(x3 - x1) + f3*x1*x2*(x1 - x2);
      function_[0] /= (x1 - x2)*(x1 - x3)*(x2 - x3);

    }

    std::vector<double> dij_;
    int llocal_;

  };

}

#endif
