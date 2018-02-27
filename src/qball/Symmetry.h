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
// Symmetry.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <math/D3vector.h>
#include <string>
using namespace std;

class Symmetry {
  private:
  
  double s11_, s12_, s13_, s21_, s22_, s23_, s31_, s32_, s33_;
  double ftrans1_, ftrans2_, ftrans3_;
  int ft1_,ft2_,ft3_;
  int np1_,np2_,np3_;

  public:

  Symmetry (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33);
  Symmetry (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33, double ftrans1, double ftrans2, double ftrans3);

  int setGrid(int np1, int np2, int np3);
  void applyToGridPoint(int i, int j, int k, int &outi, int &outj, int &outk);
  D3vector applyToVector(const D3vector& v, const bool applyfractrans) const;
  void applyToTensor(const double* v, double* vsym);

};

ostream& operator << ( ostream &os, Symmetry &s );
#endif

// Local Variables:
// mode: c++
// End:
