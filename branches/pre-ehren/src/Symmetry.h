////////////////////////////////////////////////////////////////////////////////
//
// Symmetry.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Symmetry.h,v 1.6 2007/01/23 18:11:17 draeger1 Exp $

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "D3vector.h"
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
