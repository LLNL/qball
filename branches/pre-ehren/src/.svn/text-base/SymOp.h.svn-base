////////////////////////////////////////////////////////////////////////////////
//
// SymOp.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SymOp.h,v 1.5 2006/10/31 22:21:48 draeger1 Exp $

#ifndef SYMOP_H
#define SYMOP_H

#include "D3vector.h"
#include "UnitCell.h"
#include <vector>
using namespace std;

class SymOp {
  private:
  
  double s11_, s12_, s13_, s21_, s22_, s23_, s31_, s32_, s33_;  // Cartesian coords
  double x11_, x12_, x13_, x21_, x22_, x23_, x31_, x32_, x33_;  // crystal coords
  vector<vector<double > > symcub_;
  vector<vector<double > > symxtal_;

  double ftrans1_, ftrans2_, ftrans3_;

  public:

  SymOp (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33);
  SymOp (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33, double ftrans1, double ftrans2, double ftrans3);

  double symxtal(const int i, const int j) const { return symxtal_[i][j]; }
  double symcub(const int i, const int j) const { return symcub_[i][j]; }

  void setFractionalTranslation(D3vector ft) { ftrans1_ = ft.x; ftrans2_ = ft.y; ftrans3_ = ft.z; };
  void convertToXtal(const UnitCell& uc);
  D3vector applyToVector(const D3vector& v, const bool applyfractrans) const;
  void applyToTensor(const double* v, double* vsym);
  void print(ostream &os);

};

ostream& operator << ( ostream &os, SymOp &s );
#endif
