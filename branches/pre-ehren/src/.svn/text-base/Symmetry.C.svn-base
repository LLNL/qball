////////////////////////////////////////////////////////////////////////////////
//
// Symmetry.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Symmetry.C,v 1.6 2007/01/23 18:11:17 draeger1 Exp $

#include "Symmetry.h"
#include <iomanip>
using namespace std;

Symmetry::Symmetry (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33) {
  s11_ = s11;  s12_ = s12;  s13_ = s13;
  s21_ = s21;  s22_ = s22;  s23_ = s23;
  s31_ = s31;  s32_ = s32;  s33_ = s33;
  ftrans1_ = 0.;  ftrans2_ = 0.;  ftrans3_ = 0.;
  np1_ = -1; np2_ = -1; np3_ = -1; 
}

Symmetry::Symmetry (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33, double ftrans1, double ftrans2, double ftrans3) {
  s11_ = s11;  s12_ = s12;  s13_ = s13;
  s21_ = s21;  s22_ = s22;  s23_ = s23;
  s31_ = s31;  s32_ = s32;  s33_ = s33;
  ftrans1_ = ftrans1;  ftrans2_ = ftrans2;  ftrans3_ = ftrans3;
  np1_ = -1; np2_ = -1; np3_ = -1; 
}

int Symmetry::setGrid(int np1, int np2, int np3) {
   const double threshold = 1.E-5;
   np1_ = np1;
   np2_ = np2;
   np3_ = np3;
   // assuming grid points between 0 and np0,np1,np2
   ft1_ = (int)(ftrans1_*(np1_ + threshold) );
   ft2_ = (int)(ftrans2_*(np2_ + threshold) );
   ft3_ = (int)(ftrans3_*(np3_ + threshold) );

  // check that fractional translations match grid
  if (ft1_ != 0 && abs(ftrans1_*np1 - (double)ft1_) > threshold)
    return 1;
  if (ft2_ != 0 && abs(ftrans2_*np2 - (double)ft2_) > threshold)
    return 1;
  if (ft3_ != 0 && abs(ftrans3_*np3 - (double)ft3_) > threshold)
    return 1;
  return 0;
}
void Symmetry::applyToGridPoint(int i, int j, int k, int &outi, int &outj, int &outk) {
  assert (np1_ > 0);
  assert (np2_ > 0);
  assert (np3_ > 0);

  outi = (int)s11_*i + (int)s21_*j*np1_/np2_ + (int)s31_*k*np1_/np3_ - ft1_;
  outj = (int)s12_*i*np2_/np1_ + (int)s22_*j + (int)s32_*k*np2_/np3_ - ft2_;
  outk = (int)s13_*i*np3_/np1_ + (int)s23_*j*np3_/np2_ + (int)s33_*k - ft3_;

  outi = (outi%np1_);
  outj = (outj%np2_);
  outk = (outk%np3_);

  if (outi < 0) outi += np1_;
  if (outj < 0) outj += np2_;
  if (outk < 0) outk += np3_;
  return;
}

D3vector Symmetry::applyToVector(const D3vector& v, const bool applyfractrans) const {
  D3vector vsym;
  // we use this to identify symmetry-equivalent atoms
  if (applyfractrans) { 
    vsym.x = s11_*v.x + s21_*v.y + s31_*v.z - ftrans1_;
    vsym.y = s12_*v.x + s22_*v.y + s32_*v.z - ftrans2_;
    vsym.z = s13_*v.x + s23_*v.y + s33_*v.z - ftrans3_;
  }
  // we use this to average forces over symmetry equivalent atoms
  else {
    vsym.x = s11_*v.x + s21_*v.y + s31_*v.z;
    vsym.y = s12_*v.x + s22_*v.y + s32_*v.z;
    vsym.z = s13_*v.x + s23_*v.y + s33_*v.z;
  }
  return vsym;
}

void Symmetry::applyToTensor(const double* v, double* vsym) {
  // v[0]=v_xx, v[1]=v_yy, v[2]=v_zz, v[3]=v_xy, v[4]=v_yz, v[5]=v_xz, 
  //vsym_ij = sum_kl s(i,k)*v(k,l)*s(j,l) = 
  //  si1_*sj1*v[0] + si1_*sj2*v[3] + si1_*sj3*v[5] + 
  //  si2_*sj1*v[3] + si2_*sj2*v[2] + si2_*sj3*v[4] + 
  //  si3_*sj1*v[5] + si3_*sj2*v[4] + si3_*sj3*v[3];

  vsym[0] = 
    s11_*s11_*v[0] + s11_*s12_*v[3] + s11_*s13_*v[5] + 
    s12_*s11_*v[3] + s12_*s12_*v[1] + s12_*s13_*v[4] + 
    s13_*s11_*v[5] + s13_*s12_*v[4] + s13_*s13_*v[2];
  vsym[1] = 
    s21_*s21_*v[0] + s21_*s22_*v[3] + s21_*s23_*v[5] + 
    s22_*s21_*v[3] + s22_*s22_*v[1] + s22_*s23_*v[4] + 
    s23_*s21_*v[5] + s23_*s22_*v[4] + s23_*s23_*v[2];
  vsym[2] = 
    s31_*s31_*v[0] + s31_*s32_*v[3] + s31_*s33_*v[5] + 
    s32_*s31_*v[3] + s32_*s32_*v[1] + s32_*s33_*v[4] + 
    s33_*s31_*v[5] + s33_*s32_*v[4] + s33_*s33_*v[2];
  vsym[3] = 
    s11_*s21_*v[0] + s11_*s22_*v[3] + s11_*s23_*v[5] + 
    s12_*s21_*v[3] + s12_*s22_*v[1] + s12_*s23_*v[4] + 
    s13_*s21_*v[5] + s13_*s22_*v[4] + s13_*s23_*v[2];
  vsym[4] = 
    s21_*s31_*v[0] + s21_*s32_*v[3] + s21_*s33_*v[5] + 
    s22_*s31_*v[3] + s22_*s32_*v[1] + s22_*s33_*v[4] + 
    s23_*s31_*v[5] + s23_*s32_*v[4] + s23_*s33_*v[2];
  vsym[5] = 
    s11_*s31_*v[0] + s11_*s32_*v[3] + s11_*s33_*v[5] + 
    s12_*s31_*v[3] + s12_*s32_*v[1] + s12_*s33_*v[4] + 
    s13_*s31_*v[5] + s13_*s32_*v[4] + s13_*s33_*v[2];

  return;
}

ostream& operator << ( ostream &os, Symmetry &sym ) {
  os.setf(ios::left,ios::adjustfield);
  os << "<symmetry>output not implemented yet</symmetry>\n";

  //os << "  <symmetry fractional_translation = \"" << ftrans1_ << "  " << ftrans2_ << "  " << ftrans3_ << "\">\n";
  //os.setf(ios::fixed,ios::floatfield);
  //os.setf(ios::right,ios::adjustfield);
  //os << "    " << s11_ << "  " << s12_ << "  " << s13_ << "\n";
  //os << "    " << s21_ << "  " << s22_ << "  " << s23_ << "\n";
  //os << "    " << s31_ << "  " << s32_ << "  " << s33_ << "\n";
  //os << "  </symmetry>\n";

  return os;
}
