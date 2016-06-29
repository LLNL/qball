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
//  TorsionConstraint.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TorsionConstraint.C,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#include "TorsionConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void TorsionConstraint::setup(const AtomSet& atoms)
{
  // find position in tau array corresponding to atom name1
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  assert(is1_>=0);
  assert(ia1_>=0);
  m1_    = atoms.species_list[is1_]->mass() * 1822.89;
  assert(m1_>0.0);
  m1_inv_ = 1.0 / m1_;

  is2_ = atoms.is(name2_);
  ia2_ = atoms.ia(name2_);
  assert(is2_>=0);
  assert(ia2_>=0);
  m2_    = atoms.species_list[is2_]->mass() * 1822.89;
  assert(m2_>0.0);
  m2_inv_ = 1.0 / m2_;

  is3_ = atoms.is(name3_);
  ia3_ = atoms.ia(name3_);
  assert(is3_>=0);
  assert(ia3_>=0);
  m3_    = atoms.species_list[is3_]->mass() * 1822.89;
  assert(m3_>0.0);
  m3_inv_ = 1.0 / m3_;

  is4_ = atoms.is(name4_);
  ia4_ = atoms.ia(name4_);
  assert(is4_>=0);
  assert(ia4_>=0);
  m4_    = atoms.species_list[is4_]->mass() * 1822.89;
  assert(m4_>0.0);
  m4_inv_ = 1.0 / m4_;
}

////////////////////////////////////////////////////////////////////////////////
void TorsionConstraint::update(double dt)
{
  angle_ += velocity_ * dt;
  if ( angle_ >  180.0 ) angle_ -= 360.0;
  if ( angle_ < -180.0 ) angle_ += 360.0;
  sin_angle_ = sin((M_PI/180.0)*angle_);
  cos_angle_ = cos((M_PI/180.0)*angle_);
}

////////////////////////////////////////////////////////////////////////////////
bool TorsionConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  const double* pr4 = &r0[is4_][3*ia4_];
  double* pr1p  = &rp[is1_][3*ia1_];
  double* pr2p  = &rp[is2_][3*ia2_];
  double* pr3p  = &rp[is3_][3*ia3_];
  double* pr4p  = &rp[is4_][3*ia4_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);
  D3vector r4(pr4);

  D3vector r1p(pr1p);
  D3vector r2p(pr2p);
  D3vector r3p(pr3p);
  D3vector r4p(pr4p);

  const double h = 0.001;
  const double fac = 0.5 / h;
  D3vector dx(h,0,0), dy(0,h,0), dz(0,0,h);

  // compute gradient at r
  D3vector g1,g2,g3,g4;
  grad_sigma(r1,r2,r3,r4,g1,g2,g3,g4);
  const double a = torsion_angle(r1,r2,r3,r4);
  double ng = g1*g1 + g2*g2 + g3*g3 + g4*g4;
  assert(ng>=0.0);

  // compute gradient at rp
  D3vector g1p,g2p,g3p,g4p;
  grad_sigma(r1p,r2p,r3p,r4p,g1p,g2p,g3p,g4p);
  const double ap = torsion_angle(r1p,r2p,r3p,r4p);
  const double ngp = g1p*g1p + g2p*g2p + g3p*g3p + g4p*g4p;
  assert(ngp>=0.0);

  const double err = fabs(ap-angle_);

  if ( ng == 0.0 )
  {
    // gradient is ill defined at r and is likely to change
    // rapidly between r and rp, invalidating the Taylor expansion
    // use the gradient at rp for the correction
    if ( ngp == 0.0 )
    {
      // gradient is undefined at r and rp
      // Insufficient information to correct rp
      return err < tol_;
    }
    else
    {
      // gradient is defined at rp but not at r
      // use gradient at rp only
      g1 = g1p;
      g2 = g2p;
      g3 = g3p;
      g4 = g4p;
      ng = ngp;
    }
  }
  else
  {
    // gradient is defined at r
    if ( ngp == 0.0 )
    {
      // gradient is defined at r but not at rp
      return err < tol_;
    }
    else
    {
      // gradient is defined at r and rp
      // no action necessary
    }
  }

  // test alignment of the gradient at r and at rp
  // compute scalar product of normalized gradients
  const double sp = ( g1*g1p + g2*g2p + g3*g3p + g4*g4p ) / sqrt( ng * ngp );

  if ( fabs(sp) < 0.5*sqrt(3.0) )
  {
    g1 = g1p;
    g2 = g2p;
    g3 = g3p;
    g4 = g4p;
  }
#if DEBUG_CONSTRAINTS
  cout << " TorsionConstraint::enforce_r: "
       << name1_ << " " << name2_ << " " << name3_ << " " << name4_ << endl;
  cout << " TorsionConstraint::enforce_r: tol = " << tol_ << endl;
  cout << " TorsionConstraint::enforce_r: ap = " << ap << endl;
  cout << " TorsionConstraint::enforce_r: err = " << err << endl;
#endif
  if ( err < tol_ ) return true;

  // make one shake iteration
  const double den = m1_inv_ * g1 * g1p +
                     m2_inv_ * g2 * g2p +
                     m3_inv_ * g3 * g3p +
                     m4_inv_ * g4 * g4p;
  assert(fabs(den)>0.0);

  const double lambda = - sigma(r1p,r2p,r3p,r4p) / den;

  pr1p[0] += m1_inv_ * lambda * g1.x;
  pr1p[1] += m1_inv_ * lambda * g1.y;
  pr1p[2] += m1_inv_ * lambda * g1.z;

  pr2p[0] += m2_inv_ * lambda * g2.x;
  pr2p[1] += m2_inv_ * lambda * g2.y;
  pr2p[2] += m2_inv_ * lambda * g2.z;

  pr3p[0] += m3_inv_ * lambda * g3.x;
  pr3p[1] += m3_inv_ * lambda * g3.y;
  pr3p[2] += m3_inv_ * lambda * g3.z;

  pr4p[0] += m4_inv_ * lambda * g4.x;
  pr4p[1] += m4_inv_ * lambda * g4.y;
  pr4p[2] += m4_inv_ * lambda * g4.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool TorsionConstraint::enforce_v(const vector<vector<double> > &r0,
  vector<vector<double> > &v0) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  const double* pr4 = &r0[is4_][3*ia4_];
  double* pv1 = &v0[is1_][3*ia1_];
  double* pv2 = &v0[is2_][3*ia2_];
  double* pv3 = &v0[is3_][3*ia3_];
  double* pv4 = &v0[is4_][3*ia4_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);
  D3vector r4(pr4);

  D3vector v1(pv1);
  D3vector v2(pv2);
  D3vector v3(pv3);
  D3vector v4(pv4);

  D3vector g1,g2,g3,g4;

  grad_sigma(r1,r2,r3,r4,g1,g2,g3,g4);

  const double norm2 = g1*g1 + g2*g2 + g3*g3 +g4*g4;

  // if the gradient is too small, do not attempt correction
  if ( norm2 < 1.e-6 ) return true;

  const double proj = v1*g1 + v2*g2 + v3*g3 + v4*g4;
  const double err = fabs(proj)/sqrt(norm2);
#if DEBUG_CONSTRAINTS
  cout << " TorsionConstraint::enforce_v: "
       << name1_ << " " << name2_ << " " << name3_ << " " << name4_<< endl;
  cout << " TorsionConstraint::enforce_v: tol = " << tol_ << endl;
  cout << " TorsionConstraint::enforce_v: err = " << err
  cout << " TorsionConstraint::enforce_v: g1  = " << g1 << endl;
  cout << " TorsionConstraint::enforce_v: g2  = " << g2 << endl;
  cout << " TorsionConstraint::enforce_v: g3  = " << g3 << endl;
  cout << " TorsionConstraint::enforce_v: g4  = " << g4 << endl;
#endif
  if ( err < tol_ ) return true;

  // make one shake iteration
  const double den = m1_inv_ * g1 * g1 +
                     m2_inv_ * g2 * g2 +
                     m3_inv_ * g3 * g3 +
                     m4_inv_ * g4 * g4;
  assert(den>0.0);

  const double eta = -proj / den;

  pv1[0] += m1_inv_ * eta * g1.x;
  pv1[1] += m1_inv_ * eta * g1.y;
  pv1[2] += m1_inv_ * eta * g1.z;

  pv2[0] += m2_inv_ * eta * g2.x;
  pv2[1] += m2_inv_ * eta * g2.y;
  pv2[2] += m2_inv_ * eta * g2.z;

  pv3[0] += m3_inv_ * eta * g3.x;
  pv3[1] += m3_inv_ * eta * g3.y;
  pv3[2] += m3_inv_ * eta * g3.z;

  pv4[0] += m4_inv_ * eta * g4.x;
  pv4[1] += m4_inv_ * eta * g4.y;
  pv4[2] += m4_inv_ * eta * g4.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void TorsionConstraint::compute_force(const vector<vector<double> > &r0,
  const vector<vector<double> > &f)
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  const double* pr4 = &r0[is4_][3*ia4_];
  const double* pf1 = &f[is1_][3*ia1_];
  const double* pf2 = &f[is2_][3*ia2_];
  const double* pf3 = &f[is3_][3*ia3_];
  const double* pf4 = &f[is4_][3*ia4_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);
  D3vector r4(pr4);

  D3vector f1(pf1);
  D3vector f2(pf2);
  D3vector f3(pf3);
  D3vector f4(pf4);

  const double h = 0.001;
  const double fac = 0.5 / h;
  D3vector dx(h,0,0), dy(0,h,0), dz(0,0,h);

  // compute gradient at r

  D3vector g1,g2,g3,g4;

  grad_sigma(r1,r2,r3,r4,g1,g2,g3,g4);

  const double norm2 = g1*g1 + g2*g2 + g3*g3 +g4*g4;
  assert(norm2>0.0);
  const double proj = f1*g1 + f2*g2 + f3*g3 + f4*g4;
  if ( norm2 == 0.0 )
  {
    force_ = 0.0;
    return;
  }
  force_ = -proj/norm2;
  // compute weight
  const double z = m1_inv_ * g1 * g1 +
                   m2_inv_ * g2 * g2 +
                   m3_inv_ * g3 * g3 +
                   m4_inv_ * g4 * g4;
  assert(z > 0.0);
  weight_ = 1.0 / sqrt(z);
}

////////////////////////////////////////////////////////////////////////////////
ostream& TorsionConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << name1_ << " ";
  os << name2_ << " " << name3_ << " " << name4_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  value=\"" << setprecision(6) << angle_;
  os << "\" velocity=\"" << setprecision(6) << velocity_ << "\"\n";
  os << "  force=\"" << setprecision(6) << force_;
  os << "\" weight=\"" << setprecision(6) << weight_ << "\"/>";
  return os;
}

////////////////////////////////////////////////////////////////////////////////
double TorsionConstraint::sigma(D3vector a, D3vector b,
 D3vector c, D3vector d) const
{
  // compute the constraint function
  return torsion_angle(a,b,c,d) - angle_;
}

////////////////////////////////////////////////////////////////////////////////
void TorsionConstraint::grad_sigma(const D3vector &r1, const D3vector &r2,
                const D3vector &r3, const D3vector &r4,
                D3vector &g1, D3vector &g2, D3vector &g3, D3vector &g4) const
{
  const double h = 0.001;
  const double fac = 0.5 / h;
  D3vector dx(h,0,0), dy(0,h,0), dz(0,0,h);

  // compute gradient at r

  g1.x = fac * ( sigma(r1+dx,r2,r3,r4) - sigma(r1-dx,r2,r3,r4) );
  g1.y = fac * ( sigma(r1+dy,r2,r3,r4) - sigma(r1-dy,r2,r3,r4) );
  g1.z = fac * ( sigma(r1+dz,r2,r3,r4) - sigma(r1-dz,r2,r3,r4) );

  g2.x = fac * ( sigma(r1,r2+dx,r3,r4) - sigma(r1,r2-dx,r3,r4) );
  g2.y = fac * ( sigma(r1,r2+dy,r3,r4) - sigma(r1,r2-dy,r3,r4) );
  g2.z = fac * ( sigma(r1,r2+dz,r3,r4) - sigma(r1,r2-dz,r3,r4) );

  g3.x = fac * ( sigma(r1,r2,r3+dx,r4) - sigma(r1,r2,r3-dx,r4) );
  g3.y = fac * ( sigma(r1,r2,r3+dy,r4) - sigma(r1,r2,r3-dy,r4) );
  g3.z = fac * ( sigma(r1,r2,r3+dz,r4) - sigma(r1,r2,r3-dz,r4) );

  g4.x = fac * ( sigma(r1,r2,r3,r4+dx) - sigma(r1,r2,r3,r4-dx) );
  g4.y = fac * ( sigma(r1,r2,r3,r4+dy) - sigma(r1,r2,r3,r4-dy) );
  g4.z = fac * ( sigma(r1,r2,r3,r4+dz) - sigma(r1,r2,r3,r4-dz) );
}

////////////////////////////////////////////////////////////////////////////////
double TorsionConstraint::torsion_angle(D3vector a, D3vector b,
 D3vector c, D3vector d) const
{
  // compute the torsion angle defined by vectors a,b,c,d
  D3vector e12(normalized(a-b));
  D3vector e32(normalized(c-b));
  D3vector e23(-e32);
  D3vector e43(normalized(d-c));
  const double sin123 = length(e12^e32);
  const double sin234 = length(e23^e43);

  double an = 0;
  if ( sin123 != 0.0 && sin234 != 0.0 )
  {
    D3vector e123 = normalized(e12^e32);
    D3vector e234 = normalized(e23^e43);
    double cc = max(min(e123*e234,1.0),-1.0);
    double ss = max(min((e123^e234)*e32,1.0),-1.0);
    an = (180.0/M_PI) * atan2(ss,cc);
  }
  return an;
}
