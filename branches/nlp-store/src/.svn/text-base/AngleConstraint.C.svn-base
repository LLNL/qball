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
//  AngleConstraint.C
//
////////////////////////////////////////////////////////////////////////////////

#include "AngleConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void AngleConstraint::setup(const AtomSet& atoms)
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
}

////////////////////////////////////////////////////////////////////////////////
void AngleConstraint::update(double dt)
{
  //cout << " AngleConstraint::update" << endl;
  angle_ += velocity_ * dt;
  if ( angle_ <   0.0 ) angle_ =   0.0;
  if ( angle_ > 180.0 ) angle_ = 180.0;
}

////////////////////////////////////////////////////////////////////////////////
bool AngleConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  double* pr1p  = &rp[is1_][3*ia1_];
  double* pr2p  = &rp[is2_][3*ia2_];
  double* pr3p  = &rp[is3_][3*ia3_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);
  D3vector g1,g2,g3;
  grad_sigma(r1,r2,r3,g1,g2,g3);
  const double a = bond_angle(r1,r2,r3);
  double ng = g1*g1 + g2*g2 + g3*g3;
  assert(ng>=0.0);

  D3vector r1p(pr1p);
  D3vector r2p(pr2p);
  D3vector r3p(pr3p);
  D3vector g1p,g2p,g3p;
  grad_sigma(r1p,r2p,r3p,g1p,g2p,g3p);
  const double ap = bond_angle(r1p,r2p,r3p);
  const double ngp = g1p*g1p + g2p*g2p + g3p*g3p;
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
  const double sp = ( g1*g1p + g2*g2p + g3*g3p ) / sqrt( ng * ngp );

  if ( fabs(sp) < 0.5*sqrt(3.0) )
  {
    g1 = g1p;
    g2 = g2p;
    g3 = g3p;
  }

  const double den = m1_inv_ * g1 * g1p +
                     m2_inv_ * g2 * g2p +
                     m3_inv_ * g3 * g3p;

#if DEBUG_CONSTRAINTS
  cout << " AngleConstraint::enforce_r: "
       << name1_ << " " << name2_ << " " << name3_
       << " angle = " << ap << endl;
  cout << " AngleConstraint::enforce_r:  r1  = " << r1 << endl;
  cout << " AngleConstraint::enforce_r:  r2  = " << r2 << endl;
  cout << " AngleConstraint::enforce_r:  r3  = " << r3 << endl;
  cout << " AngleConstraint::enforce_r: r1p  = " << r1p << endl;
  cout << " AngleConstraint::enforce_r: r2p  = " << r2p << endl;
  cout << " AngleConstraint::enforce_r: r3p  = " << r3p << endl;
  cout << " AngleConstraint::enforce_r: ap = " << ap << endl;
  cout << " AngleConstraint::enforce_r: err = " << err << endl;
  cout << " AngleConstraint::enforce_r: tol = " << tol_ << endl;
  cout << " AngleConstraint::enforce_r: g1  = " << g1 << endl;
  cout << " AngleConstraint::enforce_r: g2  = " << g2 << endl;
  cout << " AngleConstraint::enforce_r: g3  = " << g3 << endl;
  cout << " AngleConstraint::enforce_r: g1p = " << g1p << endl;
  cout << " AngleConstraint::enforce_r: g2p = " << g2p << endl;
  cout << " AngleConstraint::enforce_r: g3p = " << g3p << endl;
  cout << " AngleConstraint::enforce_r: den = " << den << endl;
#endif
  if ( err < tol_ ) return true;

  // make one shake iteration

  assert(fabs(den)>0.0);

  const double lambda = - sigma(r1p,r2p,r3p) / den;

  pr1p[0] += m1_inv_ * lambda * g1.x;
  pr1p[1] += m1_inv_ * lambda * g1.y;
  pr1p[2] += m1_inv_ * lambda * g1.z;

  pr2p[0] += m2_inv_ * lambda * g2.x;
  pr2p[1] += m2_inv_ * lambda * g2.y;
  pr2p[2] += m2_inv_ * lambda * g2.z;

  pr3p[0] += m3_inv_ * lambda * g3.x;
  pr3p[1] += m3_inv_ * lambda * g3.y;
  pr3p[2] += m3_inv_ * lambda * g3.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool AngleConstraint::enforce_v(const vector<vector<double> > &r0,
vector<vector<double> > &v0) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  double* pv1 = &v0[is1_][3*ia1_];
  double* pv2 = &v0[is2_][3*ia2_];
  double* pv3 = &v0[is3_][3*ia3_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);

  D3vector g1,g2,g3;

  grad_sigma(r1,r2,r3,g1,g2,g3);

  D3vector v1(pv1);
  D3vector v2(pv2);
  D3vector v3(pv3);

  const double proj = v1*g1 + v2*g2 + v3*g3;
  const double norm2 = g1*g1 + g2*g2 + g3*g3;

  // if the gradient is too small, do not attempt correction
  if ( norm2 < 1.e-6 ) return true;

  const double err = fabs(proj)/sqrt(norm2);

#if DEBUG_CONSTRAINTS
  cout << " AngleConstraint::enforce_v: "
       << name1_ << " " << name2_ << " " << name3_ << endl;
  cout << " AngleConstraint::enforce_v: tol = " << tol_ << endl;
  cout << " AngleConstraint::enforce_v: err = " << err << endl;
  cout << " AngleConstraint::enforce_v: g1  = " << g1 << endl;
  cout << " AngleConstraint::enforce_v: g2  = " << g2 << endl;
  cout << " AngleConstraint::enforce_v: g3  = " << g3 << endl;
#endif
  if ( err < tol_ ) return true;

  // make one shake iteration
  const double den = m1_inv_ * g1 * g1 +
                     m2_inv_ * g2 * g2 +
                     m3_inv_ * g3 * g3;
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

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void AngleConstraint::compute_force(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  const double* pr3 = &r0[is3_][3*ia3_];
  const double* pf1 = &f[is1_][3*ia1_];
  const double* pf2 = &f[is2_][3*ia2_];
  const double* pf3 = &f[is3_][3*ia3_];

  D3vector r1(pr1);
  D3vector r2(pr2);
  D3vector r3(pr3);
  D3vector g1,g2,g3;

  grad_sigma(r1,r2,r3,g1,g2,g3);

  D3vector f1(pf1);
  D3vector f2(pf2);
  D3vector f3(pf3);

  const double norm2 = g1*g1 + g2*g2 + g3*g3;
  if ( norm2 == 0.0 )
  {
    force_ = 0.0;
    return;
  }

  const double proj = f1*g1 + f2*g2 + f3*g3;
  force_ = -proj/norm2;

  // compute weight
  const double z = m1_inv_ * g1 * g1 +
                   m2_inv_ * g2 * g2 +
                   m3_inv_ * g3 * g3;
  assert(z > 0.0);
  weight_ = 1.0 / sqrt(z);
#if DEBUG_CONSTRAINTS
  // check value of z
  const double r12s = norm(r1-r2);
  const double r32s = norm(r3-r2);
  const double fac = 180.0/M_PI;
  const double cos_theta = normalized(r1-r2)*normalized(r3-r2);
  const double zcheck = fac*fac *
    ( m1_inv_ / r12s +
      m2_inv_ * ( (r12s+r32s-2*sqrt(r12s*r32s)*cos_theta) /(r12s*r32s) ) +
      m3_inv_ / r32s
    );
  cout << " AngleConstraint: z=" << z << " zcheck=" << zcheck << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
double AngleConstraint::sigma(D3vector a, D3vector b, D3vector c) const
{
  // compute the constraint function
  return bond_angle(a,b,c) - angle_;
}

////////////////////////////////////////////////////////////////////////////////
void AngleConstraint::grad_sigma(const D3vector &r1, const D3vector &r2,
                const D3vector &r3,
                D3vector &g1, D3vector &g2, D3vector &g3) const
{
  D3vector r12(r1-r2);
  D3vector r32(r3-r2);
  assert(norm(r12) > 0.0);
  assert(norm(r32) > 0.0);
  D3vector e12(normalized(r12));
  D3vector e32(normalized(r32));
  const double ss = length(e12^e32);

  // use simplified expression if the angle is smaller than 0.2 degrees
  if ( ss < sin(0.2*(M_PI/180.0) ) )
  {
    const double eps = 1.e-8;
    if ( ss < eps )
    {
      g1 = D3vector(0.0,0.0,0.0);
      g2 = D3vector(0.0,0.0,0.0);
      g3 = D3vector(0.0,0.0,0.0);
      //cout << " ======== grad_sigma returning zero gradient" << endl;
    }
    else
    {
      // choose direction e as e12+e32
      D3vector e(e12+e32);
      assert(norm(e)>0.0);
      e = normalized(e);
      const double r12_inv = 1.0/length(r12);
      const double r32_inv = 1.0/length(r32);
      g1 = -(180.0/M_PI) * r12_inv * e;
      g2 =  (180.0/M_PI) * ( r12_inv + r32_inv ) * e;
      g3 = -(180.0/M_PI) * r32_inv * e;
      //cout << " ========= grad_sigma returning e12+e32 e =" << e << endl;
    }
  }
  else
  {
    // angle is large enough. Use finite differences
    //cout << " ========= grad_sigma using finite differences" << endl;

    const double r12_inv = 1.0/length(r12);
    const double r32_inv = 1.0/length(r32);
    const double a = bond_angle(r1,r2,r3);

    const double l12 = length(r1-r2);
    const double l32 = length(r3-r2);
    // displacement h causes changes in angle of less than 0.05 degrees
    const double h = 0.05 * (M_PI/180.0) * min(l12,l32);
    const double fac = 0.5 / h;
    D3vector dx(h,0,0), dy(0,h,0), dz(0,0,h);

    // compute gradient at r

    g1.x = fac * ( sigma(r1+dx,r2,r3) - sigma(r1-dx,r2,r3) );
    g1.y = fac * ( sigma(r1+dy,r2,r3) - sigma(r1-dy,r2,r3) );
    g1.z = fac * ( sigma(r1+dz,r2,r3) - sigma(r1-dz,r2,r3) );

    g2.x = fac * ( sigma(r1,r2+dx,r3) - sigma(r1,r2-dx,r3) );
    g2.y = fac * ( sigma(r1,r2+dy,r3) - sigma(r1,r2-dy,r3) );
    g2.z = fac * ( sigma(r1,r2+dz,r3) - sigma(r1,r2-dz,r3) );

    g3.x = fac * ( sigma(r1,r2,r3+dx) - sigma(r1,r2,r3-dx) );
    g3.y = fac * ( sigma(r1,r2,r3+dy) - sigma(r1,r2,r3-dy) );
    g3.z = fac * ( sigma(r1,r2,r3+dz) - sigma(r1,r2,r3-dz) );
  }
}

////////////////////////////////////////////////////////////////////////////////
double AngleConstraint::bond_angle(D3vector a, D3vector b, D3vector c) const
{
  // compute the bond angle defined by vectors a,b,c
  D3vector e12(normalized(a-b));
  D3vector e32(normalized(c-b));
  const double ss = length(e12^e32);
  const double cc = e12*e32;
  double an = (180.0/M_PI) * atan2(ss,cc);
  return an;
}

////////////////////////////////////////////////////////////////////////////////
ostream& AngleConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << name1_ << " ";
  os << name2_ << " " << name3_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  value=\"" << setprecision(6) << angle_;
  os << "\" velocity=\"" << setprecision(6) << velocity_ << "\"\n";
  os << "  force=\"" << setprecision(6) << force_;
  os << "\" weight=\"" << setprecision(6) << weight_ << "\"/>";
  return os;
}
