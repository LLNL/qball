////////////////////////////////////////////////////////////////////////////////
//
//  DistanceConstraint.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: DistanceConstraint.C,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#include "DistanceConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::setup(const AtomSet& atoms)
{
  // find position in tau array corresponding to atom name1
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  assert(is1_>=0);
  assert(ia1_>=0);
  m1_    = atoms.species_list[is1_]->mass() * 1822.89;
  assert(m1_>0.0);
  m1_inv_ = 1.0/m1_;

  is2_ = atoms.is(name2_);
  ia2_ = atoms.ia(name2_);
  assert(is2_>=0);
  assert(ia2_>=0);
  m2_    = atoms.species_list[is2_]->mass() * 1822.89;
  assert(m2_>0.0);
  m2_inv_ = 1.0/m2_;
}

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::update(double dt)
{
  // cout << " DistanceConstraint::update" << endl;
  if ( distance_ + velocity_ * dt > 0.0 )
    distance_ += velocity_ * dt;
}

////////////////////////////////////////////////////////////////////////////////
bool DistanceConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  D3vector r1(pr1);
  D3vector r2(pr2);
  double* pr1p  = &rp[is1_][3*ia1_];
  double* pr2p  = &rp[is2_][3*ia2_];
  D3vector r1p(pr1p);
  D3vector r2p(pr2p);

  // compute gradient at r
  D3vector r12(r1-r2);
  D3vector g1,g2;
  g1 = 2.0 * r12;
  g2 = -g1;
  double ng = g1*g1 + g2*g2;
  assert(ng>=0.0);

  // compute gradient at rp
  D3vector r12p(r1p-r2p);
  D3vector g1p,g2p;
  g1p = 2.0 * r12p;
  g2p = -g1p;
  const double ngp = g1p*g1p + g2p*g2p;
  assert(ngp>=0.0);

  // test alignment of the gradient at r and at rp
  // compute scalar product of normalized gradients
  const double sp = ( g1*g1p + g2*g2p ) / sqrt( ng * ngp );
  if ( fabs(sp) < 0.5*sqrt(3.0) )
  {
    g1 = g1p;
    g2 = g2p;
#if DEBUG_CONSTRAINTS
    cout << " g and gp nearly orthogonal, use gp only" << endl;
#endif
  }

  const double sigma = r12p*r12p - distance_*distance_;
#if DEBUG_CONSTRAINTS
  cout << " DistanceConstraint::enforce_r: " << name1_ << " " << name2_ << endl;
  cout << " DistanceConstraint::enforce_r: r1 = " << r1 << endl;
  cout << " DistanceConstraint::enforce_r: r2 = " << r2 << endl;
  cout << " DistanceConstraint::enforce_r: tol = " << tol_ << endl;
  cout << " DistanceConstraint::enforce_r: err = " << sigma << endl;
  cout << " DistanceConstraint::enforce_r: g1  = " << g1 << endl;
  cout << " DistanceConstraint::enforce_r: g2  = " << g2 << endl;
  cout << " DistanceConstraint::enforce_r: g1p = " << g1p << endl;
  cout << " DistanceConstraint::enforce_r: g2p = " << g2p << endl;
#endif
  if ( fabs(sigma) < tol_ ) return true;

  // make one shake iteration
  const double den = m1_inv_ * g1 * g1p + m2_inv_ * g2 * g2p;

  const double lambda = -sigma / den;

  pr1p[0] += m1_inv_ * lambda * g1.x;
  pr1p[1] += m1_inv_ * lambda * g1.y;
  pr1p[2] += m1_inv_ * lambda * g1.z;

  pr2p[0] += m2_inv_ * lambda * g2.x;
  pr2p[1] += m2_inv_ * lambda * g2.y;
  pr2p[2] += m2_inv_ * lambda * g2.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool DistanceConstraint::enforce_v(const vector<vector<double> > &r0,
vector<vector<double> > &v0) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  D3vector r1(pr1);
  D3vector r2(pr2);
  double* pv1 = &v0[is1_][3*ia1_];
  double* pv2 = &v0[is2_][3*ia2_];
  D3vector v1(pv1);
  D3vector v2(pv2);

  // compute gradient at r
  D3vector r12(r1-r2);
  D3vector g1,g2;
  g1 = 2.0 * r12;
  g2 = -g1;
  const double norm2 = g1*g1 + g2*g2;
  assert(norm2>=0.0);

  // if the gradient is too small, do not attempt correction
  if ( norm2 < 1.e-6 ) return true;
  const double proj = v1 * g1 + v2 * g2;
  const double err = fabs(proj)/sqrt(norm2);
#if DEBUG_CONSTRAINTS
  cout << " DistanceConstraint::enforce_v: " << name1_ << " " << name2_ << endl;
  cout << " DistanceConstraint::enforce_v: r1 = " << r1 << endl;
  cout << " DistanceConstraint::enforce_v: r2 = " << r2 << endl;
  cout << " DistanceConstraint::enforce_v: v1 = "
  << v1[0] << " " << v1[1] << " " << v1[2] << endl;
  cout << " DistanceConstraint::enforce_v: v2 = "
  << v2[0] << " " << v2[1] << " " << v2[2] << endl;
  cout << " DistanceConstraint::enforce_v: tol = " << tol_ << endl;
  cout << " DistanceConstraint::enforce_v: err = " << err
       << endl;
#endif
  if ( err < tol_ ) return true;

  // make one shake iteration

  const double den = m1_inv_ * g1 * g1 +
                     m2_inv_ * g2 * g2;
  assert(den>0.0);
  const double eta = -proj / den;

  pv1[0] += m1_inv_ * eta * g1.x;
  pv1[1] += m1_inv_ * eta * g1.y;
  pv1[2] += m1_inv_ * eta * g1.z;

  pv2[0] += m2_inv_ * eta * g2.x;
  pv2[1] += m2_inv_ * eta * g2.y;
  pv2[2] += m2_inv_ * eta * g2.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::compute_force(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  const double* pr1 = &r0[is1_][3*ia1_];
  const double* pr2 = &r0[is2_][3*ia2_];
  D3vector r1(pr1);
  D3vector r2(pr2);
  const double* pf1 = &f[is1_][3*ia1_];
  const double* pf2 = &f[is2_][3*ia2_];
  D3vector f1(pf1);
  D3vector f2(pf2);

  // compute gradient at r
  D3vector r12(r1-r2);
  D3vector g1,g2;
  g1 = 2.0 * r12;
  g2 = -g1;
  const double norm2 = g1*g1;
  assert(norm2>=0.0);
  if ( norm2 == 0.0 )
  {
    force_ = 0.0;
    return;
  }
  const double norm = sqrt(norm2);

  D3vector e1(g1/norm);
  D3vector e2(-e1);

  const double proj1 = f1*e1;
  const double proj2 = f2*e2;

  force_ = -0.5*(proj1+proj2);
}

////////////////////////////////////////////////////////////////////////////////
ostream& DistanceConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << name1_ << " " << name2_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  value=\"" << setprecision(6) << distance_;
  os << "\" velocity=\"" << setprecision(6) << velocity_ << "\"\n";
  os << "  force=\"" << setprecision(6) << force_;
  os << "\" weight=\"" << setprecision(6) << weight_ << "\"/>";
  return os;
}


