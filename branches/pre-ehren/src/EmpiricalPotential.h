////////////////////////////////////////////////////////////////////////////////
//
// EmpiricalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EmpiricalPotential.h,v 1.6 2007/04/19 18:36:11 draeger1 Exp $

#ifndef EMPIRICALPOTENTIAL_H
#define EMPIRICALPOTENTIAL_H

#include "Context.h"
#include "D3vector.h"
#include <string>
#include <vector>
using namespace std;

class EmpiricalPotential {
  private:

  const Context& ctxt_;
  string pottype_;
  string spname1_, spname2_;
  string filename_;
  int npts_;
  vector<double> r_;
  vector<double> pot_;
  vector<double> pot_spl_;
  const double param1_;
  const double param2_;
  const double param3_;

  public:

  int is1,is2;
  EmpiricalPotential(const Context& ctxt, string pottype, string spname1, string spname2, double param1, double param2);
  EmpiricalPotential(const Context& ctxt, string pottype, string spname1, string spname2, double param1, double param2, double param3);
  EmpiricalPotential(const Context& ctxt, string spname1, string spname2, string filename);
  ~EmpiricalPotential();

  double r(int i);
  double pot(double rval);
  D3vector force(D3vector r12);  // returns force on r1 (f2 = -f1)
  string spname1(void) { return spname1_; }
  string spname2(void) { return spname2_; }
  double param1(void) { return param1_; }
  double param2(void) { return param2_; }
  double param3(void) { return param3_; }
  string pottype(void) { return pottype_; }
  string filename(void) { return filename_; }
  int npts(void) { return npts_; }
  void printsys(ostream& os) const;
  
};

#endif
