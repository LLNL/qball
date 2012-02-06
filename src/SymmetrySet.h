////////////////////////////////////////////////////////////////////////////////
//
// SymmetrySet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SymmetrySet.h,v 1.2 2006/03/24 19:29:24 draeger1 Exp $

#ifndef SYMMETRYSET_H
#define SYMMETRYSET_H

#include "Context.h"
#include "Symmetry.h"
#include <string>
using namespace std;

class SymmetrySet {
  private:
  
  const Context& ctxt_;
  int nsym_;
  
  public:
  
  SymmetrySet(const Context& ctxt) : ctxt_(ctxt) {
    nsym_ = 0;
  }

  vector<Symmetry *> symlist; // symlist[nsym], symmetry operators

  const Context& context(void) const { return ctxt_; }
  bool addSymmetry(Symmetry *sym);
  bool delSymmetry(int i);
  bool reset(void); // remove all symmetries
  int nsym(void) const { return nsym_; };
  int size(void);
};
ostream& operator << ( ostream &os, SymmetrySet &ss );
#endif
