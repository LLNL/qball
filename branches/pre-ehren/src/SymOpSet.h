////////////////////////////////////////////////////////////////////////////////
//
// SymOpSet.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SymOpSet.h,v 1.5 2006/10/31 22:21:48 draeger1 Exp $

#ifndef SYMOPSET_H
#define SYMOPSET_H

#include "D3vector.h"
#include "SymOp.h"
#include "UnitCell.h"
#include <vector>
using namespace std;

class SymOpSet {
  private:
  
  int nsym_;
  vector<SymOp*> symset_;

  public:

  SymOpSet();
  ~SymOpSet();

  void generateOps(char* name);
  void convertOpsToXtal(const UnitCell& uc);
  void printXtal(ostream& os);  
  vector<SymOp*> returnAllOps(void) const {return symset_; }
  SymOp* returnOp(int n) const {return symset_[n]; }
  int nsym(void) { return nsym_; }
  
  void addOp(SymOp* addsym);
  void removeOp(SymOp* delsym);
  void clear(void);
};

#endif
