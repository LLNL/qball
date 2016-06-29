////////////////////////////////////////////////////////////////////////////////
//
// MMSpecies.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MMSpecies.h,v 1.1 2007/04/10 23:30:12 draeger1 Exp $

#ifndef MMSPECIES_H
#define MMSPECIES_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;
#include "Context.h"

class MMSpecies {
  private:
  
  const Context& ctxt_;
  string name_;         // name used to refer to species in current application
  double mass_;        // mass in a.m.u (Carbon = 12.0)

  public:

  MMSpecies(const Context& ctxt, string name, double mass);
  
  const Context& context(void) const { return ctxt_; }
  const string& name(void) const { return name_; }
  double mass(void) const { return mass_; }

  void info(ostream& os);
  void printsys(ostream& os) const;
    
};
ostream& operator << ( ostream &os, MMSpecies &a );
#endif
