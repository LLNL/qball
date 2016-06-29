////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.h,v 1.2 2010/01/07 18:01:48 draeger1 Exp $

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Sample;
class EnergyFunctional;
class Wavefunction;

#include <vector>
#include <valarray>
using namespace std;

class Preconditioner
{
  private:
  
  const Sample& s_;
  const EnergyFunctional& ef_;
  const Wavefunction& wf_;
  vector<vector<valarray<double> > > diag_; // diag_[ispin][ikp][ig]

  public:

  void update(void);
  
  const valarray<double>& diag(int ispin, int ikp) const
  { return diag_[ispin][ikp]; }

  Preconditioner(const Sample& s, const Wavefunction& wf, const EnergyFunctional& ef);
  //~Preconditioner();
};
#endif
