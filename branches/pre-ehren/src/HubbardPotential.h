////////////////////////////////////////////////////////////////////////////////
//
// HubbardPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HubbardPotential.h,v 1.3 2010/08/26 18:39:00 draeger1 Exp $

#ifndef HUBBARDPOTENTIAL_H
#define HUBBARDPOTENTIAL_H

#include "AtomSet.h"
#include "SymmetrySet.h"
#include "Wavefunction.h"
#include "Context.h"
#include "Matrix.h"
#include "D3vector.h"

class HubbardPotential
{
  private:
  
  const Context& ctxt_;
  AtomSet& atoms_;
  const SymmetrySet& symset_;
  const Wavefunction& wf_;
  
  int nsp;   // number of species
  vector<int> na; // number of atoms of each species

  // DFT+U parameters (alpha unused for now, should be zero)
  vector<int> hub_l_;
  vector<double> hub_u_;
  vector<double> hub_alpha_;
  vector<int> lmsize_;
  vector<vector<vector<vector<double> > > > phiylm;
  vector<vector<vector<double> > > ylmsym;
  vector<D3vector> rmrand;
  vector<double> rmrandlen;
  vector<vector<vector<double> > > rmsym;
  vector<vector<double> > rmsymlen;
  
  void init(void);
   
  public:
  
  HubbardPotential(AtomSet& as, const SymmetrySet& symset, const Wavefunction& wf) :  
      wf_(wf), ctxt_(wf.context()), atoms_(as), symset_(symset) { init(); }
  ~HubbardPotential(void);
               
  void update_phiylm(void);
  double energy(bool compute_hpsi, Wavefunction& dwf);

};
#endif
