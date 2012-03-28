////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.h,v 1.7 2009/09/25 23:18:11 draeger1 Exp $

#ifndef CHARGEDENSITY_H
#define CHARGEDENSITY_H

#include <vector>
#include <valarray>
#include <complex>
#include <string>
#include <map>
#include "Timer.h"
#include "Context.h"

class Sample;
class Wavefunction;
class FourierTransform;
class Basis;
class AtomSet;

typedef map<string,Timer> TimerMap;

class ChargeDensity {
  private:
  
  const Context& ctxt_;
  Context vcontext_;
  Wavefunction& wf_;
  AtomSet& atoms_;
  Basis* vbasis_;
  FourierTransform* vft_;
  vector<vector<FourierTransform*> > ft_; // ft_[ispin][ikp];
  valarray<complex<double> > rhotmp;
  vector<int> symindexloc;
  vector<int> symmultloc;
  int nsym_;
  int nsymgrp_;
  bool ultrasoft_;
  bool highmem_;
  bool nlcc_;
  vector<vector<complex<double> > > qnmg_;
  vector<double> rhornlcc_; 
  
  public:
  
  mutable TimerMap tmap;

  vector<vector<double> > rhor; // rhor[ispin][i]
  vector<vector<complex<double> > > rhog; // rhog[ispin][ig]
  vector<vector<double> > xcrhor; 
  vector<complex<double> > rhognlcc; 
  vector<vector<complex<double> > > xcrhog; 

  void update_density();
  void update_rhor(void);
  
  Basis* vbasis(void) const { return vbasis_; }
  const Context& vcontext(void) const { return vcontext_; }
  FourierTransform* vft(void) const { return vft_; }
  FourierTransform* ft(int ispin, int ikp) const { return ft_[ispin][ikp]; }
  void reshape_rhor(const Context& oldvctxt,const Context& newvctxt);
  void update_usfns();
  void set_nlcc(bool nlcc) { nlcc_ = nlcc; }
  bool nlcc() { return nlcc_; }
  void update_nlcc();
  void add_nlccden();
  void nlcc_forceden(int is, vector<complex<double> > &rhog);
  void print_memory(ostream&os, double& totsum, double& locsum) const;
  void print_timing();
  
  ChargeDensity(Sample& s);
  ~ChargeDensity();
};
#endif
