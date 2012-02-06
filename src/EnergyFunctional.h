////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.h,v 1.6 2010/08/26 17:44:16 draeger1 Exp $

#ifndef ENERGYFUNCTIONAL_H
#define ENERGYFUNCTIONAL_H

#include <complex>
#include <vector>
#include <valarray>
#include <map>
#include <string>
#include "ChargeDensity.h"
#include "StructureFactor.h"
#include "Timer.h"
using namespace std;

class Sample;
class Basis;
class AtomSet;
class Wavefunction;
class UnitCell;
class FourierTransform;
class XCPotential;
class NonLocalPotential;
class ConfinementPotential;
class EnthalpyFunctional;
class HubbardPotential;

typedef map<string,Timer> TimerMap;

class EnergyFunctional
{
  private:
  
  const Sample& s_;
  const Wavefunction& wf_;
  ChargeDensity& cd_;
  Basis* vbasis_;
  FourierTransform *vft;
  vector<vector<FourierTransform*> > ft;
  StructureFactor sf;
  XCPotential* xcp;
  EnthalpyFunctional* epvf;
  vector<vector<NonLocalPotential*> > nlp;
  vector<vector<ConfinementPotential*> > cfp; // cfp[ispin][ikp]
  HubbardPotential* hubp_;
  
  vector<vector<double> > vps, dvps, rhops;
  vector<complex<double> > tmp_r, vion_local_g, dvion_local_g, vlocal_g,
      rhopst, rhogt, rhoelg, vtemp;
  vector<double> ftmp;
  
  vector<vector<double> > tau0, taum, fion_esr;
  vector<double> zv_, rcps_;
  vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, econf_, eps_, enl_, ehub_, ehart_, 
         ecoul_, exc_, esr_, eself_, ets_, epv_, etotal_;
  valarray<double> sigma_ekin,sigma_econf,sigma_eps,sigma_ehart,sigma_exc,
    sigma_enl, sigma_esr, sigma;

  //ewd DEBUG
  double pwscf_ewald_;
  
  public:

  vector<vector<double> > v_r;
  vector<vector<complex<double> > > veff_g;
  mutable TimerMap tmap;
  
  double energy(bool compute_hpsi, Wavefunction& dwf,
    bool compute_forces, vector<vector<double> >& fion,
    bool compute_stress, valarray<double>& sigma);
  
  double etotal(void) const { return etotal_; }
  double ekin(void) const { return ekin_; }
  double econf(void) const { return econf_; }
  double eps(void) const { return eps_; }
  double enl(void) const { return enl_; }
  double ehart(void) const { return ehart_; }
  double ecoul(void) const { return ecoul_; }
  double exc(void) const { return exc_; }
  double esr(void) const { return esr_; }
  double eself(void) const { return eself_; }
  double ets(void) const { return ets_; }
  double epv(void) const { return epv_; }
  double ehub(void) const { return ehub_; }
  
  const ConfinementPotential *confpot(int ispin, int ikp) const { return cfp[ispin][ikp]; }
  
  void update_vhxc(void);
  
  void atoms_moved(void);
  void cell_moved(void);
  
  void print(ostream& os) const;
  void print_memory(ostream&os, double& totsum, double& locsum) const;
  
  EnergyFunctional(const Sample& s, const Wavefunction& wf, ChargeDensity& cd);
  ~EnergyFunctional();
};
ostream& operator << ( ostream& os, const EnergyFunctional& e );
#endif
