////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef NONLOCALPOTENTIAL_H
#define NONLOCALPOTENTIAL_H

#include "AtomSet.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "Context.h"
#include "Matrix.h"

class StructureFactor;

class NonLocalPotential
{
  private:
  
  const Context& ctxt_;
  AtomSet& atoms_;
  SlaterDet& sd_;
  const Basis& basis_;
  Basis* cdbasis_;  // pointer to charge density basis, needed for ultrasoft
  
  int nsp;          // number of species
  int nspnl;        // number of non-local species
  bool ultrasoft_;  // there are ultrasoft potentials
  bool highmem_; 
  
  vector<int>             lmax;     // lmax[is]
  vector<int>             lloc;     // lloc[is]
  vector<int>             na;       // na[is]
  vector<int>             npr;      // npr[is]
  vector<int>             nprna;    // nprna[is]
  vector<vector<int> >    lproj;    // lproj[is][ipr]
  vector<vector<double> > wt;       // wt[is][ipr]
  vector<vector<double> > twnl;     // twnl[is][npr*ngwl]
  vector<vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl], ij=0,..,5
  
  vector<int>             nquad;    // nquad[is]
  vector<vector<double> > rquad;    // rquad[is][iquad], iquad = 0, nquad[is]-1
  vector<vector<double> > wquad;    // wquad[is][iquad], iquad = 0, nquad[is]-1

  vector<vector<complex<double> > > qnmg_;   // ultrasoft function Q_nm^I(G), wf basis
  vector<vector<complex<double> > > sfactcd_;  // structure factor of local atoms, cd basis
  vector<vector<complex<double> > > sfactwf_;  // structure factor of local atoms, wf basis
  
  
  mutable TimerMap tmap;
  void init(void);
   
  public:
  
  NonLocalPotential(AtomSet& as, SlaterDet& sd) :  
    ctxt_(sd.context()), atoms_(as), sd_(sd), basis_(sd.basis()) { init(); }
  ~NonLocalPotential(void);
               
  void update_twnl(void);
  void update_usfns(Basis* cdbasis);  // update Q_nm^I(G), beta^I(G) when atoms
                                      // move or basis changes
  void use_highmem(void) { highmem_ = true; }  // use extra memory to speed calculation
  double energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl,
    vector<complex<double> >& veff);

  void print_memory(ostream&os, int kmult, int kmultloc, double& totsum, double& locsum) const;
};
#endif
