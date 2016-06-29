////////////////////////////////////////////////////////////////////////////////
//
// ParallelOptimizer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ParallelOptimizer.h,v 1.3 2010/01/07 18:01:48 draeger1 Exp $

#ifndef PARALLELOPTIMIZER_H
#define PARALLELOPTIMIZER_H

#include "Sample.h"
#include "Timer.h"
#include <map>
using namespace std;

typedef map<string,Timer> TimerMap;

class ParallelOptimizer {
  private:

  Sample& s_;
  int niter_, nitscf_, nite_;
  vector<vector<double> > fion;
  valarray<double> sigma_eks, sigma_kin, sigma_ext, sigma;
  
  // Do not allow construction of ParallelOptimizer unrelated to a Sample
  ParallelOptimizer(void);

  public:

  mutable TimerMap tmap;
  
  void optimize(int niter, int nitscf, int nite);
  //double runtime(bool print_timing);
  double runtime(int nrowmax, int npark, int nspin, bool reshape, bool print_timing);
  list<int> factorize(int n);
  list<int> factorize(int n1, int n2);
  ParallelOptimizer(Sample& s);
  ~ParallelOptimizer();
};
#endif
