////////////////////////////////////////////////////////////////////////////////
//
// StructureFactor.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StructureFactor.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef STRUCTUREFACTOR_H
#define STRUCTUREFACTOR_H

#include <vector>
#include <complex>
#include <iostream>
using namespace std;

class Basis;

class StructureFactor
{
  private:
  
  int _nsp, _ng;
  vector<int> _na;
  
  int _k0max, _k1max, _k2max, 
      _k0min, _k1min, _k2min,
      _k0range, _k1range, _k2range;
   
  public:
  
  // convenience pointer access functions:
  // double *c0 = cos0_ptr(is,ia);
  // c0[ kx ] == cos(-i gx*tau[is][ia].x)
    
  double *cos0_ptr(int is, int ia) 
    { return &cos0[is][ia*_k0range-_k0min]; }
    
  double *cos1_ptr(int is, int ia) 
    { return &cos1[is][ia*_k1range-_k1min]; }
    
  double *cos2_ptr(int is, int ia) 
    { return &cos2[is][ia*_k2range-_k2min]; }
    
  double *sin0_ptr(int is, int ia) 
    { return &sin0[is][ia*_k0range-_k0min]; }
    
  double *sin1_ptr(int is, int ia) 
    { return &sin1[is][ia*_k1range-_k1min]; }
    
  double *sin2_ptr(int is, int ia) 
    { return &sin2[is][ia*_k2range-_k2min]; }

  // kx in [k0min, k0max]
  // ky in [k1min, k1max]
  // kz in [k2min, k2max]

  vector<vector<double> > cos0;  // cos0[is][ia*k0range-k0min+kx]
  vector<vector<double> > cos1;  // cos1[is][ia*k1range-k1min+ky]
  vector<vector<double> > cos2;  // cos2[is][ia*k2range-k2min+ky]
  vector<vector<double> > sin0;  // sin0[is][ia*k0range-k0min+kx]
  vector<vector<double> > sin1;  // sin1[is][ia*k1range-k1min+ky]
  vector<vector<double> > sin2;  // sin2[is][ia*k2range-k2min+ky]
  vector<vector<complex<double> > > sfac;  // sfac[is][ig]
  
  void init(const vector<vector<double> >& tau, const Basis& basis);
  void update(const vector<vector<double> >& tau, const Basis& basis);

};
#endif
