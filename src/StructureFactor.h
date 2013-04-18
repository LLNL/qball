////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// StructureFactor.h
//
////////////////////////////////////////////////////////////////////////////////

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
