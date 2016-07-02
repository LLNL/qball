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
// EmpiricalPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Context.h"
#include "Spline.h"
#include "EmpiricalPotential.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EmpiricalPotential::EmpiricalPotential(const Context& ctxt, string pottype, string spname1, 
string spname2, double param1, double param2): ctxt_(ctxt), pottype_(pottype), 
spname1_(spname1), spname2_(spname2), param1_(param1), param2_(param2), param3_(0.0) {

  // constructor for two-parameter empirical potentials
  is1 = -1;
  is2 = -1;

  // Lennard-Jones potential, defined by:
  //    V(r) = 4*epsilon*[ (sigma/r)^12 - (sigma/r)^6]
  // where
  //    param1 = epsilon
  //    param2 = sigma
  if (pottype_ == "L-J" || pottype_ == "l-j" || pottype_ == "Lennard-Jones" || 
      pottype_ == "LENNARD-JONES" || pottype_ == "lennard-jones") {
    pottype_ = "lennard-jones";
  }
  // Repulsive-only term from Lennard-Jones potential, defined by:
  //    V(r) = 4*epsilon*[ (sigma/r)^12 ]
  // where
  //    param1 = epsilon
  //    param2 = sigma
  else if (pottype_ == "r12" || pottype_ == "R12") {
    pottype_ = "r12";
  }
  else {
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
EmpiricalPotential::EmpiricalPotential(const Context& ctxt, string pottype, string spname1
, string spname2, double param1, double param2, double param3): ctxt_(ctxt), pottype_(pottype), 
spname1_(spname1), spname2_(spname2), param1_(param1), param2_(param2), param3_(param3) {

  // constructor for three-parameter empirical potentials
  is1 = -1;
  is2 = -1;

  // Morse potential, defined by:
  //    V(r) = D_e*[1-exp(alpha*(r-r_e))]^2
  // where
  //    param1 = D_e
  //    param2 = alpha
  //    param3 = r_e
  if (pottype_ == "Morse" || pottype_ == "morse" || pottype_ == "MORSE") {
    pottype_ = "morse";
  }
  else {
    assert(false);
  }

}

////////////////////////////////////////////////////////////////////////////////
EmpiricalPotential::EmpiricalPotential(const Context& ctxt, string filename, string spname1, 
string spname2): ctxt_(ctxt), filename_(filename), spname1_(spname1), spname2_(spname2), 
param1_(0.0), param2_(0.0), param3_(0.0) {

  // constructor for empirical potentials input on a linear grid
  is1 = -1;
  is2 = -1;

  pottype_ = "file";
  ifstream is;
  int status;
  npts_ = 0;

  if (ctxt_.oncoutpe()) {
    is.open(filename.c_str());    // text input
    status = !is.is_open();
    if (status) {
      cout << "File " << filename << " not found!" << endl; 
      filename_ = "";
    }
  }

#if USE_MPI
  MPI_Bcast(&status,1,MPI_INT,0,ctxt_.comm());
#endif

  if (ctxt_.oncoutpe()) {
    if ( !status ) {
      while (!is.eof()) {
        char readin[256]; 
        is.getline(readin,256);
        string sline(readin);

        string::size_type comment_loc = sline.find("#",0);
        if (comment_loc == string::npos) {   // no comment character in line
          istringstream iss(sline);
          vector<string> words;
          string tok;
          while (iss >> tok)
            words.push_back(tok);

          if (words.size() == 2) {
            r_.push_back(atof(words[0].c_str()));
            pot_.push_back(atof(words[1].c_str()));
            npts_++;
          }
        }
      }
      if (npts_ == 0) 
        cout << "<!-- EmpiricalPotential ERROR:  no values read from " << filename << " -->" << endl;
      else 
        cout << "<!-- EmpiricalPotential:  " << npts_ << " values read from " << filename << ":  rmin = " << r_[0] << ", rmax = " << r_[npts_-1] << " -->" << endl;
    }
  }

#if USE_MPI
  if ( !status ) {
    MPI_Bcast(&npts_,1,MPI_INT,0,ctxt_.comm());
    double tmpr[npts_];
    double tmppot[npts_];
    if (ctxt_.oncoutpe()) {
      for (int i=0; i<npts_; i++) {
        tmpr[i] = r_[i];
        tmppot[i] = pot_[i];
      }
    }
    MPI_Bcast(&tmpr[0],npts_,MPI_DOUBLE,0,ctxt_.comm());
    MPI_Bcast(&tmppot[0],npts_,MPI_DOUBLE,0,ctxt_.comm());
    if (!ctxt_.oncoutpe()) {
      r_.resize(npts_);
      pot_.resize(npts_);
      for (int i=0; i<npts_; i++) {
        r_[i] = tmpr[i];
        pot_[i] = tmppot[i];
      }
    }
  }
#endif

  // spline potential 
  if ( !status ) {
    pot_spl_.resize(npts_);
    spline(&r_[0],&pot_[0],npts_,
         SPLINE_NATURAL_BC,SPLINE_NATURAL_BC,&pot_spl_[0]);
  }
}

////////////////////////////////////////////////////////////////////////////////
EmpiricalPotential::~EmpiricalPotential(void) {

}

////////////////////////////////////////////////////////////////////////////////
double EmpiricalPotential::r(int i) {
  if (i < npts_ && i >= 0) {
    return r_[i];
  }
  else {
    return 0.0;
  }
}
////////////////////////////////////////////////////////////////////////////////
double EmpiricalPotential::pot(double rval) {

  double v = 0.0;
  if (pottype_ == "lennard-jones") {
    if (rval == 0.0) {
      v = 1.E+30;
    }
    else {
      const double x = param2_/rval;
      const double x6 = x*x*x*x*x*x;
      v = 4.*param1_*(x6*x6-x6);
    }
  }
  else if (pottype_ == "r12") {
    if (rval == 0.0) {
      v = 1.E+30;
    }
    else {
      const double x = param2_/rval;
      const double x6 = x*x*x*x*x*x;
      v = 4.*param1_*(x6*x6);
    }
  }
  else if (pottype_ == "morse") {
    const double x = 1. - exp(param2_*(rval-param3_));
    v = param1_*x*x;
  }
  else if (pottype_ == "file") {
    if (npts_ == 0 || rval > r_[npts_-1]) { 
      if (ctxt_.oncoutpe()) 
        cout << "<!-- EmpiricalPotential ERROR:  pot not defined at r = " << rval << " -->" << endl;
    }
    else {
      splint(&r_[0],&pot_[0],&pot_spl_[0],npts_,rval,&v);
    }
  }
  else {
    assert(false);
  }
  return v;
}
////////////////////////////////////////////////////////////////////////////////
D3vector EmpiricalPotential::force(D3vector r12) {

  // input vector r1 - r2, returns force on r1

  double dvdr = 0.0;

  D3vector f1;
  double rval = length(r12);

  if (rval == 0.0) {
    f1.x = 0.0; f1.y = 0.0; f1.z = 0.0;
  }
  else {
    D3vector rhat = normalized(r12);  // vector pointing from r2 to r1

    if (pottype_ == "lennard-jones") {
      const double x = param2_/rval;
      const double x6 = x*x*x*x*x*x;
      dvdr = -4.*param1_*(12.*x6*x6-6.*x6)/rval;
    }
    else if (pottype_ == "r12") {
      const double x = param2_/rval;
      const double x6 = x*x*x*x*x*x;
      dvdr = -4.*param1_*(12.*x6*x6)/rval;
    }
    else if (pottype_ == "morse") {
      const double x = exp(param2_*(rval-param3_));
      dvdr = -2.*param2_*param1_*x*(1.-x);
    }
    else if (pottype_ == "file") {
      double v;
      splintd(&r_[0],&pot_[0],&pot_spl_[0],npts_,rval,&v,&dvdr);
    }

    f1 = -dvdr*rhat;
  }
  return f1;
}

////////////////////////////////////////////////////////////////////////////////
void EmpiricalPotential::printsys(ostream& os) const {
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);

  if (pottype_ == "lennard-jones") {
    os << "empirical_potential " << spname1_ << " " << spname2_ << " lennard-jones " 
       << param1_ << " " << param2_ << endl; 
  }
  else if (pottype_ == "r12") {
    os << "empirical_potential " << spname1_ << " " << spname2_ << " r12 " 
       << param1_ << " " << param2_ << endl; 
  }
  else if (pottype_ == "morse") {
    os << "empirical_potential " << spname1_ << " " << spname2_ << " morse " 
       << param1_ << " " << param2_ << " " << param3_ << endl; 
  }    
  else if (pottype_ == "file") {
    os << "empirical_potential " << spname1_ << " " << spname2_ << " " << filename_ << endl; 
  }    
  else {
    assert(false);
  }
}
