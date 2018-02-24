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
// qbLink.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef QBLINK_H
#define QBLINK_H

#include <valarray>
#include <vector>
using namespace std;

class Context;
class UserInterface;
class Sample;
class SampleStepper;

class qbLink {

  private:
  
  Context* ctxt;
  UserInterface* ui;
  Sample* s;
  SampleStepper* stepper;
  bool active_;
  bool onfirstpe_;

  ofstream* qboxlog;
  streambuf* saved_cout;
  streambuf* empty_cout;

  void init(void); // intialize context, sample, UserInterface variables
  void cout_to_qboxlog(void);  // redirect cout to qboxlog filestream
  void restore_cout(void);     // restore normal cout output

  public:

  void processInputLine(string inputline); // parses string as input, e.g. "set ecut 50"
  void processInputFile(string filename); // parses commands in file

  void set_ecut(double ecut);
  void set_cell(double a0x, double a0y, double a0z, double a1x, double a1y, double a1z, double a2x, double a2y, double a2z);
  void get_cell(double &a0x, double &a0y, double &a0z, double &a1x, double &a1y, double &a1z, double &a2x, double &a2y, double &a2z);
  void set_stress(string onoff);
  void set_th_time(double time);
  void set_th_temp(double temp);
  void set_dt(double dt);
  void set_fermi_temp(double temp);
  
  void runBOSteps(int niter, int nitscf, int nite);
  void runCPSteps(int niter);
  int nsp(void);
  int nsp_mm(void);
  int na(int isp);
  int na_mm(int isp);
  double mass(int isp);
  double mass_mm(int isp);
  int atomic_number(int isp);

  void get_positions(vector<vector<double> > &r);
  void set_positions(vector<vector<double> > &r);
  void get_velocities(vector<vector<double> > &v);
  void set_velocities(vector<vector<double> > &v);
  void get_forces(vector<vector<double> > &f);
  void get_fion_ext(vector<vector<double> > &f);
  void set_fion_ext(vector<vector<double> > &f);

  double get_etotal(void);
  double get_ekin(void);
  double get_econf(void);
  double get_eps(void);
  double get_enl(void);
  double get_ehart(void);
  double get_ecoul(void);
  double get_exc(void);
  double get_esr(void);
  double get_eself(void);
  double get_ets(void);
  valarray<double> get_stress_tot(void);
  valarray<double> get_stress_kin(void);
  valarray<double> get_stress_ext(void);
  valarray<double> get_stress_eks(void);

  bool active(void) const { return active_; }
  qbLink();
  qbLink(string logfilename);
  qbLink(string logfilename, int firstpe, int lastpe);
  ~qbLink();
};
#endif

// Local Variables:
// mode: c++
// End:
