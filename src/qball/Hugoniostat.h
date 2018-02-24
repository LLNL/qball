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
// Hugoniostat.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef HUGONIOSTAT_H
#define HUGONIOSTAT_H

#include <vector>
using namespace std;

class Hugoniostat {
  private:

  double ref_etot_;
  double ref_pressure_;
  double ref_volume_;
  double th_temp_;
  double deltatemp_;
  double lasthugavg_;
  bool firststep_;

  vector<double> temp_history;
  vector<double> etot_history;
  vector<double> press_history;
  vector<double> vol_history;
  vector<double> hug_history;

  int updatefreq_;
  bool oncoutpe_;
  bool sameSign(double a, double b);

  public:

  bool updatenow;

  Hugoniostat(double etotref, double vref, double pref, double temp, bool oncoutpe);
  ~Hugoniostat();
  void set_updatefreq(int updatefreq);
  void set_deltatemp(double deltatemp);
  void addValues(double etot, double volume, double pressure, double temp);
  void updateTemp(double& temperature);
};

#endif

// Local Variables:
// mode: c++
// End:
