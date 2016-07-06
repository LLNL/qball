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
// Hugoniostat.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Hugoniostat.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Hugoniostat::Hugoniostat(double etotref, double vref, double pref, double temp, bool oncoutpe): 
ref_etot_(etotref), ref_volume_(vref), ref_pressure_(pref), th_temp_(temp), oncoutpe_(oncoutpe) {
  if (oncoutpe_) 
    cout << "<!-- Hugoniot object initialized:  Eref = " << ref_etot_ << ", Pref = " << ref_pressure_ << ", Vref = " << ref_volume_ << " -->" << endl;

  updatenow = false;
  updatefreq_ = 5;    // default initial updatefrequency
  deltatemp_ = 10.0;  // default initial deltatemp
  firststep_ = true;

  if (oncoutpe_) 
    cout << "<!-- Hugoniotstat:  temperature will update every " << updatefreq_ << " steps. -->" << endl;
}
////////////////////////////////////////////////////////////////////////////////
Hugoniostat::~Hugoniostat(void) {

}
////////////////////////////////////////////////////////////////////////////////
void Hugoniostat::set_updatefreq(int updatefreq) { 
  updatefreq_ = updatefreq;
  if (oncoutpe_) 
    cout << "<!-- Hugoniostat:  temperature will update every " << updatefreq_ << " steps. -->" << endl;
}
////////////////////////////////////////////////////////////////////////////////
void Hugoniostat::set_deltatemp(double deltatemp) { 
  deltatemp_ = deltatemp;
  if (oncoutpe_) 
    cout << "<!-- Hugoniostat:  deltaT set to " << deltatemp << " K. -->" << endl;
}
////////////////////////////////////////////////////////////////////////////////
void Hugoniostat::addValues(double etot, double volume, double pressure, double temp) {

  etot_history.push_back(etot);
  vol_history.push_back(volume);
  press_history.push_back(pressure);
  temp_history.push_back(temp);
  double hug = etot-ref_etot_ + 0.5*(volume-ref_volume_)*(pressure-ref_pressure_);
  hug_history.push_back(hug);

  if (oncoutpe_) 
    cout << "<!-- Hugoniostat:  Hugoniot minimization factor =  " << hug << "  -->" << endl;

  if (hug_history.size() > updatefreq_)
    updatenow = true;

}
////////////////////////////////////////////////////////////////////////////////
void Hugoniostat::updateTemp(double& temperature) {

  double hugendpt = hug_history[hug_history.size()-1];
  double hugavg = 0.0;
  for (int i = 0; i < hug_history.size(); i++)
    hugavg += hug_history[i];
  hugavg /= hug_history.size();

  if (firststep_) {
    firststep_ = false;
  }
  else {
    // if most recent value of Hugoniot minimization value changed sign, do nothing
    if (! sameSign(hugavg,hugendpt)) {
      if (oncoutpe_) 
        cout << "<!-- Hugoniostat:  thermostat temperature unchanged, th_temp =  " << temperature << "  -->" << endl;
    }
    // if average value is different sign than previous average value, decrease deltaT and change sign
    else if (! sameSign(hugavg,lasthugavg_)) {
      deltatemp_ *= -0.5;
      temperature += deltatemp_;
      if (oncoutpe_) 
        cout << "<!-- Hugoniostat:  thermostat temperature changed, th_temp =  " << temperature << ", deltatemp = " << deltatemp_ << "  -->" << endl;
    }
    // haven't hit zero yet, keep changing temperature
    else if ( sameSign(hugavg,lasthugavg_)) {
      if (fabs(hugavg) > fabs(lasthugavg_))
        deltatemp_ *= -1.0;
      
      temperature += deltatemp_;
      if (oncoutpe_) 
        cout << "<!-- Hugoniostat:  thermostat temperature changed, th_temp =  " << temperature << ", deltatemp = " << deltatemp_ << "  -->" << endl;
    }
  }  

  hug_history.clear();
  etot_history.clear();
  press_history.clear();
  vol_history.clear();
  temp_history.clear();

  lasthugavg_ = hugavg;
  updatenow = false;
}

////////////////////////////////////////////////////////////////////////////////
bool Hugoniostat::sameSign(double a, double b) {
  return (a < 0.0) == (b < 0.0);
}
