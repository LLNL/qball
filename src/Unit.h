////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade.
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
// Unit.h
//
////////////////////////////////////////////////////////////////////////////////
#ifndef UNIT_H
#define UNIT_H

#include <config.h>

#include <string>

#include "Dimensions.h"

using namespace std;

class Unit {

 public:

  // Unit conversion factors are obtained from the GNU Units program
  // executed with 15 digit precision (units -d 15).
  
  Unit(Dimensions dims, string unit_name = "unknown"):
    factor_(0.0),
    init_(false),
    dims_(dims){
    
    //convert to lower case
    std::transform(unit_name.begin(), unit_name.end(), unit_name.begin(), ::tolower);

    switch(dims_){
    case(Dimensions::one):
      set_(1.0);
      break;
      
    case(Dimensions::energy):
      if(unit_name == "hartree"      || unit_name == "ha") set_(1.0);
      if(unit_name == "rydberg"      || unit_name == "ry") set_(0.5);
      if(unit_name == "electronvolt" || unit_name == "ev") set_(0.0367493224862429);
      if(unit_name == "kelvin"       || unit_name == "k")  set_(3.1668105153288e-06);
      break;

    case(Dimensions::length):
      if(unit_name == "bohr"      || unit_name == "bohrradius" || unit_name == "b") set_(1.0);
      if(unit_name == "picometer" || unit_name == "pm")                             set_(0.0188972612544037);
      if(unit_name == "angstrom"  || unit_name == "a")                              set_(1.88972612544037);
      if(unit_name == "nanometer" || unit_name == "nm")                             set_(18.8972612544037);
      break;
      
    case(Dimensions::time):
      if(unit_name == "atomictime"  || unit_name == "hbar/hartree" || unit_name == "hbar/ha") set_(1.0);
      if(unit_name == "attosecond"  || unit_name == "as")                                     set_(0.0413413733348933);
      if(unit_name == "femtosecond" || unit_name == "fs")                                     set_(41.3413733348933);
      break;

    case(Dimensions::velocity):
      if(unit_name == "atomicvelocity"  || unit_name == "bohr*hbar/hartree" || unit_name == "b*hbar/ha") set_(1.0);
      break;
    }

  }

  template <typename T>
  T to_atomic(const T & value) const {
    return value*factor_;
  }

  bool exists() const {
    return init_;
  }
  
 private:

  void set_(const double factor){
    init_ = true;
    factor_ = factor;
  }
  
  double factor_;
  bool init_;
  Dimensions dims_;
  
};

// Local Variables:
// mode: c++
// End:

#endif
