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

using namespace std;

class Unit {

 public:

  // Unit conversion factors are obtained from the GNU Units program
  // executed with 16 digit precision (units -d 16).
  
  static Unit Energy(string unit_name){
    //convert to lower case
    std::transform(unit_name.begin(), unit_name.end(), unit_name.begin(), ::tolower);
    
    if(unit_name == "hartree"      || unit_name == "ha") return Unit(1.0);
    if(unit_name == "rydberg"      || unit_name == "ry") return Unit(0.5);
    if(unit_name == "electronvolt" || unit_name == "ev") return Unit(0.0367493224862429);

    return Unit();
  }
  
  static Unit Length(string unit_name){
    //convert to lower case
    std::transform(unit_name.begin(), unit_name.end(), unit_name.begin(), ::tolower);

    if(unit_name == "bohr"      || unit_name == "bohrradius" || unit_name == "b") return Unit(1.0);
    if(unit_name == "angstrom"  || unit_name == "a")                              return Unit(1.88972612544037);
    if(unit_name == "nanometer" || unit_name == "nm")                             return Unit(18.8972612544037);
    if(unit_name == "picometer" || unit_name == "pm")                             return Unit(0.0188972612544037);
    
    return Unit();
  }

  template <typename T>
  T to_atomic(const T & value) const {
    return value*factor_;
  }

  bool exists() const {
    return init_;
  }
  
 private:

  Unit():factor_(0.0), init_(false){
  }
  
  Unit(double factor):factor_(factor), init_(true){
  }
  
  double factor_;
  bool init_;
  
};

// Local Variables:
// mode: c++
// End:

#endif
