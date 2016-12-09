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
// UserInterface.h:
//
////////////////////////////////////////////////////////////////////////////////
// $ Id: $

#include <config.h>

#ifndef STANDARDVAR_H
#define STANDARDVAR_H

#include <string>
#include "Dimensions.h"

using namespace std;

#include "UserInterface.h"

class StandardVar : public Var {
public:

  virtual const char * name() const { return name_.c_str(); }
  const Dimensions & dimensions() const { return dims_; }

  typedef unsigned Attributes;

  static const Attributes any_value    = 0;
  static const Attributes non_negative = 1 << 0;
  
protected:
  
  StandardVar(const string & name, Dimensions dims = Dimensions::one, const string & default_unit_name = "none"):
    name_(name),
    dims_(dims),
    default_unit_name_(default_unit_name){
  }

  int parse(int argc, char **argv, double & value, Attributes attr = 0){
    
    string unit_name;
    
    if ( argc == 2 ) {
      unit_name = default_unit_name_;
      ui->warning("Units missing for variable '" + string(name()) + "', assuming '" + default_unit_name_ + "'.");
    } else if ( argc != 3 ) {
      ui->error("The variable '" + string(name()) + "' takes two arguments: the value followed by its units.");
      return 1;
    } else {
      unit_name = argv[2];
    }
    
    Unit unit(dimensions(), unit_name);

    if(!unit.exists()) {
      ui->error("Unknown energy unit '" + unit_name + "'.");
      return 1;
    }
    
    value = unit.to_atomic(atof(argv[1]));

    if ( attr & non_negative && value < 0.0 ) {
      ui->error("The variable '" + string(name()) + "' must be non-negative.");
      return 1;
    }
    
    return 0;
  }  

private:
  string name_;
  Dimensions dims_;
  string default_unit_name_;

  
};

// Local Variables:
// mode: c++
// End:

#endif
