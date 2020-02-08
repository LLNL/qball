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
// VectorPotentialDynamics.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef VECTORPOTENTIALDYNAMICS_H
#define VECTORPOTENTIALDYNAMICS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>
#include <qball/VectorPotential.h>

class VectorPotentialDynamics : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "vector_potential_dynamics"; };

  int set ( int argc, char **argv ) {
    if (argc != 2){
      ui->error("vector_potential_dynamics takes only one value.");
      return 1;
    }
    
    string v = argv[1];

    if (v == "none"){
      s->ctrl.vector_potential_dynamics = VectorPotential::Dynamics::NONE;
    } else if (v == "polarization"){
      s->ctrl.vector_potential_dynamics = VectorPotential::Dynamics::POLARIZATION;
    } else {
      ui->error("vector_potential_dynamics must be \"none\" or \"polarization\"");
      return 1;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     if (s->ctrl.vector_potential_dynamics == VectorPotential::Dynamics::NONE)
       st << setw(10) << "none";
     else if (s->ctrl.vector_potential_dynamics == VectorPotential::Dynamics::POLARIZATION)
       st << setw(10) << "polarization";
     return st.str();
  }

  VectorPotentialDynamics(Sample *sample) : s(sample) { s->ctrl.vector_potential_dynamics = VectorPotential::Dynamics::NONE; };
};
#endif

// Local Variables:
// mode: c++
// End:
