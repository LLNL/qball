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
// WfDyn.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef WFDYN_H
#define WFDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>
#include <qball/Wavefunction.h>
#include <qball/SlaterDet.h>

class WfDyn : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "wf_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      // AK: for PETSC, accept integrator type as later args 
      string v = argv[1];
      if ( v != "PETSC" )
      {
        if ( ui->oncoutpe() )
          cout << " <ERROR> wf_dyn takes only one value </ERROR>" << endl;
        return 1;
      }
    }
    
    string v = argv[1];
    // AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>
    // AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>
    if ( !( v == "LOCKED"  ||
	    v == "SD"      ||
	    v == "PSD"     ||
            v == "PSDA"    ||
	    v == "RMMDIIS" ||
	    v == "JD"      ||
	    v == "MD"      || 
            v == "TDEULER" ||
	    v == "SOTD"    ||
	    v == "SORKTD"  ||
            v == "FORKTD"  ||
	    v == "ETRS"    ||
	    v == "AETRS"   ||
            v == "PETSC") )
    {
       if ( ui->oncoutpe() )
          cout << " wf_dyn must be in [LOCKED,SD,PSD,PSDA,RMMDIIS,JD,MD,TDEULER,SOTD,SORKTD,FORKTD,ETRS,AETRS,PETSC]" << endl;
       return 1;
    }

    if (v == "TDEULER" || v == "SOTD" || v == "SORKTD" || v == "FORKTD" || v == "ETRS" || v == "AETRS" || v == "PETSC") {
       s->ctrl.tddft_involved = true;
       if (!( s->wf.force_complex_set() )) {
          cout << "WfDyn::wave functions must be complex to propagate them in time" << endl
               << "WfDyn::wf_dyn will NOT be set to " << v << " and remains " << s->ctrl.wf_dyn << " instead!" << endl;
       }
       // AS: delete wavefunction velocities if present, since SOTD uses wfv to store |psi(t-tddt)>
       if ( s->wfv != 0 )
          delete s->wfv;
       s->wfv = 0;
    }

    s->ctrl.wf_dyn = v;


    // AK: set PETSC time stepper type and subtype
    if (v=="PETSC")
    {
      if (argc < 3)
      {
        if ( ui->oncoutpe() )
          cout << " <ERROR> wf_dyn PETSC needs more values to set the time stepper type</ERROR>" << endl;
        return 1;        
      }

      string ts_type = argv[2];

      // AK: check for valid time stepper type. not comprehensive
      if ( !( ts_type == "rk" || ts_type == "euler" || ts_type == "beuler" || ts_type == "pseudo" ||
              ts_type == "cn" || ts_type == "sundials" || ts_type == "python" || ts_type == "theta" || ts_type == "alpha" ||
              ts_type == "alpha2" || ts_type == "glle" || ts_type == "glee" || ts_type == "ssp" || ts_type == "arkimex" ||
              ts_type == "rosw" || ts_type == "eimex" || ts_type == "mimex" || ts_type == "bdf" || ts_type == "radau5" ) )
      {
        if ( ui->oncoutpe() )
          cout << " <ERROR> invalid PETSc time stepper type</ERROR>" << endl;
        return 1;
      }

      s->ctrl.petsc_ts_type = ts_type;

      // AK: if there is an additional value, use it as subtype 
      // To do: input validation
      if (argc > 3)
        s->ctrl.petsc_ts_subtype = argv[3];
      else
        s->ctrl.petsc_ts_subtype="";
    }

    
//     if ( v == "MD" )
//     {
//       if ( s->wfv == 0 )
//       {
//         s->wfv = new Wavefunction(s->wf);
//         s->wfv->clear();
//       }
//     }
//     else
//     {
//       delete s->wfv;
//       s->wfv = 0;
//     }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_dyn;
     return st.str();
  }

  WfDyn(Sample *sample) : s(sample) { s->ctrl.wf_dyn = "SD"; s->ctrl.tddft_involved = false;};
};
#endif

// Local Variables:
// mode: c++
// End:
