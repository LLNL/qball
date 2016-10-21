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
// Nparallelkpts.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef NPARALLELKPTS_H
#define NPARALLELKPTS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Nparallelkpts : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "nparallelkpts"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> nparallelkpts takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> nparallelkpts must be positive </ERROR>" << endl;
      return 1;
    }

    bool reshape_rhog = false;
    if (s->rhog_last.size() > 0)
      if (s->rhog_last[0].size() > 0)
        if (s->wf.hasdata())
          reshape_rhog = true;

    ChargeDensity* cd_ = 0;
    Context oldvctxt;
    if (reshape_rhog) {
      // load mixed charge density into cd_
      cd_ = new ChargeDensity(*s);
      oldvctxt = cd_->vcontext();
      for (int ispin = 0; ispin < s->wf.nspin(); ispin++) {
        for (unsigned i=0; i < s->rhog_last.size(); i++ )
          cd_->rhog[ispin][i] = s->rhog_last[ispin][i];
        cd_->update_rhor();
      }
    }
    s->wf.set_nparallelkpts(v);
    if ( s->wfv != 0 ) {
      s->wfv->set_nparallelkpts(v);
    }
    if (reshape_rhog) {
      const Context& newvctxt = s->wf.sdloc(0)->basis().context();
      cd_->reshape_rhor(oldvctxt,newvctxt);
      s->rhog_last.resize(s->wf.nspin());
      for (int ispin = 0; ispin < s->wf.nspin(); ispin++) {
        int rhogsize = cd_->rhog[ispin].size();
        s->rhog_last[ispin].resize(rhogsize);
        complex<double> *rhogpold = &cd_->rhog[ispin][0];
        for (int j = 0; j < rhogsize; j++)
          s->rhog_last[ispin][j] = rhogpold[j];
      }
      if (cd_ != 0)
        delete cd_;
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nparallelkpts();
     return st.str();
  }

  Nparallelkpts(Sample *sample) : s(sample) {};
};
#endif
