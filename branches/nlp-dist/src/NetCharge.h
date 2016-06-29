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
// NetCharge.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NetCharge.h,v 1.4 2008-09-08 15:56:18 fgygi Exp $

#ifndef NETCHARGE_H
#define NETCHARGE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class NetCharge : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "net_charge"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " net_charge takes only one value" << endl;
      return 1;
    }

    const int v = atoi(argv[1]);

    // compute the current netcharge
    // Definition: wf_nel = atoms_nel - netcharge
    // Note: the sign of netcharge is negative if extra electrons are present
    const int netcharge_before = s->atoms.nel() - s->wf.nel();
    if ( v == netcharge_before )
      return 0;

    // set new netcharge to v
    if ( s->atoms.nel() - v < 0 )
    {
      if ( ui->oncoutpe() )
        cout << " net_charge: cannot remove more than "
             << s->atoms.nel() << " electrons" << endl;
      return 1;
    }

    s->wf.set_nel(s->atoms.nel() - v);
    s->wf.update_occ(0.0,0);
    if ( s->wfv != 0 )
    {
      s->wfv->set_nel(s->atoms.nel() - v);
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->atoms.nel() - s->wf.nel();
     return st.str();
  }

  NetCharge(Sample *sample) : s(sample) {};
};
#endif
