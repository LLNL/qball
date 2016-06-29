////////////////////////////////////////////////////////////////////////////////
//
// KpointCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: KpointCmd.h,v 1.16 2009/12/18 18:40:13 draeger1 Exp $

#ifndef KPOINTCMD_H
#define KPOINTCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class KpointCmd : public Cmd {
  public:

  Sample *s;

  KpointCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "kpoint"; }
  char *help_msg(void) const {
    return 
    "\n kpoint\n\n"
    " syntax: kpoint kx ky kz weight\n\n"
    "   The kpoint command defines a new k-point and adds it to the k-point list.\n\n";
  }

  int action(int argc, char **argv) {
    double kx,ky,kz,wt;
  
    // kpoint must be defined with either 3, 4, 5 or 6 arguments
    if ( argc < 4 && argc > 7) {
      if ( ui->oncoutpe() )
        cout << "<!-- use:  kpoint kx ky kz weight [unit keyword] -->" << endl;
      return 1;
    }
  
    kx = atof(argv[1]);
    ky = atof(argv[2]);
    kz = atof(argv[3]);

    wt = 1.0;
    if (argc >= 5)
      wt = atof(argv[4]);

    string units = "crystal";
    if (argc >= 6)
      units = argv[5];

    if ( !(units == "crystal" || units == "CRYSTAL" || units == "cartesian" || units == "CARTESIAN" || units == "Cartesian") ) {
      if (ui->oncoutpe() )
        cout << "<ERROR> kpoint command:  units keyword must be crystal or cartesian </ERROR>" << endl;
      return 1;
    }

    UnitCell wfcell;
    if (s->wf.refcell().volume() > 0.0) 
      wfcell = s->wf.refcell();
    else
      wfcell = s->wf.cell();

    D3vector addkpt;
    if (units == "crystal" || units == "CRYSTAL") {
      addkpt = D3vector(kx,ky,kz);
    }
    else if (units == "cartesian" || units == "CARTESIAN" || units == "Cartesian") {
      double tkx = kx;
      double tky = ky;
      double tkz = kz;
      if (argc == 7) {
        double a = atof(argv[6]);
        //double twopia = 2.0*3.14159265358979/a;
        double twopia = 2.0*M_PI/a;
        tkx *= twopia;
        tky *= twopia;
        tkz *= twopia;
      }
      D3vector a0 = wfcell.a(0);
      D3vector a1 = wfcell.a(1);
      D3vector a2 = wfcell.a(2);
      const double twopiinv = 0.5/M_PI;
      D3vector cartv(tkx,tky,tkz);
      addkpt.x = cartv*a0*twopiinv;
      addkpt.y = cartv*a1*twopiinv;
      addkpt.z = cartv*a2*twopiinv;
    }
    s->wf.add_kpoint(addkpt,wt);

    return 0;
  }
};
#endif
