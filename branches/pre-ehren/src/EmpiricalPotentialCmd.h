////////////////////////////////////////////////////////////////////////////////
//
// EmpiricalPotentialCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EmpiricalPotentialCmd.h,v 1.6 2008/04/07 22:00:37 draeger1 Exp $

#ifndef EMPIRICALPOTENTIALCMD_H
#define EMPIRICALPOTENTIALCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"
#include "EmpiricalPotential.h"

class EmpiricalPotentialCmd : public Cmd {
  public:

  Sample *s;

  EmpiricalPotentialCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "empirical_potential"; }

  char *help_msg(void) const {
    return 
    "\n empirical_potential\n\n"
    " syntax: empirical_potential [species1] [species2] file [filename]\n"
    " syntax: empirical_potential [species1] [species2] lennard-jones [epsilon] [sigma]\n"
    " syntax: empirical_potential [species1] [species2] r12 [epsilon] [sigma]\n"
    " syntax: empirical_potential [species1] [species2] morse [D_e] [alpha] [r_e]\n\n"
    "   The empirical_potential command defines an empirical interaction between two ionic species.\n"
    "   The interaction can be read from file, with lines of the form \" r  V(r)\" on a linear grid,\n"
    "   or one of the predefined analytic forms can be used.  Current potentials are:\n\n"
    "      Lennard-Jones:    V(r) = 4*epsilon*[ (sigma/r)^12 - (sigma/r)^6]\n"
    "      r12:              V(r) = 4*epsilon*[ (sigma/r)^12 ]\n"
    "      Morse:            V(r) = D_e*[1-exp(alpha*(r-r_e))]^2\n\n"
    "   Note that potentials read from file use cubic spline interpolation with natural boundary\n"
    "   conditions (y\"=0), which can cause oscillations near r=0 for highly repulsive potentials.\n"
    "   Use a fine grid to ensure potential and forces are well-defined.\n\n";
  }

  int action(int argc, char **argv) {

    if ( argc < 5 || argc > 7 ) {
      if ( ui->oncoutpe() ) {
        cout << "<!-- use: empirical_potential [species1] [species2] file [filename] -->" << endl;
        cout << "<!-- use: empirical_potential [species1] [species2] lennard-jones [epsilon] [sigma] -->" << endl;
        cout << "<!-- use: empirical_potential [species1] [species2] r12 [epsilon] [sigma] -->" << endl;
        cout << "<!-- use: empirical_potential [species1] [species2] morse [D_e] [alpha] [r_e] -->" << endl;
      }
      return 1;
    }

    string name1 = argv[1];
    string name2 = argv[2];
    string pottype = argv[3];


    if (pottype == "L-J" || pottype == "l-j" || pottype == "Lennard-Jones" || 
        pottype == "LENNARD-JONES" || pottype == "lennard-jones") {
      double eps = atof(argv[4]);
      double sigma = atof(argv[5]);
      EmpiricalPotential* ep = new EmpiricalPotential(s->ctxt_,pottype,name1,name2,eps,sigma);
      if (!s->atoms.addEmpiricalPotential(ep)) {
        if ( ui->oncoutpe() )
          cout << "<ERROR> EmpiricalPotentialCmd: could not add " << pottype << 
            " potential between " << name1 << " and " << name2 << " </ERROR>" << endl;
        return 1;
      }
    }
    else if (pottype == "r12" || pottype == "R12") {
      double eps = atof(argv[4]);
      double sigma = atof(argv[5]);
      EmpiricalPotential* ep = new EmpiricalPotential(s->ctxt_,pottype,name1,name2,eps,sigma);
      if (!s->atoms.addEmpiricalPotential(ep)) {
        if ( ui->oncoutpe() )
          cout << "<ERROR> EmpiricalPotentialCmd: could not add " << pottype << 
            " potential between " << name1 << " and " << name2 << " </ERROR>" << endl;
        return 1;
      }
    }
    else if (pottype == "Morse" || pottype == "morse" || pottype == "MORSE") {
      double de = atof(argv[4]);
      double alpha = atof(argv[5]);
      double re = atof(argv[6]);
      EmpiricalPotential* ep = new EmpiricalPotential(s->ctxt_,pottype,name1,name2,de,alpha,re);
      if (!s->atoms.addEmpiricalPotential(ep)) {
        if ( ui->oncoutpe() )
          cout << "<ERROR> EmpiricalPotentialCmd: could not add " << pottype << 
            " potential between " << name1 << " and " << name2 << " </ERROR>" << endl;
        return 1;
      }
    }
    else if (pottype == "File" || pottype == "file" || pottype == "FILE") {
      string filename = argv[4];
      EmpiricalPotential* ep = new EmpiricalPotential(s->ctxt_,name1,name2,filename);
      if (!s->atoms.addEmpiricalPotential(ep)) {
        if ( ui->oncoutpe() )
          cout << "<ERROR> EmpiricalPotentialCmd: could not add " << pottype << 
            " potential between " << name1 << " and " << name2 << " </ERROR>" << endl;
        return 1;
      }
    }
    else {
      if ( ui->oncoutpe() ) {
        cout << "<!-- use: empirical_potential [species1] [species2] file [filename] -->" << endl;
        cout << "<!-- use: empirical_potential [species1] [species2] lennard-jones [epsilon] [sigma] -->" << endl;
        cout << "<!-- use: empirical_potential [species1] [species2] morse [D_e] [alpha] [r_e] -->" << endl;
      }
      return 1;
    }
    return 0;
  }

};
#endif
