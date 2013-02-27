////////////////////////////////////////////////////////////////////////////////
//
// WfDyn.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WfDyn.h,v 1.4 2009/10/07 18:53:22 draeger1 Exp $

#ifndef WFDYN_H
#define WFDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class WfDyn : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "wf_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> wf_dyn takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    // AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>
    // AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>
    if ( !( v == "LOCKED" || v == "SD" || v == "PSD" ||
            v == "PSDA" || v == "JD" || v == "MD" || 
            v == "TDEULER" || v == "SOTD" || v == "SORKTD" ||
            v == "FORKTD" ) )
    {
       if ( ui->oncoutpe() )
          cout << " wf_dyn must be in [LOCKED,SD,PSD,PSDA,JD,MD,TDEULER,SOTD,SORKTD,FORKTD]" << endl;
       return 1;
    }

    if (v == "TDEULER" || v == "SOTD" || v == "SORKTD" || v == "FORKTD") {
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
