////////////////////////////////////////////////////////////////////////////////
//
// Nrowmax.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Nrowmax.h,v 1.7 2009/12/18 18:40:13 draeger1 Exp $

#ifndef NROWMAX_H
#define NROWMAX_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Context.h"

class Nrowmax : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "nrowmax"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> nrowmax takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> nrowmax must be positive </ERROR>" << endl;
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
      for (int ispin = 0; ispin < s->wf.nspin(); ispin++)
        for ( int i=0; i < s->rhog_last.size(); i++ )
          cd_->rhog[ispin][i] = s->rhog_last[ispin][i];
      cd_->update_rhor();
    }
    s->wf.set_nrowmax(v);
    if ( s->wfv != 0 ) {
      s->wfv->set_nrowmax(v);
      //s->wfv->clear();
    }
    if (reshape_rhog) {
      const Context& newvctxt = s->wf.sdloc(0)->basis().context();
      cd_->reshape_rhor(oldvctxt,newvctxt);
      s->rhog_last.resize(s->wf.nspin());
      for (int ispin = 0; ispin < s->wf.nspin(); ispin++) {
        int rhogsize = cd_->rhog[ispin].size();
        s->rhog_last[ispin].resize(rhogsize);
        complex<double> *rhogp = &s->rhog_last[ispin][0];
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
     st << setw(10) << s->wf.nrowmax();
     return st.str();
  }

  Nrowmax(Sample *sample) : s(sample) {};
};
#endif
