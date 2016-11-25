////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade (xavier@llnl.gov), Erik Draeger
// (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
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
// Extrapolator.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SlaterDet.h"
#include "ExtrapolatorSimple.h"

using namespace std;

ExtrapolatorSimple::ExtrapolatorSimple(){
}

void ExtrapolatorSimple::extrapolate_wavefunction(string extrap, Wavefunction & wf, Wavefunction* wfv, Wavefunction* wfmm, int iter, double dt, const Context& ctxt){

  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      if (wf.spinactive(ispin))
	{
	  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
            {
              if (wf.kptactive(ikp))
		{
		  assert(wf.sd(ispin,ikp) != 0);
		  if (extrap == "SIMPLE")
		    {
		      if ( ctxt.mype()==0 )
			cout << "Extrapolating wavefunction using simple algorithm." << endl;
                   
		      double* c = (double*) wf.sd(ispin,ikp)->c().cvalptr();
		      double* cv = (double*) wfv->sd(ispin,ikp)->c().cvalptr();
		      const int mloc = wf.sd(ispin,ikp)->c().mloc();
		      const int nloc = wf.sd(ispin,ikp)->c().nloc();
		      const int len = 2*mloc*nloc;
		      if ( iter == 0 )
			{
			  // copy c to cv
			  for ( int i = 0; i < len; i++ )
			    {
			      const double x = c[i];
			      const double v = cv[i];
			      c[i] = x + dt * v;
			      cv[i] = x;
			    }
			  //tmap["lowdin"].start();
			  //wf.sd(ispin,ikp)->lowdin();
			  //tmap["lowdin"].stop();
			  if (wf.ultrasoft()) {
			    //tmap["usfns"].start();
			    wf.sd(ispin,ikp)->update_usfns();
			    //tmap["usfns"].stop();
			  }
			  //tmap["gram"].start();
			  wf.sd(ispin,ikp)->gram();
			  //tmap["gram"].stop();
			}
		      else
			{
			  //tmap["align"].start();
			  wfv->align(wf);
			  //tmap["align"].stop();
                  
			  // linear extrapolation
			  for ( int i = 0; i < len; i++ )
			    {
			      const double x = c[i];
			      const double xm = cv[i];
			      c[i] = 2.0 * x - xm;
			      cv[i] = x;
			    }
			  //tmap["ortho_align"].start();
			  //wf.sd(ispin,ikp)->ortho_align(*wfv->sd(ispin,ikp));
			  //tmap["ortho_align"].stop();
                    
			  //tmap["riccati"].start();
			  //wf.sd(ispin,ikp)->riccati(*wfv->sd(ispin,ikp));
			  //tmap["riccati"].stop();
                    
			  if (wf.ultrasoft()) {
			    //tmap["usfns"].start();
			    wf.sd(ispin,ikp)->update_usfns();
			    //tmap["usfns"].stop();
			  }
			  //ewd: lowdin doesn't yet work correctly with ultrasoft
			  //tmap["lowdin"].start();
			  //wf.sd(ispin,ikp)->lowdin();
			  //tmap["lowdin"].stop();
			  //tmap["gram"].start();
			  wf.sd(ispin,ikp)->gram();
			  //tmap["gram"].stop();
			}
		    }
		}
            }
	}
    }
}


