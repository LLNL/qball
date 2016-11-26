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
#include "ExtrapolatorASP.h"

using namespace std;

ExtrapolatorASP::ExtrapolatorASP(){
}

void ExtrapolatorASP::extrapolate_wavefunction(Wavefunction & wf, Wavefunction* wfv, Wavefunction* wfmm, int iter, double dt, const Context& ctxt){

  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      if (wf.spinactive(ispin))
	{
	  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
            {
              if (wf.kptactive(ikp))
		{
		  assert(wf.sd(ispin,ikp) != 0);
		      if ( ctxt.mype()==0 )
			cout << "Extrapolating wavefunction using ASP algorithm." << endl;
                   
		      double* c = (double*) wf.sd(ispin,ikp)->c().cvalptr();
		      double* cv = (double*) wfv->sd(ispin,ikp)->c().cvalptr();
		      double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();
		      const int mloc = wf.sd(ispin,ikp)->c().mloc();
		      const int nloc = wf.sd(ispin,ikp)->c().nloc();
		      const int len = 2*mloc*nloc;
		      if ( iter == 0 )
			{
			  for ( int i = 0; i < len; i++ )
			    {
			      const double x = c[i];
			      const double v = cv[i];
			      // extrapolation using velocity in cv
			      c[i] = x + dt * v;
			      cv[i] = x;
			    }
			  //ewd: 10-5-12b, uncomment this
			  if (wf.ultrasoft()) {
			    //tmap["usfns"].start();
			    wf.sd(ispin,ikp)->update_usfns();
			    //tmap["usfns"].stop();
			  }
			  //tmap["gram"].start();
			  wf.sd(ispin,ikp)->gram();
			  //tmap["gram"].stop();
			}
		      else if ( iter == 1 )
			{
			  //wfv->align(wf);
			  for ( int i = 0; i < len; i++ )
			    {
			      const double x = c[i];
			      const double xm = cv[i];
			      c[i] = 2.0 * x - xm;
			      cv[i] = x;
			      cmm[i] = xm;
			    }
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
			  // align wf with wfmm before extrapolation
			  // wf.align(*wfmm);
			  // wfmm->align(wf);
                  
			  // extrapolate
			  for ( int i = 0; i < len; i++ )
			    {
			      const double x = c[i];   // current wf (scf converged) at t
			      const double xm = cv[i]; // extrapolated wf at t
			      const double xmm = cmm[i]; // extrapolated wf at t-dt
			      const double asp_a1 = 0.5;
			      c[i] = 2.0 * x - xm +
				asp_a1 * ( x - 2.0 * xm + xmm );
			      //c[i] = 2.5 * x - 2.0 * xm + 0.5 * xmm;
			      cmm[i] = xm;
			      cv[i] = x;
			    }
                    
			  // orthogonalize the extrapolated value
			  if (wf.ultrasoft()) {
			    //tmap["usfns"].start();
			    wf.sd(ispin,ikp)->update_usfns();
			    //tmap["usfns"].stop();
			  }
			  //tmap["gram"].start();
			  wf.sd(ispin,ikp)->gram();
			  //tmap["gram"].stop();
			  //tmap["lowdin"].start();
			  //wf.sd(ispin,ikp)->lowdin();
			  //tmap["lowdin"].stop();
                  
			  // c[i] now contains the extrapolated value
			}
		      // c[i] is now ready for electronic iterations
		}
	    }
	}
    }
}



