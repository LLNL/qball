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
// testChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <qball/Context.h>
#include "Wavefunction.h"
#include "ChargeDensity.h"
#include "SlaterDet.h"
#include <qball/FourierTransform.h>
#include "Sample.h"
#include <qball/Timer.h>

#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

int main(int argc, char **argv)
{
   int mype;
#if USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);  
#endif
  {
    // use: 
    // testChargeDensity a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel nempty nspin nkp
    assert(argc==15);
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    if (mype == 0)
       cout << " volume: " << cell.volume() << endl;
    double ecut = atof(argv[10]);
    int nel = atoi(argv[11]);
    int nempty = atoi(argv[12]);
    int nspin = atoi(argv[13]);
    int nkp = atoi(argv[14]);
    
    Timer tm;
    
    Context ctxt;
    Sample* s = new Sample(ctxt);

    // need to set default values for any s->ctrl parameters used by test code
    s->ctrl.ecutden = 0.0;
    s->ctrl.ultrasoft = false;
    s->ctrl.nlcc = false;
    s->ctrl.tddft_involved = false;
    s->ctrl.extra_memory = 3;

    s->wf.set_cell(cell);
    s->wf.set_ecut(ecut);    
    //s->wf.resize(cell,cell,ecut);
    s->wf.set_nel(nel);
    s->wf.set_nspin(nspin);

    /*
    for ( int ikp = 0; ikp < nkp-1; ikp++ )
    {
      s->wf.add_kpoint(D3vector((0.5*(ikp+1))/(nkp-1),0,0),1.0);
    }
    */
    
    /*
    for ( int ispin = 0; ispin < s->wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < s->wf.nkp(); ikp++ )
      {
        if ( s->wf.sd(ispin,ikp) != 0 && s->wf.sdcontext(ispin,ikp)->active() )
        {
        cout << "s->wf.sd(ispin=" << ispin << ",ikp=" << ikp << "): "
             << s->wf.sd(ispin,ikp)->c().m() << "x"
             << s->wf.sd(ispin,ikp)->c().n() << endl;
        cout << ctxt.mype() << ":"
             << " sdcontext[" << ispin << "][" << ikp << "]: " 
             << s->wf.sd(ispin,ikp)->context();
        }
      }
    }
    */
    
    s->wf.randomize(0.1,false);

    tm.reset();
    tm.start();
    s->wf.gram();
    if (mype == 0)
       cout << " s->wf.gram: CPU/Real: " 
            << tm.cpu() << " / " << tm.real() << endl;
         
    s->wf.update_occ(0.0,0);
        
    // compute charge density in real space
    Timer tmrho;
    tmrho.reset();
    tmrho.start();
    ChargeDensity cd(*s);
    tmrho.stop();
    if (mype == 0)
       cout << " ChargeDensity::constructor: CPU/Real: " 
            << tmrho.cpu() << " / " << tmrho.real() << endl;
         
    tmrho.reset();
    tmrho.start();
#ifdef HPM  
  HPM_Start("update_density");
#endif
    cd.update_density();
#ifdef HPM  
  HPM_Stop("update_density");
#endif
    tmrho.stop();
    if (mype == 0)
       cout << " ChargeDensity::update_density: CPU/Real: " 
            << tmrho.cpu() << " / " << tmrho.real() << endl;
         
    // print the first few Fourier coefficients of the charge density
    /*
    for ( int ispin = 0; ispin < s->wf.nspin(); ispin++ )
    {
      for ( int i = 0; i < cd.vbasis()->localsize(); i++ )
      {
        cout << ctxt.mype() << ": rho(ispin=" << ispin << ",ig=" << i << " (" 
        << cd.vbasis()->idx(3*i) << " "
        << cd.vbasis()->idx(3*i+1) << " " << cd.vbasis()->idx(3*i+2) << ") "
        << cd.rhog[ispin][i] << endl;
      }
    }
    */

    cd.print_timing();
    
#if 0
    // integral of rho in r space
    double sum = 0.0;
    for ( int i = 0; i < rho.size(); i++ )
      sum += rho[i];
#if USE_MPI
    double tsum;
    int mycol = sd.context().mycol();
    MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,sd.context().comm());
    sum = tsum;
#endif
    cout << ctxt.mype() << ": " << " rho: " << sum / ft.np012() << endl;

#endif

  }
#if USE_MPI
  MPI_Finalize();
#endif
}
