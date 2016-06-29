////////////////////////////////////////////////////////////////////////////////
//
// testWavefunction.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testWavefunction.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#include "Context.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Timer.h"

#include <iostream>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    // use: 
    // testWavefunction a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel nempty nspin nkp
    assert(argc==15);
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    cout << " volume: " << cell.volume() << endl;
    double ecut = atof(argv[10]);
    int nel = atoi(argv[11]);
    int nempty = atoi(argv[12]);
    int nspin = atoi(argv[13]);
    int nkp = atoi(argv[14]);
    
    Context ctxt;
    
    Wavefunction wf(ctxt);
    Timer tm;
    
    tm.reset(); tm.start();
    wf.resize(cell,cell,ecut);
    tm.stop();
    cout << " wf.resize: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset(); tm.start();
    wf.set_nel(nel);
    tm.stop();
    cout << " wf.set_nel: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
         
    tm.reset(); tm.start();
    wf.set_nspin(nspin);
    tm.stop();
    cout << " wf.set_nspin: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
    
    for ( int ikp = 0; ikp < nkp-1; ikp++ )
    {
      wf.add_kpoint(D3vector((0.5*(ikp+1))/(nkp-1),0,0),1.0);
    }
    
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        if ( wf.sd(ispin,ikp) != 0 && wf.sdcontext(ispin,ikp)->active() )
        {
          cout << "wf.sd(ispin=" << ispin << ",ikp=" << ikp << "): "
               << wf.sd(ispin,ikp)->c().m() << "x"
               << wf.sd(ispin,ikp)->c().n() << endl;
          cout << ctxt.mype() << ":"
               << " sdcontext[" << ispin << "][" << ikp << "]: "
               << wf.sd(ispin,ikp)->context();
        }
      }
    }

    tm.reset();
    tm.start();
    wf.randomize(0.1);
    tm.stop();
    cout << " wf.randomize: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset();
    tm.start();
    wf.gram();
    cout << " wf.gram: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
           
    cout << " copy constructor...";
    Wavefunction wfm(wf);
    cout << "done" << endl;
    wfm.gram();
    
#if 0
    wf.randomize(0.1);
    wf.update_occ(false);
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      if ( wf.sd[ikp] != 0 )
        cout << " ekin[" << ikp << "]: " << wf.sd[ikp]->ekin() << endl;
    }
    

#endif
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
