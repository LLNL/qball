////////////////////////////////////////////////////////////////////////////////
//
// testChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testChargeDensity.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#include "Context.h"
#include "Wavefunction.h"
#include "ChargeDensity.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Timer.h"

#include <iostream>
#include <iomanip>
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
    // testChargeDensity a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel nempty nspin nkp
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
    
    Timer tm;
    
    Context ctxt;
    Wavefunction wf(ctxt);
    
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
         
    wf.update_occ();
        
    // compute charge density in real space
    Timer tmrho;
    tmrho.reset();
    tmrho.start();
    cout << " ChargeDensity::ctor..." << endl;
    ChargeDensity cd(wf);
    
    cout << "done" << endl;
    
    tmrho.stop();
    cout << " ChargeDensity::ctor: CPU/Real: " 
         << tmrho.cpu() << " / " << tmrho.real() << endl;
         
    tmrho.reset();
    tmrho.start();
    cout << " ChargeDensity::update_density..." << endl;
    cd.update_density();
    tmrho.stop();
    cout << " ChargeDensity::update_density: CPU/Real: " 
         << tmrho.cpu() << " / " << tmrho.real() << endl;
         
    // print the first few Fourier coefficients of the charge density
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int i = 0; i < cd.vbasis()->localsize(); i++ )
      {
        cout << ctxt.mype() << ": rho(ispin=" << ispin << ",ig=" << i << " (" 
        << cd.vbasis()->idx(3*i) << " "
        << cd.vbasis()->idx(3*i+1) << " " << cd.vbasis()->idx(3*i+2) << ") "
        << cd.rhog[ispin][i] << endl;
      }
    }
         
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
