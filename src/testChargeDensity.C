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
#include "Sample.h"
#include "Timer.h"

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
