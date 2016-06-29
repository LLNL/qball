////////////////////////////////////////////////////////////////////////////////
//
// testSlaterDet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testSlaterDet.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#include "Context.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Timer.h"

#include <iostream>
#include <fstream>
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
    // use: testSlaterDet a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nst kx ky kz npr npc
    double err;
    assert(argc==17);
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    double ecut = atof(argv[10]);
    int nst = atoi(argv[11]);
    D3vector kpoint(atof(argv[12]),atof(argv[13]),atof(argv[14]));
    
    int npr = atoi(argv[15]);
    int npc = atoi(argv[16]);
    
    Timer tm;
    
    Context ctxt(npr,npc,'c');
    
    SlaterDet sd(ctxt,kpoint);
    cout << sd.context();
    
    sd.resize(cell,cell,ecut,nst);
    if ( ctxt.myproc() == 0 )
    {
      cout << " kpoint: " << sd.basis().kpoint() << endl;
      cout << " ngw:    " << sd.basis().size() << endl;
      cout << " nst:    " << sd.nst() << endl;
    }
    
    vector<double> occ(nst);
    for ( int i = 0; i < nst; i++ )
      occ[i] = 1.0;
    
    err = sd.ortho_error();
    if ( ctxt.myproc() == 0 )
      cout << " Initial orthogonality error before rand " << err << endl;
      
    sd.randomize(0.2/sd.basis().size());
    if ( ctxt.myproc() == 0 )
      cout << " sd.randomize: done" << endl;
    err = sd.ortho_error();
    if ( ctxt.myproc() == 0 )
      cout << " Orthogonality error after rand " << err << endl;
     
    tm.reset();
    tm.start();
    sd.gram();
    tm.stop();
    cout << " Gram: CPU/Real: " << tm.cpu() << " / " << tm.real() << endl;
    err = sd.ortho_error();
    if ( ctxt.myproc() == 0 )
      cout << " Gram orthogonality error " << err << endl;

    SlaterDet sdm(ctxt,kpoint);
    sdm.resize(cell,cell,ecut,nst);
         
    sdm.c() = sd.c();
    sd.randomize(0.1/sd.basis().size());
    
    tm.reset();
    tm.start();
    sd.riccati(sdm);
    tm.stop();
    cout << " riccati: CPU/Real: " << tm.cpu() << " / " << tm.real() << endl;
    err = sd.ortho_error();
    if ( ctxt.myproc() == 0 )
      cout << " Riccati orthogonality error " << err << endl;
    
    sd.riccati(sdm);
    err = sd.ortho_error();
    if ( ctxt.myproc() == 0 )
    {
      cout << " Riccati orthogonality error " << err << endl;
      cout << " sd.total size (MB): " << setprecision(4)
           << sd.memsize() / 1048576.0 << endl;
      cout << " sd.local size (MB): " << setprecision(4)
           << sd.localmemsize() / 1048576.0 << endl;
    }
    
    cout << " ekin: " << setprecision(8) << sd.ekin(occ) << endl;
        
    // compute charge density in real space
    FourierTransform ft(sd.basis(),
      2*sd.basis().np(0), 2*sd.basis().np(1), 2*sd.basis().np(2));
    
    vector<complex<double> > f(ft.np012loc());
    vector<double> rho(ft.np012loc());
    
    Timer tmrho;
    tmrho.reset();
    tmrho.start();
    
    cout << " compute_density..." << endl;
    double weight = 1.0;
    sd.compute_density(ft,occ,weight,&rho[0]);
    
    tmrho.stop();
    cout << " compute_density: CPU/Real: " 
         << tmrho.cpu() << " / " << tmrho.real() << endl;
         
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
    cout << ctxt.mype() << ": " 
         << " rho: " << sum * cell.volume() / ft.np012() << endl;
 
    SlaterDet sdp(sd);
    //SlaterDet sdp(ctxt,kpoint);
    //sdp.resize(cell,cell,ecut,nst);
    //sdp = sd;
    sdp.randomize(0.001);
    sdp.gram();
    err = sdp.ortho_error();
    if ( ctxt.myproc() == 0 )
      cout << " Gram orthogonality error " << err << endl;
      
    Timer tmv;
    tmv.reset();
    tmv.start();
    sd.rs_mul_add(ft,&rho[0],sdp);
    tmv.stop();
    cout << " rs_mul_add: CPU/Real: " 
         << tmv.cpu() << " / " << tmv.real() << endl;
         
#if 0
    ofstream outfile("sd.dat");
    tm.reset();
    tm.start();
    outfile << sd;
    tm.stop();
    cout << " write: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
#endif
      
    for ( TimerMap::iterator i = sd.tmap.begin(); i != sd.tmap.end(); i++ )
      cout << ctxt.mype() << ": "
           << setw(15) << (*i).first 
           << " : " << setprecision(3) << (*i).second.real() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
