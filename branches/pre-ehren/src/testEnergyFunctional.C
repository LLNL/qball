////////////////////////////////////////////////////////////////////////////////
//
// testEnergyFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testEnergyFunctional.C,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#include "Context.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
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
    // use: testEnergyFunctional a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel
    assert(argc==12);
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    double ecut = atof(argv[10]);
    int nel = atoi(argv[11]);
    
    Timer tm;
    
    Context ctxt;
    cout << " initial context: " << ctxt;
    Sample s(ctxt);
    s.atoms.addAtom(
      new Atom("G","gaussium",D3vector(0.2,0.3,0),D3vector(0,0,0)));
    s.atoms.listAtoms();
    s.atoms.listSpecies();
    
    s.wf.resize(cell,cell,ecut);
    s.wf.set_nel(nel);
    
    if ( ctxt.oncoutpe() ) cout << " nel: " << s.wf.nel() << endl;

    s.wf.update_occ();
    //s.wf.randomize(0.05);
    
    tm.reset();
    tm.start();
    s.wf.gram();
    tm.stop();
    cout << " Gram: CPU/Real: " << tm.cpu() << " / " << tm.real() << endl;

    tm.reset();
    tm.start();
    EnergyFunctional ef(s);
    tm.stop();
    cout << " EnergyFunctional:ctor: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
    
    tm.reset();
    tm.start();
    cout << " ef.energy(): " << ef.energy() << endl;
    tm.stop();
    cout << " EnergyFunctional:energy: CPU/Real: " 
         << tm.cpu() << " / " << tm.real() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
