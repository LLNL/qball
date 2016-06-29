////////////////////////////////////////////////////////////////////////////////
//
// ParOptCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ParOptCmd.C,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#include <fstream>
#include<iostream>
#include<ctime>
#include<cassert>
#include "ParOptCmd.h"
#include "ParallelOptimizer.h"
#include "Timer.h"
using namespace std;


int ParOptCmd::action(int argc, char **argv) {

  if ( argc < 3 || argc > 5) {
    if ( ui->oncoutpe() )
      cout << " use: paropt filename niter_ionic" << endl;
      cout << " use: paropt filename niter_ionic niter_scf" << endl;
      cout << " use: paropt filename niter_ionic niter_scf niter_nonscf" << endl;
    return 1;
  }
  
  if ( s->wf.nst() == 0 ) {
    if ( ui->oncoutpe() )
      cout << " <!-- ParOptCmd: no states, cannot run -->" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 ) {
    if ( ui->oncoutpe() )
      cout << " <!-- ParOptCmd: ecut = 0.0, cannot run -->" << endl;
    return 1;
  }
  
  char* parcmdfile = argv[1];
  
  int niter = atoi(argv[2]);
  int nite = 1;
  int nitscf = 1;
  if ( argc == 4 ) {
    nitscf = atoi(argv[3]);
  }
  else if ( argc == 5 ) {
    nitscf = atoi(argv[3]);
    nite = atoi(argv[4]);
  }

  ParallelOptimizer* popt = new ParallelOptimizer(*s);
  
  if (ui->oncoutpe() )
    cout << "<parallel_optimize niter_ionic=\"" << niter << "\" niter_scf=\"" << nitscf << "\" niter_nonscf=\"" << nite << "\">" << endl;


  popt->optimize(niter,nitscf,nite);

  if (ui->oncoutpe() )
    cout << "</parallel_optimize>" << endl;

  // write out input commands to parcmdfile
  if (ui->oncoutpe() ) {
    ofstream os;
    os.open(parcmdfile,ofstream::out);
    os << "set nrowmax " << s->wf.nrowmax() << endl;
    if (s->wf.nkp() > 1) 
      os << "set nparallelkpts " << s->wf.nparallelkpts() << endl;
    if (s->ctrl.reshape_context)
      os << "set reshape_context ON" << endl;    
    os.close();
  }

  delete popt;
  return 0;
}
