////////////////////////////////////////////////////////////////////////////////
//
// testSample.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testSample.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#include <iostream>
using namespace std;

#include "Context.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "D3vector.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  {
    Context ctxt;
 
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << ctxt.mype() << " on " << processor_name << endl;
 
    Sample s(ctxt);
    
    D3vector cell(18,18,18);
    double ecut = 25.0;
    s.wf.resize(cell,cell,ecut);
    s.wf.set_nel(12*54);
    
    s.wf.randomize(1.e-4);
    s.wf.gram();
    cout << " ortho_error: " << s.wf.sd[0][0]->ortho_error() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
