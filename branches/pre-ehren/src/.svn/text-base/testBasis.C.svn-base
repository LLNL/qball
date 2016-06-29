////////////////////////////////////////////////////////////////////////////////
//
// testBasis.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testBasis.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#include "Basis.h"
#include "Context.h"
#include <iostream>
#include <new>
#include <cstdlib>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
  // use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z ecut kx ky kz npr npc
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    if ( argc !=16 )
    {
      cout <<
      " use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z ecut kx ky kz npr npc"
      << endl;
      return 1;
    }
    const D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    const D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    const D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    
    double ecut = atof(argv[10]);
    D3vector kpoint(atof(argv[11]),atof(argv[12]),atof(argv[13]));
    int npr = atoi(argv[14]);
    int npc = atoi(argv[15]);
    
    Context ctxt(npr,npc);
    UnitCell cell(a0,a1,a2);
    
    Basis basis(ctxt,kpoint);
    try
    {
      basis.resize(cell,cell,ecut);
    }
    catch ( bad_alloc )
    {
      cout << " bad_alloc caught in Basis::resize" << endl;
      throw;
    }
    
    //cout << basis;
    
    //Basis b2(basis);
    //cout << b2;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
