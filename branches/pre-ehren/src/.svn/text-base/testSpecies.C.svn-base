//
// testSpecies.C
//
// use: testSpecies uri
//

#include "Species.h"
#include "Context.h"
#include "SpeciesReader.h"
#include <iostream>
#include <cassert>
#include <string>
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
  
  Context ctxt;
  if ( argc != 2 )
  {
    cerr << "use: testSpecies uri" << endl;
    return 1;
  }
  
  Species s(ctxt,"unknown_name");
  
  SpeciesReader rdr(ctxt);
  
  string uri(argv[1]);
  
  try
  {
    cout << " s.uri() = " << s.uri() << endl;
    cout << " testSpecies: invoking readSpecies: uri=" << uri << endl;
    rdr.readSpecies(s,uri);
  }
  catch ( const SpeciesReaderException& e )
  {
    cout << " SpeciesReaderException caught in testSpecies" << endl;
    throw;
  }
  cerr << " SpeciesReader::readSpecies done" << endl;
  
  const double rcps = 1.0;
  
  try
  {
    s.initialize(rcps);
  }
  catch ( SpeciesInitException& e )
  {
    cout << " Exception in Species initialization: " << e.msg << endl;
    throw;
  }
  cerr << s;
  
#if 1
  
  if ( ctxt.oncoutpe() )
  {
  
  double dr = 0.01;
  double dg = 0.02;
  double v,dv;

  int n = 1000;
  
  for ( int l = 0; l <= s.lmax(); l++ )
  {
    cout << n << " Vps(l=" << l << ",r) " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.vpsr(l,r,v);
      cout << r << " " << v << endl;
    }

    cout << n << " dVps(l=" << l << ",r)/dr " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.dvpsr(l,r,v,dv);
      cout << r << " " << dv << endl;
    }
  }
  
  cout << n << " Vloc(g) " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.vlocg(g,v);
    cout << g << " " << v << endl;
  }
  
  cout << n << " dVloc(g)/dg " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.dvlocg(g,v,dv);
    cout << g << " " << dv << endl;
  }
  
  for ( int l = 0; l <= s.lmax(); l++ )
  {
    if ( l != s.llocal() )
    {
      cout << n << " Vnl(l=" << l << ",g) " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.vnlg(l,g,v);
        cout << g << " " << v << endl;
      }
   
      cout << n << " dVnl(l=" << l << ",g)/dg " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.dvnlg(l,g,v,dv);
        cout << g << " " << dv << endl;
      }
    }
  }
  
  } // oncoutpe
#endif

  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
