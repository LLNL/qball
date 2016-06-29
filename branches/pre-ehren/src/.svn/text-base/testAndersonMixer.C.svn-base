////////////////////////////////////////////////////////////////////////////////
//
// testAndersonMixer.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testAndersonMixer.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $


#include <iostream>
using namespace std;

#include "Context.h"
#include "AndersonMixer.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    Context ctxt;
 
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << ctxt.mype() << " on " << processor_name << endl;

    vector<double> x,f,xlast,fbar;
    double theta;
    x.resize(3);
    xlast.resize(3);
    f.resize(3);
    fbar.resize(3);
    
    AndersonMixer mixer(x.size(),ctxt);
    
    x[0] = 1.0 + 0.2 * ctxt.mype();
    x[1] = 2.0 + 0.2 * ctxt.mype();
    x[2] = 3.0 + 0.2 * ctxt.mype();
    
    for ( int iter = 0; iter < 100; iter++ )
    {
      f[0] = - 2.0 * x[0] - 3.0 * x[0]*x[0]*x[0];
      f[1] = - 1.0 * x[1] - 1.5 * x[1]*x[1]*x[1];
      f[2] = - 4.0 * x[2];
      
      double resnorm = 0.0;
      for ( int i = 0; i < x.size(); i++ )
      {
        resnorm += f[i]*f[i];
        f[i] *= 0.1;
      }
      ctxt.dsum(1,1,&resnorm,1);

      cout << ctxt.mype() << ": x: " << x[0] << " " << x[1] << " " << x[2] 
           << " resnorm: " << resnorm << endl;
    
      mixer.update(&f[0],&theta,&fbar[0]);
      
      for ( int i = 0; i < x.size(); i++ )
      {
        // xbar = x + theta * ( x - xlast )
        // x = xbar + fbar
        double xtmp = x[i];
        x[i] += theta * ( x[i] - xlast[i] ) + fbar[i];
        xlast[i] = xtmp;
      }
    }
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
