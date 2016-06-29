#include <complex>
#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main() {
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {

    cout << setprecision(10);
  
    // create a complex<double> array
    int zsize = 3;
    vector<complex<double> > z(zsize);
    z[0] = complex<double>(0.111111111,-0.444222444);
    z[1] = complex<double>(3.53,2.42);
    z[2] = complex<double>(0.8877889,9.876543);

    for (int i=0; i<zsize; i++) 
      cout << "z[" << i << "] = " << z[i] << endl;

    double* dp = (double*)&z[0];

    for (int i=0; i<2*zsize; i++) 
      cout << "dp[" << i << "] = " << dp[i] << endl;


    const double dtmp = 1.23456789;
    for (int i=0; i<2*zsize; i++)
      dp[i] += dtmp;
    cout << "Added constant " << dtmp << " to dp" << endl;

    for (int i=0; i<zsize; i++) 
      cout << "z[" << i << "] = " << z[i] << endl;



  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}  
  
  
  
