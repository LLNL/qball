////////////////////////////////////////////////////////////////////////////////
//
//
// Test of Intel SSE2 C++ Class library
//
// Compute the real part of a complex scalar product
//
// Compile on Xeon with: icc -O3 -xW thisfile.C
//
// Run with: a.out 100
//
////////////////////////////////////////////////////////////////////////////////


#include <dvec.h>
#include <complex>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
long long readTSC(void)
{
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
  return ticks.complete;
}

////////////////////////////////////////////////////////////////////////////////
double f(int n, complex<double> *a, complex<double> *b)
{
  double sum = 0.0;
  for ( int i = 0; i < n; i++ )
  {
    sum += a[i].real()*b[i].real() + a[i].imag()*b[i].imag();
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double ff(int n, complex<double> *a, complex<double> *b)
{
  F64vec2 *aa = (F64vec2*) a;
  F64vec2 *bb = (F64vec2*) b;
  F64vec2 cc(0.0); 
  #pragma unroll(4)
  for ( int i = 0; i < n; i++ )
  {
    cc += aa[i] * bb[i];
  }
  return add_horizontal(cc);
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  int n = atoi(argv[1]);
  
  // allocate 16-byte aligned arrays
  complex<double> *a = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *b = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *aa = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *bb = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  
  // initialize arrays
  for ( int i = 0; i < n; i++ )
  {
    a[i] = (double) i;
    b[i] = complex<double>(0.5*i,i);
    aa[i] = a[i];
    bb[i] = b[i];
  }
  
  long long t0;
  long long tf[10],tff[10];
  
  double sumf,sumff;
  for ( int i = 0; i < 10; i++ )
  {
    t0 = readTSC();
    sumf = f(n,a,b);
    tf[i] = readTSC()-t0;
    
    t0 = readTSC();
    sumff = ff(n,aa,bb);
    tff[i] = readTSC()-t0;
  }
  
  // check correctness
  assert(sumf==sumff);
    
  cout << " cycles per point:" << endl;
  for ( int i = 0; i < 10; i++ )
    cout << "f: "     << ((double)tf[i])/n 
         << "   ff: " << ((double)tff[i])/n << endl;

  return 0;
}
