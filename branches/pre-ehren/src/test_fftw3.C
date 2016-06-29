////////////////////////////////////////////////////////////////////////////////
//
// test_fftw.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Timer.h"

#include <iostream>
#include <complex>
#include <valarray>
using namespace std;
#include <cassert>

#include "fftw3.h"
extern "C" {
  void zdscal_(int *,double *,complex<double> *,int *);
}

int main(int argc, char **argv)
{
  const int niter = 10;
  const int np2_ = atoi(argv[1]);
  const int nvec_ = atoi(argv[2]);
  const int ldz = np2_ + 1;

  fftw_plan fwplan2, bwplan2;

  // resize array zvec holding columns
  valarray<complex<double> > zvec_(nvec_ * ldz);
  fftw_complex *in = (fftw_complex*) &zvec_[0];
  fftw_complex *out = in;
  
  // initialization of FFT libs

  int rank = 1;
  const int n[1] = { np2_ };
  const int inembed[1] = { ldz };
  const int onembed[1] = { ldz };
  int howmany = nvec_;
  int stride = 1;
  int dist = ldz;
#if FFTWMEASURE
  // FFTWMEASURE
  fwplan2 = fftw_plan_many_dft(rank,n,howmany,in,NULL,stride,dist,
                               out,NULL,stride,dist,FFTW_FORWARD,FFTW_MEASURE);
  bwplan2 = fftw_plan_many_dft(rank,n,howmany,in,NULL,stride,dist,
                               out,NULL,stride,dist,FFTW_BACKWARD,FFTW_MEASURE);
#else
  // FFTW_ESTIMATE
  fwplan2 = fftw_plan_many_dft(rank,n,howmany,in,NULL,stride,dist,
                               out,NULL,stride,dist,FFTW_FORWARD,FFTW_ESTIMATE);
  bwplan2 = fftw_plan_many_dft(rank,n,howmany,in,NULL,stride,dist,
                               out,NULL,stride,dist,FFTW_BACKWARD,FFTW_ESTIMATE);
#endif

  Timer t_fwd,t_bwd;

  for ( int iter = 0; iter < niter; iter++ )
  {
    t_bwd.start();
    fftw_execute(bwplan2);
    t_bwd.stop();
 

    t_fwd.start();
    fftw_execute(fwplan2);
    t_fwd.stop();

    int len = zvec_.size();
    int inc1 = 1;
    double fac = 3.14;
    zdscal_(&len,&fac,&zvec_[0],&inc1);
  }

  fftw_destroy_plan(fwplan2);
  fftw_destroy_plan(bwplan2);

  cout << " fwd: " << t_fwd.real()/niter << endl;
  cout << " fwd: time per transform (microseconds): " 
       << 1.e6*t_fwd.real()/(niter*nvec_)
       << endl;
  cout << " bwd: " << t_bwd.real()/niter << endl;
  cout << " bwd: time per transform (microseconds): " 
       << 1.e6*t_bwd.real()/(niter*nvec_)
       << endl;

  return 0;
}
