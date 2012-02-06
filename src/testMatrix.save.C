// $Id: testMatrix.C,v 1.3 2006/01/10 01:15:38 draeger1 Exp $
//
// test Matrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testMatrix input_file [-check] [-ortho]
// input_file:
// nprow npcol
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <map>
using namespace std;

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_APC
#include "apc.h"
#endif

#include "Context.h"
#include "Matrix.h"

// mpi_trace functions
//extern "C" {
//  extern  void    trace_start();
//  extern  void    trace_stop();
//}
//

double aa(int i, int j) { return 1.0/(i+1)+2.0/(j+1); }
double bb(int i, int j) { return i-j-3; }

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.0001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{

  // set up map of timers
  map<string,Timer> tmap;
  tmap["total"].start();

  // choose random number seed based on system clock
  srand((unsigned)time(NULL));

  int mype;
  int npes;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#else
  npes=1;
  mype=0;
#endif

#if USE_APC
  ApcInit();
#endif

  char* infilename = argv[1];
  ifstream infile(infilename);

  assert(argc == 2 || argc == 3);
  bool tcheck = false;
  bool tortho = false;
  bool tortholoop = false;
  bool tchol = false;
  bool teigencomplex = false;
  bool teigen = false;
  bool teigensweep = false;
  if ( argc == 3 )
  {
    if ( !strcmp(argv[2],"-check") )
      tcheck = true;
    else if ( !strcmp(argv[2],"-ortho") )
      tortho = true;
    else if ( !strcmp(argv[2],"-ortholoop") )
      tortholoop = true;
    else if ( !strcmp(argv[2],"-cholesky") )
      tchol = true;
    else if ( !strcmp(argv[2],"-eigencomplex") )
      teigencomplex = true;
    else if ( !strcmp(argv[2],"-eigen") )
      teigen = true;
    else if ( !strcmp(argv[2],"-eigensweep") )
      teigensweep = true;
    else
    {
      cerr << " invalid argv[2]" << endl;
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
  }
  Timer tm;
  int nprow, npcol;
  int m_a, n_a, mb_a, nb_a;
  int m_b, n_b, mb_b, nb_b;
  int m_c, n_c, mb_c, nb_c;
  char ta, tb;
  if(mype == 0)
  {
    infile >> nprow >> npcol;
    cout<<"nprow="<<nprow<<", npcol="<<npcol<<endl;
    infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
    infile >> m_b >> n_b >> mb_b >> nb_b >> tb;
    cout<<"m_b="<<m_b<<", n_b="<<n_a<<endl;
    infile >> m_c >> n_c >> mb_c >> nb_c;
    cout<<"m_c="<<m_c<<", n_c="<<n_c<<endl;
  }
#ifdef USE_MPI
  MPI_Bcast(&nprow, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&npcol, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&ta, 1, MPI_CHAR, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&tb, 1, MPI_CHAR, 0, MPI_COMM_WORLD);    
#endif
  {  
    if ( ta == 'N' ) ta = 'n';
    if ( tb == 'N' ) tb = 'n';

    Context ctxt(nprow,npcol);

    if ( mype == 0 )
    {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }

    DoubleMatrix a(ctxt,m_a,n_a,mb_a,nb_a);
    DoubleMatrix b(ctxt,m_b,n_b,mb_b,nb_b);
    DoubleMatrix c(ctxt,m_c,n_c,mb_c,nb_c);

    if ( mype == 0 )
    {
      cout << " m_a x n_a / mb_a x nb_a / ta = "
           << a.m() << "x" << a.n() << " / "
           << a.mb() << "x" << a.nb() << " / " << ta << endl;
      cout << " m_b x n_b / mb_b x nb_b / tb = "
           << b.m() << "x" << b.n() << " / "
           << b.mb() << "x" << b.nb() << " / " << tb << endl;
      cout << " m_c x n_c / mb_c x nb_c      = "
           << c.m() << "x" << c.n() << " / "
           << c.mb() << "x" << c.nb() << endl;
    }

    for ( int m = 0; m < a.nblocks(); m++ )
      for ( int l = 0; l < a.mblocks(); l++ )
        for ( int y = 0; y < a.nbs(m); y++ )  
          for ( int x = 0; x < a.mbs(l); x++ )
          {
            int i = a.i(l,x);
            int j = a.j(m,y);
            // double aij = a.i(l,x) * 10 + a.j(m,y);
            //double aij = aa(i,j);

	    double drand =  rand()*nrandinv*maxrand;
	    double aij;
	    if (i == j)
	      aij = 1.0 + drand;
	    else
	      aij = drand;

            int iii = x + l*a.mb();
            int jjj = y + m*a.nb();
            int ival = iii + jjj * a.mloc();
            a[ival] = aij;
          }

    for ( int m = 0; m < b.nblocks(); m++ )
      for ( int l = 0; l < b.mblocks(); l++ )
        for ( int y = 0; y < b.nbs(m); y++ )  
          for ( int x = 0; x < b.mbs(l); x++ )
          {
            int i = b.i(l,x);
            int j = b.j(m,y);
            // double bij = b.i(l,x) * 10 + b.j(m,y);
            double bij = bb(i,j);
            int iii = x + l*b.mb();
            int jjj = y + m*b.nb();
            int ival = iii + jjj * b.mloc();
            b[ival] = bij;
          }

    if (tortho) {

#ifdef USE_APC
  ApcStart(1);
#endif
      if (mype == 0) cout << "Running gemm..." << endl;
      tmap["gemm"].start();
      c.gemm(ta,tb,1.0,a,b,0.0);
      tmap["gemm"].stop();

    // Gram-Schmidt orthogonalization from SlaterDet.C
      //if ( mype == 0 ) cout << "Starting Gram-Schmidt orthogonalization..." << endl;
      tmap["gram"].start();

      DoubleMatrix a_proxy(a);
      DoubleMatrix s(ctxt,m_c,n_c,mb_c,nb_c);


      if ( mype == 0 ) cout << "Starting syrk..." << endl;
      tmap["syrk"].start();
      //s.syrk('l','t',2.0,a_proxy,0.0);
      s.syrk('l','t',2.0,a,0.0);
      tmap["syrk"].stop();

      //tmap["syr"].start();
      //s.syr('l',-1.0,a_proxy,0,'r');
      //tmap["syr"].stop();

      double dnpes = (double) npes;
      int npsq = (int) sqrt(dnpes);
      if (mype == 0) cout << "Reshaping context from " << nprow << " x " << npcol << " to " << npsq << " x " << npsq << endl;
      Context ctxtsq(npsq,npsq);

      // check both cases where m_c is/isn't multiple of npsq
      assert(m_c == n_c);
      int mbsq, nbsq;
      if (m_c%npsq == 0) {
	mbsq = m_c/npsq;
	nbsq = mbsq;
      }
      else {
	mbsq = (int) m_c/npsq + 1;
	nbsq = (int) n_c/npsq + 1;
      }
      if (mype == 0) cout << "mbsq, nbsq = " << mbsq << ", " << nbsq << endl;

      if ( mype == 0 ) cout << "Starting potrf..." << endl;
      tmap["Cholesky_orig"].start();
      s.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_orig"].stop();

      DoubleMatrix ssq(ctxtsq,m_c,n_c,mbsq,nbsq);
      tmap["getsub1"].start();
      ssq.getsub(s,s.m(),s.n(),0,0);
      tmap["getsub1"].stop();
      
      tmap["Cholesky_remap"].start();
      ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_remap"].stop();

      tmap["getsub2"].start();
      s.getsub(ssq,ssq.m(),ssq.n(),0,0);
      tmap["getsub2"].stop();

      // solve triangular system X * L^T = C
      if ( mype == 0 ) cout << "Starting trsm..." << endl;
      tmap["trsm"].start();
      a_proxy.trsm('r','l','t','n',1.0,s);
      tmap["trsm"].stop();
      tmap["gram"].stop();

      if (mype==0) cout << "Finished with Gram-Schmidt." << endl;

#ifdef USE_APC
  ApcStop(1);
#endif

    }
    if (tortholoop) {

#ifdef USE_APC
  ApcStart(1);
#endif

    // Gram-Schmidt orthogonalization from SlaterDet.C
      //if ( mype == 0 ) cout << "Starting Gram-Schmidt orthogonalization..." << endl;

    for (int i=0; i<20; i++) {
    
      DoubleMatrix a_proxy(a);
      DoubleMatrix s(ctxt,m_c,n_c,mb_c,nb_c);

      if ( mype == 0 ) cout << "Starting syrk..." << endl;
      //s.syrk('l','t',2.0,a_proxy,0.0);
      s.syrk('l','t',2.0,a,0.0);

      //s.syr('l',-1.0,a_proxy,0,'r');

      //if ( mype == 0 ) cout << "Starting potrf..." << endl;
      //s.potrf('l'); // Cholesky decomposition: S = L * L^T

      // solve triangular system X * L^T = C
      if ( mype == 0 ) cout << "Starting trsm..." << endl;
      a_proxy.trsm('r','l','t','n',1.0,s);

      if (mype==0) cout << "Finished iteration " << i << endl;

    }

#ifdef USE_APC
  ApcStop(1);
#endif

    }
    else if (tchol) { 

      int npsq1 = npcol;
      int npsq2 = (int) sqrt((double) npes); // largest possible square context
      assert(npcol <= npsq2);
      int npsq3 = ((int) npsq2/npcol)*npcol; // nearest integer multiple of npcol
      assert(npcol <= npsq3);

      if (mype == 0) cout << "npsq2 = " << npsq2 << ", npsq3 = " << npsq3 << endl;

      int mbsq1 = m_c/npsq1 + (m_c%npsq1 == 0 ? 0 : 1);
      int mbsq2 = m_c/npsq2 + (m_c%npsq2 == 0 ? 0 : 1);
      int mbsq3 = m_c/npsq3 + (m_c%npsq3 == 0 ? 0 : 1);

      if (mype == 0) cout << "m_c = " << m_c << ", mbsq1 = " << mbsq1 << ", mbsq2 = " << mbsq2 << ", mbsq3 = " << mbsq3 << endl;

      // square matrix on original context
      DoubleMatrix s(ctxt,m_c,m_c,mb_c,mb_c);

      // square matrix on npcol by npcol context
      Context ctxtsq1(npsq1,npsq1);
      DoubleMatrix ssq1(ctxtsq1,m_c,m_c,mbsq1,mbsq1);

      // square matrix on largest possible square context
      Context ctxtsq2(npsq2,npsq2);
      DoubleMatrix ssq2(ctxtsq2,m_c,m_c,mbsq2,mbsq2);

      // square matrix on square context which is largest integer multiple of npcol
      Context ctxtsq3(npsq3,npsq3);
      DoubleMatrix ssq3(ctxtsq3,m_c,m_c,mbsq3,mbsq3);

      if (mype == 0) { 
	cout << "orig matrix s:  " << s.m() << " x " << s.n() << ", " << s.nbs(0) << " x " << s.mbs(0) << " blocks " << endl;
	cout << "matrix ssq1, on ncol x ncol context:  " << ssq1.m() << " x " << ssq1.n() << ", " << ssq1.nbs(0) << " x " << ssq1.mbs(0) << " blocks " << endl;
	cout << "matrix ssq2, on ncol x ncol context:  " << ssq2.m() << " x " << ssq2.n() << ", " << ssq2.nbs(0) << " x " << ssq2.mbs(0) << " blocks " << endl;
	cout << "matrix ssq3, on ncol x ncol context:  " << ssq3.m() << " x " << ssq3.n() << ", " << ssq3.nbs(0) << " x " << ssq3.mbs(0) << " blocks " << endl;
      }

      // fill all matrices in the same way
      for ( int m = 0; m < s.nblocks(); m++ )
	for ( int l = 0; l < s.mblocks(); l++ )
	  for ( int y = 0; y < s.nbs(m); y++ )  
	    for ( int x = 0; x < s.mbs(l); x++ )
	      {
		int i = s.i(l,x);
		int j = s.j(m,y);
		double drand =  rand()*nrandinv*maxrand;
		double sij;
		if (i == j)
		  sij = 1.0*i + drand;
		else
		  sij = drand;

		int iii = x + l*s.mb();
		int jjj = y + m*s.nb();
		int ival = iii + jjj * s.mloc();
		s[ival] = sij;
	      }
      for ( int m = 0; m < ssq1.nblocks(); m++ )
	for ( int l = 0; l < ssq1.mblocks(); l++ )
	  for ( int y = 0; y < ssq1.nbs(m); y++ )  
	    for ( int x = 0; x < ssq1.mbs(l); x++ )
	      {
		int i = ssq1.i(l,x);
		int j = ssq1.j(m,y);
		double drand =  rand()*nrandinv*maxrand;
		double sij;
		if (i == j)
		  sij = 1.0*i + drand;
		else
		  sij = drand;

		int iii = x + l*ssq1.mb();
		int jjj = y + m*ssq1.nb();
		int ival = iii + jjj * ssq1.mloc();
		ssq1[ival] = sij;
	      }
      for ( int m = 0; m < ssq2.nblocks(); m++ )
	for ( int l = 0; l < ssq2.mblocks(); l++ )
	  for ( int y = 0; y < ssq2.nbs(m); y++ )  
	    for ( int x = 0; x < ssq2.mbs(l); x++ )
	      {
		int i = ssq2.i(l,x);
		int j = ssq2.j(m,y);
		double drand =  rand()*nrandinv*maxrand;
		double sij;
		if (i == j)
		  sij = 1.0*i + drand;
		else
		  sij = drand;

		int iii = x + l*ssq2.mb();
		int jjj = y + m*ssq2.nb();
		int ival = iii + jjj * ssq2.mloc();
		ssq2[ival] = sij;
	      }
      for ( int m = 0; m < ssq3.nblocks(); m++ )
	for ( int l = 0; l < ssq3.mblocks(); l++ )
	  for ( int y = 0; y < ssq3.nbs(m); y++ )  
	    for ( int x = 0; x < ssq3.mbs(l); x++ )
	      {
		int i = ssq3.i(l,x);
		int j = ssq3.j(m,y);
		double drand =  rand()*nrandinv*maxrand;
		double sij;
		if (i == j)
		  sij = 1.0*i + drand;
		else
		  sij = drand;

		int iii = x + l*ssq3.mb();
		int jjj = y + m*ssq3.nb();
		int ival = iii + jjj * ssq3.mloc();
		ssq3[ival] = sij;
	      }

      // now time Cholesky times for all four matrices
      // note:  if s, on original context, is only distributed on 
      // square subset of processors, time should be very similar to ssq1

      tmap["Cholesky_orig"].start();
      s.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_orig"].stop();
      tmap["Cholesky_sq2"].start();
      ssq2.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_sq2"].stop();
      tmap["Cholesky_sq1"].start();
      ssq1.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_sq1"].stop();
      tmap["Cholesky_sq3"].start();
      ssq3.potrf('l'); // Cholesky decomposition: S = L * L^T
      tmap["Cholesky_sq3"].stop();

    }
    else if (teigencomplex) { 

      ComplexMatrix s(ctxt,m_c,n_c,mb_c,nb_c);

      for ( int m = 0; m < s.nblocks(); m++ )
	for ( int l = 0; l < s.mblocks(); l++ )
	  for ( int y = 0; y < s.nbs(m); y++ )  
	    for ( int x = 0; x < s.mbs(l); x++ ) {
              int i = s.i(l,x);
              int j = s.j(m,y);

              double drand1 =  rand()*nrandinv*maxrand;
              double drand2 =  rand()*nrandinv*maxrand;

              double sij;
              if (i == j)
                sij = 1.0*i + drand1;
              else
                sij = drand1;

              double zij;
              if (i == j)
                zij = 1.0*i + drand2;
              else 
                zij = drand2;
              
              int iii = x + l*s.mb();
              int jjj = y + m*s.nb();
              int ival = iii + jjj * s.mloc();
              complex<double> cij = (sij,zij);
              s[ival] = cij;
            }


      assert(m_c == n_c);
      valarray<double> w(s.m());

      if (mype == 0) 
	cout << "Calculating complex eigenvalues on rectangular context..." << endl;

      ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());

      tmap["heev-rect"].start();
      //s.heev('l',w);
      //s.heev('l',w,z);
      s.heevd('l',w,z);
      tmap["heev-rect"].stop();

      if (mype == 0) {
        cout << "Eigenvalues:  " << endl;
        for (int i=0; i<s.m(); i++) {
          cout << "  " << w[i];
          if (i%8 == 0) cout << endl;
        }
        cout << endl;
      }

      if (mype == 0) 
	cout << "Done." << endl;

    }
    else if (teigen) { 

      DoubleMatrix s(ctxt,m_c,n_c,mb_c,nb_c);

      for ( int m = 0; m < s.nblocks(); m++ )
	for ( int l = 0; l < s.mblocks(); l++ )
	  for ( int y = 0; y < s.nbs(m); y++ )  
	    for ( int x = 0; x < s.mbs(l); x++ )
	      {
		int i = s.i(l,x);
		int j = s.j(m,y);

		double drand =  rand()*nrandinv*maxrand;
		//double drand = 0.0;

		double sij;
		if (i == j)
		  sij = 1.0*i + drand;
		else
		  sij = drand;

		int iii = x + l*s.mb();
		int jjj = y + m*s.nb();
		int ival = iii + jjj * s.mloc();
		s[ival] = sij;
	      }

      //tmap["syrk"].start();
      //s.syrk('l','t',2.0,a,0.0);
      //tmap["syrk"].stop();

      assert(m_c == n_c);
      //assert(mb_c == nb_c);
      valarray<double> w(s.m());

      //trace_start();

      if (mype == 0) 
	cout << "Calculating eigenvalues on rectangular context..." << endl;

      tmap["syevd-rect"].start();
      s.syevd('l',w,c);
      tmap["syevd-rect"].stop();

      if (mype == 0) 
	cout << "Mapping matrix onto square context..." << endl;
      // remap onto square-ish context
      //int tmpcol = npes;
      //int tmprow = 1;
      //while (tmpcol > tmprow) {
      //tmprow *= 2;
      //tmpcol /= 2;
      //}
      int tmpcol = 1;
      int tmprow = npes;
      while (tmpcol < tmprow) {
	tmprow /= 2;
	tmpcol *= 2;
      }

      if (mype == 0) cout << "Reshaping context from " << nprow << " x " << npcol << " to " << tmprow << " x " << tmpcol << endl;
      Context ctxtsq(tmprow,tmpcol);
      int mbsq = m_c/tmprow + (m_c%tmprow == 0 ? 0 : 1);
      int nbsq = n_c/tmpcol + (n_c%tmpcol == 0 ? 0 : 1);

      if (mbsq > nbsq) 
	nbsq = mbsq;
      else
	mbsq = nbsq;

      if (mype == 0) cout << "New context:  " << tmprow << " x " << tmpcol << ", subblocks = " << mbsq << " x " << nbsq << endl;

      DoubleMatrix ssq(ctxtsq,m_c,n_c,mbsq,nbsq);
      DoubleMatrix csq(ctxtsq,m_c,n_c,mbsq,nbsq);

      tmap["eigen-getsub1"].start();
      ssq.getsub(s,s.m(),s.n(),0,0);
      tmap["eigen-getsub1"].stop();

      valarray<double> wsq(ssq.m());

      if (mype == 0) 
	cout << "Calculating eigenvalues on square context..." << endl;
      tmap["syevd-sq"].start();
      ssq.syevd('l',wsq,csq);
      tmap["syevd-sq"].stop();

      if (mype == 0) 
	cout << "Mapping matrix back onto rectangular context..." << endl;
      tmap["eigen-getsub2"].start();
      s.getsub(ssq,ssq.m(),ssq.n(),0,0);
      tmap["eigen-getsub2"].stop();

      if (mype == 0) 
	cout << "Done." << endl;
      //trace_stop();



    }
    else if (teigensweep) { 

      string timernames[17];
      timernames[0] = "1row";
      timernames[1] = "2row";
      timernames[2] = "4row";
      timernames[3] = "8row";
      timernames[4] = "16row";
      timernames[5] = "32row";
      timernames[6] = "64row";
      timernames[7] = "128row";
      timernames[8] = "256row";
      timernames[9] = "512row";
      timernames[10] = "1024row";
      timernames[11] = "2048row";
      timernames[12] = "4096row";
      timernames[13] = "8192row";
      timernames[14] = "16384row";
      timernames[15] = "32768row";
      timernames[16] = "65536row";

      assert(npes%2==0);

      int ctxtcnt = 7;
      for (int nrow=128;nrow>=8;nrow/=2) {


	//	int ncol = npes/nrow;
	int ncol = nrow;

	//int nrow = nprow;
	//int ncol = npcol;

	if (mype == 0) cout << "Creating context: " << nrow << " x " << ncol << endl;

	Context ctxtsw(nrow,ncol);
	assert (m_c == n_c);
	int mbsw = m_c/nrow + (m_c%nrow == 0 ? 0 : 1);
	int nbsw = n_c/ncol + (n_c%ncol == 0 ? 0 : 1);

	if (mbsw > nbsw) 
	  nbsw = mbsw;
	else 
	  mbsw = nbsw;

	DoubleMatrix s(ctxtsw,m_c,n_c,mbsw,nbsw);
	DoubleMatrix csw(ctxtsw,m_c,n_c,mbsw,nbsw);

	for ( int m = 0; m < s.nblocks(); m++ )
	  for ( int l = 0; l < s.mblocks(); l++ )
	    for ( int y = 0; y < s.nbs(m); y++ )  
	      for ( int x = 0; x < s.mbs(l); x++ )
		{
		  int i = s.i(l,x);
		  int j = s.j(m,y);
		  
		  double drand =  rand()*nrandinv*maxrand;
		  //double drand = 0.0;
		  
		  double sij;
		  if (i == j)
		    sij = 1.0*i + drand;
		  else
		    sij = drand;
		  
		  int iii = x + l*s.mb();
		  int jjj = y + m*s.nb();
		  int ival = iii + jjj * s.mloc();
		  s[ival] = sij;
		}

	valarray<double> w(s.m());

	//trace_start();
	string tname;

	if (mype == 0) 
	  cout << "Calculating eigenvalues with syevd..." << endl;

	tname = timernames[ctxtcnt] + "-syevd";
	if (mype == 0) 
	  cout << "timer name = " << tname << endl;

	tmap[tname].start();
	s.syevd('l',w,csw);
	tmap[tname].stop();

	double time,tmin,tmax;
	time = tmap[tname].cpu();
	tmin = time;
	tmax = time;
    
	ctxt.dmin(1,1,&tmin,1);
	ctxt.dmax(1,1,&tmax,1);
	if (mype == 0) { 
	  cout << "  timing "
	       << setw(10) << tname
	       << " : " << setprecision(4) << setw(9) << tmin
	       << " "   << setprecision(4) << setw(9) << tmax << endl;
	}



	if (mype == 0) 
	  cout << "Calculating eigenvalues with syev..." << endl;

	tname = timernames[ctxtcnt] + "-syev";
	if (mype == 0) 
	  cout << "timer name = " << tname << endl;

	tmap[tname].start();
	s.syev('l',w,csw);
	tmap[tname].stop();

	time = tmap[tname].cpu();
	tmin = time;
	tmax = time;
    
	ctxt.dmin(1,1,&tmin,1);
	ctxt.dmax(1,1,&tmax,1);
	if (mype == 0) { 
	  cout << "  timing "
	       << setw(10) << tname
	       << " : " << setprecision(4) << setw(9) << tmin
	       << " "   << setprecision(4) << setw(9) << tmax << endl;
	}

	if (mype == 0) 
	  cout << "Done." << endl;

	//trace_stop();
	
	ctxtcnt--;

      }

    }
    else {

      tm.start();
      c.gemm(ta,tb,1.0,a,b,0.0);
      tm.stop();

      if ( tcheck )
	{
	  cout << " checking results..." << endl;
	  for ( int m = 0; m < c.nblocks(); m++ )
	    for ( int l = 0; l < c.mblocks(); l++ )
	      for ( int y = 0; y < c.nbs(m); y++ )  
		for ( int x = 0; x < c.mbs(l); x++ )
		  {
		    int i = c.i(l,x);
		    int j = c.j(m,y);
		    double sum = 0.0;
		    int kmax = ( ta == 'n' ) ? a.n() : a.m();
		    
		    if ( ( ta == 'n' ) && ( tb == 'n' ) )
		      {
			for ( int k = 0; k < kmax; k++ )
			  sum += aa(i,k) * bb(k,j);
		      }
		    else if ( ( ta != 'n' ) && ( tb == 'n' ) )
		      {
			for ( int k = 0; k < kmax; k++ )
			  sum += aa(k,i) * bb(k,j);
		      }
		    else if ( ( ta == 'n' ) && ( tb != 'n' ) )
		      {
			for ( int k = 0; k < kmax; k++ )
			  sum += aa(i,k) * bb(j,k);
		      }
		    else if ( ( ta != 'n' ) && ( tb != 'n' ) )
		      {
			for ( int k = 0; k < kmax; k++ )
			  sum += aa(k,i) * bb(j,k);
		      }
		    
		    int iii = x + l*c.mb();
		    int jjj = y + m*c.nb();
		    int ival = iii + jjj * c.mloc();
		    if ( fabs( c[ival] - sum ) > 1.e-8 )
		      {
			cout << " error at element (" << i << "," << j << ") "
			     << c[ival] << " " << sum << endl;
			exit(1);
		      }
		  }
          
          
	  cout << " results checked" << endl;
	}
      
      cout << " CPU/Real: " << setw(8) << tm.cpu() 
	   << " / " << setw(8) << tm.real();
      if ( tm.real() > 0.0 )
	{
	  int kmax = ( ta == 'n' ) ? a.n() : a.m();
	  cout << "  MFlops: " 
	       << (2.0e-6*m_c*n_c*kmax) / tm.real() << endl;
	}
#if 1    
      double norma=a.nrm2();
      if(mype == 0)cout<<"Norm(a)="<<norma<<endl;
      if(mype == 0)cout<<"DoubleMatrix::matgather..."<<endl;
      double*  aa=new double[a.m()*a.n()];
      a.matgather(aa, a.m());
      if(mype == 0)cout<<"DoubleMatrix::init..."<<endl;
      b.init(aa, a.m());
      double norm=b.nrm2();
      if ( mype == 0 ) cout << "Norm(b)=" << norm << endl;
      if ( fabs(norm-norma)>0.000001 )
	cout << "DoubleMatrix: problem with matgather/init" << endl;
      
      if ( c.n() == b.m() && c.m() == b.n() )
	{
	  if(mype == 0)cout<<"DoubleMatrix::transpose..."<<endl;
	  c.transpose(1.0,b,0.0);
	  norm=c.nrm2();
	  if(mype == 0)cout<<"Norm(c)="<<norm<<endl;
	}
      
      if(mype == 0)cout<<"DoubleMatrix::scal..."<<endl;
      c.scal(0.5);
    
      if ( a.m() == b.m() && a.n() == b.n() )
	{
	  if(mype == 0)cout<<"DoubleMatrix::axpy..."<<endl;
	  a.axpy(-2., b);
	}
      
      if ( a.m() == c.m() && a.n() == c.n() )
	{
	  if(mype == 0)cout<<"DoubleMatrix::operator=..."<<endl;
	  c=a;
	}
    
      if(mype == 0)cout<<"DoubleMatrix::nrm2..."<<endl;
      norm=c.nrm2();
      if (mype == 0) cout<<"Norm="<<norm<<endl;
      
      a.identity();
      DoubleMatrix a2(a);
      a -= a2;
      norm = a.nrm2();
      if (mype == 0) cout << "Norm(a)=" << norm << endl;
      
      // Eigenvalues and eigenvectors of c if c is square
      if ( c.m() == c.n() && c.mb() == c.nb() ) {
	if (mype == 0) cout << "Eigenproblem... ";
	DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
	valarray<double> w(c.m());
	c.syev('l',w,z);
	if (mype == 0) cout << " done" << endl;
      }
    }
#endif
    tmap["total"].stop();

    for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
      double time = (*i).second.cpu();
      double tmin = time;
      double tmax = time;
    
      ctxt.dmin(1,1,&tmin,1);
      ctxt.dmax(1,1,&tmax,1);
      if (mype == 0) { 
	
	cout << "  timing "
	     << setw(10) << (*i).first
	     << " : " << setprecision(4) << setw(9) << tmin
	     << " "   << setprecision(4) << setw(9) << tmax << endl;
      }
    }

  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
#if USE_APC
  ApcFinalize();
#endif
}
