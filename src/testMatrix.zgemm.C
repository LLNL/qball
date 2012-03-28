// $Id: testMatrix.old2.C,v 1.1 2005/12/05 23:11:16 draeger1 Exp $
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
#include "blas.h"
#include "omp.h"
#ifdef BGQ
#include <bgpm/include/bgpm.h>

extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif
using namespace std;

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
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
const double maxrand = 0.000001;  // maximum random perturbation to identity matrix

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

   char* infilename = argv[1];
   ifstream infile(infilename);
   assert(argc == 2 || argc == 3);
   bool tcheck = false;
   bool tortho = false;
   bool teigen = false;
   if ( argc == 3 )
   {
      if ( !strcmp(argv[2],"-check") )
         tcheck = true;
      else if ( !strcmp(argv[2],"-ortho") )
         tortho = true;
      else if ( !strcmp(argv[2],"-eigen") )
         teigen = true;
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

      ComplexMatrix a(ctxt,m_a,n_a,mb_a,nb_a);
      ComplexMatrix b(ctxt,m_b,n_b,mb_b,nb_b);
      ComplexMatrix c(ctxt,m_c,n_c,mb_c,nb_c);

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
         cout << "a:  mloc = " << a.mloc() << ", nloc = " << a.nloc() << endl;
         cout << "b:  mloc = " << b.mloc() << ", nloc = " << b.nloc() << endl;
         cout << "c:  mloc = " << c.mloc() << ", nloc = " << c.nloc() << endl;
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

                  double drand1 =  rand()*nrandinv*maxrand;
                  double drand2 =  rand()*nrandinv*maxrand;
                  double aij;
                  if (i == j)
                     aij = 1.0 + drand1;
                  else
                     aij = drand1;

                  int iii = x + l*a.mb();
                  int jjj = y + m*a.nb();
                  int ival = iii + jjj * a.mloc();
                  a[ival] = complex<double>(aij,drand2);
               }

      for ( int m = 0; m < b.nblocks(); m++ )
         for ( int l = 0; l < b.mblocks(); l++ )
            for ( int y = 0; y < b.nbs(m); y++ )  
               for ( int x = 0; x < b.mbs(l); x++ )
               {
                  int i = b.i(l,x);
                  int j = b.j(m,y);
                  double drand1 =  rand()*nrandinv*maxrand;
                  double drand2 =  rand()*nrandinv*maxrand;
                  double bij;
                  if (i == j)
                     bij = 1.0 + drand1;
                  else
                     bij = drand1;
                  // double bij = b.i(l,x) * 10 + b.j(m,y);
                  //double bij = bb(i,j);
                  int iii = x + l*b.mb();
                  int jjj = y + m*b.nb();
                  int ival = iii + jjj * b.mloc();
                  b[ival] = complex<double>(bij,drand2);
               }

      if (tortho) {
#if 0      

         if (mype == 0) cout << "Running gemm..." << endl;
         tmap["gemm"].start();
         HPM_Start("gemm2");
         c.gemm(ta,tb,1.0,a,b,0.0);
         HPM_Stop("gemm2");
         tmap["gemm"].stop();

         // Gram-Schmidt orthogonalization from SlaterDet.C
         if ( mype == 0 ) cout << "Starting Gram-Schmidt orthogonalization..." << endl;
         tmap["gram"].start();

         DoubleMatrix a_proxy(a);
         DoubleMatrix s(ctxt,m_c,n_c,mb_c,nb_c);

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

         DoubleMatrix ssq(ctxtsq,m_c,n_c,mbsq,nbsq);
         tmap["getsub1"].start();
         ssq.getsub(s,s.m(),s.n(),0,0);
         tmap["getsub1"].stop();

         tmap["Cholesky"].start();
         ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
         tmap["Cholesky"].stop();

         tmap["getsub2"].start();
         s.getsub(ssq,ssq.m(),ssq.n(),0,0);
         tmap["getsub2"].stop();

         // solve triangular system X * L^T = C
         tmap["trsm"].start();
         a_proxy.trsm('r','l','t','n',1.0,s);
         tmap["trsm"].stop();
         tmap["gram"].stop();

         if (mype==0) cout << "Finished with Gram-Schmidt." << endl;
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
         assert(mb_c == nb_c);
         valarray<double> w(s.m());

         //trace_start();

         tmap["syev"].start();
         s.syev('l',w,c);
         tmap["syev"].stop();

         //trace_stop();
#endif
      }
      else {

         if (npes == 1) {
            complex<double> zzero = complex<double>(0.0,0.0);
            complex<double> zone = complex<double>(1.0,0.0);
            char cc='c';
            char cn='n';

            int kk = m_a;
            int mm = m_c;
            int nn = n_c;
            vector<complex<double> > avec(mm*kk);
            vector<complex<double> > bvec(nn*kk);
            vector<complex<double> > cvec(mm*nn);
            for (int ii=0; ii<avec.size(); ii++) {
               double drand1 =  rand()*nrandinv*maxrand;
               double drand2 =  rand()*nrandinv*maxrand;
               avec[ii] = complex<double>(drand1,drand2);
            }
            for (int ii=0; ii<bvec.size(); ii++) {
               double drand1 =  rand()*nrandinv*maxrand;
               double drand2 =  rand()*nrandinv*maxrand;
               bvec[ii] = complex<double>(drand1,drand2);
            }
            for (int ii=0; ii<cvec.size(); ii++) {
               double drand1 =  rand()*nrandinv*maxrand;
               double drand2 =  rand()*nrandinv*maxrand;
               cvec[ii] = complex<double>(drand1,drand2);
            }
               
            tm.start();
            HPM_Start("zgemm1");
            zgemm(&cc,&cn,&mm,&nn,&kk,&zone,&avec[0],&kk,&bvec[0],&kk,&zzero,&cvec[0],&mm);
            HPM_Stop("zgemm1");
            tm.stop();

         }
         else {
            tm.start();
            HPM_Start("gemm1");
            //c.gemm(ta,tb,1.0,a,b,0.0);
            c.gemm('c','n',1.0,a,b,0.0);
            HPM_Stop("gemm1");
            tm.stop();
         }
#if 0         
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
#endif      
         int nthreads = omp_get_max_threads();
         if (mype == 0)
            cout << "M = " << m_c << " N = " << n_c << " K = " << m_a << " zgemm time = " << setprecision(5) << setw(8) << tm.real()<< " sec, GFlops = " << (8.0e-9*m_c*n_c*m_a) / tm.real() << " on " << npes << " pes, " << nthreads << " threads" << endl;
         //cout << "M = " << m_c << " N = " << n_c << " K = " << m_a << " LDA = " << m_a << " LDB = " << m_a << " LDC = " << m_c << " zgemm time = " << setprecision(5) << setw(8) << tm.real()<< " sec, GFlops = " << (6.0e-9*m_c*n_c*m_a) / tm.real() << endl;
             

         /*
         cout << " mype = " << mype << " CPU/Real: " << setw(8) << tm.cpu() 
              << " / " << setw(8) << tm.real();
         if ( tm.real() > 0.0 )
         {
            int kmax = ( ta == 'n' ) ? a.n() : a.m();
            cout << "  GFlops: " 
                 << (6.0e-9*m_c*n_c*kmax) / tm.real() << endl;
         }
         */
#if 0    
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
#endif
         
      }
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
}
