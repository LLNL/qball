////////////////////////////////////////////////////////////////////////////////
//
// jade.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: jade.C,v 1.4 2008-09-08 15:56:20 draeger Exp $

#include <cmath>
#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits> // epsilon
#include <iostream>
#include <iomanip>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef SCALAPACK
#include "blacs.h"
#endif

#include "Context.h"
#include "Matrix.h"
#include "blas.h"
#include "Timer.h"
using namespace std;

int jade(int maxsweep, double tol, vector<DoubleMatrix*> a,
  DoubleMatrix& u, vector<vector<double> >& adiag)
{
  Timer tm_comm;
  const bool debug_diag_sum = false;

  const double eps = numeric_limits<double>::epsilon();
  assert(tol>eps);
  const Context& ctxt = u.context();
  // The input matrices are *a[k]
  // the orthogonal transformation is returned in u
  // on exit, the matrices a[k] are maximally diagonal,
  // u contains the orthogonal transformation
  // adiag[k][i] contains the diagonal elements of the a[k]'s

  for ( int k = 0; k < a.size(); k++ )
  {
    assert(a[k]->context() == u.context());
    assert(a[k]->m()==a[k]->n());
    assert(a[k]->m()==u.n());
    assert(a[k]->m()==u.m());
    assert(a[k]->mb()==u.mb());
    assert(a[k]->nb()==u.nb());
  }

  const int mloc = a[0]->mloc();
  const int nloc = a[0]->nloc();
  //cout << ctxt.mype() << ": nloc: " << nloc << endl;

  // identify the last active process column
  // process columns beyond that column do not have any elements of a[k]
  // compute num_nblocks = total number of column blocks
  // if num_nblocks >= ctxt.npcol(), all process columns are active
  // otherwise, the last active process column has index num_nblocks-1
  const int num_nblocks = u.n() / u.nb() + ( u.n()%u.nb() == 0 ? 0 : 1 );
  const int last_active_process_col = min(ctxt.npcol()-1, num_nblocks-1);

  // initialize u with the identity
  u.identity();

  // eigenvalue array
  adiag.resize(a.size());
  for ( int k = 0; k < a.size(); k++ )
    adiag[k].resize(a[k]->n());

  // check if the local number of rows is odd
  const bool nloc_odd = ( a[0]->nloc()%2 != 0 );

  // if nloc is odd, auxiliary arrays are created to host an extra column
  // for both a[k] and u
  vector<vector<double> > a_aux(a.size());
  vector<double> u_aux;
  if ( nloc_odd )
  {
    for ( int k = 0; k < a.size(); k++ )
      a_aux[k].resize(mloc);
    u_aux.resize(mloc);
  }

  // compute local number of pairs nploc
  const int nploc = (a[0]->nloc()+1)/2;
  // dimension of top and bot arrays is nploc: local number of pairs
  deque<int> top(nploc), bot(nploc);

  // compute total number of pairs np
  int np = nploc;
  ctxt.isum('r',1,1,&np,1);
  //cout << ctxt.mype() << ": np=" << np << endl;
  // initialize top and bot arrays
  // the pair i is (top[i],bot[i])
  // top[i] is the local index of the top column of pair i
  // bot[i] is the local index of the bottom column of pair i
  for ( int i = 0; i < nploc; i++ )
    top[i] = i;
  for ( int i = 0; i < nploc; i++ )
    bot[nploc-i-1] = nploc+i;
  // if top[i] or bot[i] == nloc,the data resides in the array a_aux or u_aux

  // jglobal: global column index
  // jglobal[i] is the global column index of the column residing in
  // the local vector i. If nloc_odd and i==2*nploc-1, jglobal[i] == -1
  vector<int> jglobal(2*nploc,-1);
  for ( int jblock = 0; jblock < a[0]->nblocks(); jblock++ )
    for ( int y = 0; y < a[0]->nbs(jblock); y++ )
      jglobal[y + jblock*a[0]->nb()] = a[0]->j(jblock,y);

  // store addresses of columns of a and of u in acol and ucol
  vector<vector<double*> > acol(a.size());
  vector<double*> ucol(2*nploc);
  for ( int k = 0; k < a.size(); k++ )
  {
    acol[k].resize(2*nploc);
    for ( int i = 0; i < a[k]->nloc(); i++ )
      acol[k][i] = a[k]->valptr(i*a[k]->mloc());
    // if nloc is odd, store the address of vector 2*nploc-1
    if ( nloc_odd )
      acol[k][2*nploc-1] = &a_aux[k][0];
  }
  for ( int i = 0; i < u.nloc(); i++ )
    ucol[i] = u.valptr(i*u.mloc());
  // if nloc is odd, store the address of vector 2*nploc-1
  if ( nloc_odd )
    ucol[2*nploc-1] = &u_aux[0];

  //for ( int i = 0; i < acol.size(); i++ )
  //  cout << ctxt.mype() << ": acol[" << i << "]=" << acol[i] << endl;
  //for ( int i = 0; i < ucol.size(); i++ )
  //  cout << ctxt.mype() << ": ucol[" << i << "]=" << ucol[i] << endl;

  // the vectors of the pair (top[i],bot[i]) are located at
  // addresses acol[top[i]] and acol[bot[i]]

  bool done = false;
  int nsweep = 0;
  // allocate matrix element packed array apq
  // apq[3*ipair   + k*3*nploc] = apq[k][ipair]
  // apq[3*ipair+1 + k*3*nploc] = app[k][ipair]
  // apq[3*ipair+2 + k*3*nploc] = aqq[k][ipair]
  vector<double> apq(a.size()*3*nploc);

  double diag_sum = 0.0, previous_diag_sum = 0.0;
  while ( !done )
  {
    // sweep: process local pairs and rotate 2*np-1 times
    nsweep++;
    double diag_change = 0.0;
    for ( int irot = 0; irot < 2*np-1; irot++ )
    {
      //cout << ctxt.mype() << ": top[i]: ";
      //for ( int i = 0; i < nploc; i++ )
      //  cout << setw(3) << top[i];
      //cout << endl;

      //cout << ctxt.mype() << ": bot[i]: ";
      //for ( int i = 0; i < nploc; i++ )
      //  cout << setw(3) << bot[i];
      //cout << endl;

      //cout << ctxt.mype() << ": jglobal[top[i]]: ";
      //for ( int i = 0; i < nploc; i++ )
      //  cout << setw(3) << jglobal[top[i]];
      //cout << endl;

      //cout << ctxt.mype() << ": jglobal[bot[i]]: ";
      //for ( int i = 0; i < nploc; i++ )
      //  cout << setw(3) << jglobal[bot[i]];
      //cout << endl;

      // perform Jacobi rotations for all local pairs

      // compute off-diagonal matrix elements apq for all pairs
      // skip the pair if one or both of the vectors is a dummy vector
      // i.e. a vector having jglobal==-1

      int mloc = a[0]->mloc();
      for ( int k = 0; k < a.size(); k++ )
      {
        for ( int ipair = 0; ipair < nploc; ipair++ )
        {
          const int iapq = 3*ipair + k*3*nploc;
          //cout << ctxt.mype() << ": computing apq for global pair "
          //     << jglobal[top[ipair]] << " " << jglobal[bot[ipair]] << endl;
          apq[iapq]   = 0.0;
          apq[iapq+1] = 0.0;
          apq[iapq+2] = 0.0;
          if ( jglobal[top[ipair]] >= 0 && jglobal[bot[ipair]] >= 0 )
          {
            const double *ap = acol[k][top[ipair]];
            const double *aq = acol[k][bot[ipair]];
            const double *up = ucol[top[ipair]];
            const double *uq = ucol[bot[ipair]];
            int one = 1;
            apq[iapq]   = ddot(&mloc,ap,&one,uq,&one);
            apq[iapq+1] = ddot(&mloc,ap,&one,up,&one);
            apq[iapq+2] = ddot(&mloc,aq,&one,uq,&one);
          }
        }
      } // for k
      // apq now contains partial sums of matrix elements
      tm_comm.start();
      int len = apq.size();
      ctxt.dsum('c',len,1,&apq[0],len);
      tm_comm.stop();
      // apq now contains the matrix elements

      for ( int ipair = 0; ipair < nploc; ipair++ )
      {
        if ( jglobal[top[ipair]] >= 0 &&
             jglobal[bot[ipair]] >= 0 )
        {
          // compute rotation sine and cosine
          // Cardoso-Souloumiac expressions for the rotation angle

          // compute 2x2 matrix g
          double g11 = 0.0, g12 = 0.0, g22 = 0.0;

          for ( int k = 0; k < a.size(); k++ )
          {
            const int iapq = 3*ipair + k*3*nploc;
            const double app = apq[iapq+1];
            const double aqq = apq[iapq+2];

            const double h1 = app - aqq;
            const double h2 = 2.0 * apq[iapq];

            g11 += h1 * h1;
            g12 += h1 * h2;
            g22 += h2 * h2;
          }

          // the matrix g is real, symmetric
          // compute Jacobi rotation diagonalizing g

          double c = 1.0, s = 0.0, e1 = g11, e2 = g22;
          if ( g12*g12 > eps * eps * fabs(g11*g22) )
          {
            double tau = 0.5 * ( g22 - g11 ) / g12;
            double t = 1.0 / ( fabs(tau) + sqrt(1.0 + tau*tau));
            if ( tau < 0.0 ) t *= -1.0;
            c = 1.0 / sqrt(1.0 + t*t);
            s = t * c;

            // eigenvalues
            e1 -= t * g12;
            e2 += t * g12;
          }

          // components of eigenvector associated with the largest eigenvalue
          double x,y;
          if ( e1 > e2 )
          {
            x = c;
            y = -s;
          }
          else
          {
            x = s;
            y = c;
          }

          // choose eigenvector with x positive to ensure small angle
          if ( x < 0.0 )
          {
            x = -x;
            y = -y;
          }

          // compute Jacobi rotation R(p,q)
          c = sqrt(0.5*(x+1.0));
          s = y / sqrt(2.0*(x+1.0));

          // apply the rotation R(p,q)
          //
          //           |  c   s |
          //  R(p,q) = |        |
          //           | -s   c |

          // U := U * R(p,q)^T
          // A := A * R(p,q)^T

          // apply rotation to columns of a and u

          // the drot function computes
          // c*x + s*y -> x
          //-s*x + c*y -> y
          // call drot with args c, -s

          //cout << " p=" << jglobal[top[ipair]]
          //     << " q=" << jglobal[bot[ipair]]
          //     << " g11=" << g11 << " g22=" << g22
          //     << " g12=" << g12 << endl;
          //cout << " c=" << c << " s=" << s << endl;

          int one = 1;
          for ( int k = 0; k < a.size(); k++ )
          {
            double *ap = acol[k][top[ipair]];
            double *aq = acol[k][bot[ipair]];
            int mloc = a[k]->mloc();
            drot(&mloc,ap,&one,aq,&one,&c,&s);
          }
          double *up = ucol[top[ipair]];
          double *uq = ucol[bot[ipair]];
          int mloc = u.mloc();
          drot(&mloc,up,&one,uq,&one,&c,&s);

          // new value of off-diag element apq
          for ( int k = 0; k < a.size(); k++ )
          {
            const int iapq = 3*ipair + k*3*nploc;
            const double app = apq[iapq+1];
            const double aqq = apq[iapq+2];
            const double apqnew = (c*c-s*s)*apq[iapq] - c*s*(app-aqq);

            // accumulate change in sum of squares of diag elements
            // note negative sign: decrease in offdiag is increase in diag
            diag_change -= 2.0 * ( apqnew*apqnew - apq[iapq]*apq[iapq] );
          }
        }
      } // for ipair

      // all local pairs have been processed

      // rotate top and bot arrays
      if ( nploc > 0 )
      {
        bot.push_back(top.back());
        top.pop_back();
        top.push_front(bot.front());
        bot.pop_front();

        // make rotation skip element 0 on the first process column
        // if my process column is zero, swap top[0] and top[1]
        if ( ctxt.mycol() == 0 )
        {
          if ( nploc > 1 )
          {
            int tmp = top[0];
            top[0] = top[1];
            top[1] = tmp;
          }
          else
          {
            // if there is only one local pair, exchange top[0] and bot[0]
            int tmp = top[0];
            top[0] = bot[0];
            bot[0] = tmp;
          }
        }

        // exchange columns of a[k] and u

        int rbufi_left, rbufi_right, sbufi_left, sbufi_right;
        // send buffers contain k columns of a and one of u
        int bufsize = (a.size()+1)*a[0]->mloc();
        vector<double> sbuf_left(bufsize), sbuf_right(bufsize);
        vector<double> rbuf_left(bufsize), rbuf_right(bufsize);

        // on each task except mycol==npcol-1
        // send jglobal[bot[nploc-1]] to the right
        // if jglobal != -1 send vector bot[nploc-1] to the right

        // on each task except mycol==npcol-1
        // recv jglobal from the right
        // if jglobal != -1 recv a vector from the right into bot[nploc-1]
        // set value of jglobal[bot[nploc-1]]

        // on each task except mycol==0
        // send jglobal[top[0]] to the left
        // if jglobal != -1 send vector top[0] to the left

        // on each task except mycol==0
        // recv jglobal from the left
        // if jglobal != -1 recv a vector from the left into top[0]
        // set value of jglobal[top[0]]

        // exchange jglobal values first
        tm_comm.start();
        if ( ctxt.mycol() < last_active_process_col )
        {
          sbufi_right = jglobal[bot[nploc-1]];
          ctxt.isend(1,1,&sbufi_right,1,ctxt.myrow(),ctxt.mycol()+1);
          ctxt.irecv(1,1,&rbufi_right,1,ctxt.myrow(),ctxt.mycol()+1);
          jglobal[bot[nploc-1]] = rbufi_right;
          //cout << ctxt.mype() << ": received jglobal="
          //     << jglobal[bot[nploc-1]] << " from right" << endl;
        }
        if ( ctxt.mycol() != 0 )
        {
          sbufi_left = jglobal[top[0]];
          ctxt.isend(1,1,&sbufi_left,1,ctxt.myrow(),ctxt.mycol()-1);
          ctxt.irecv(1,1,&rbufi_left,1,ctxt.myrow(),ctxt.mycol()-1);
          jglobal[top[0]] = rbufi_left;
          //cout << ctxt.mype() << ": received jglobal="
          //     << jglobal[top[0]] << " from left" << endl;
        }

        // exchange column vectors

        //ewdif ( ctxt.mycol() < last_active_process_col )
        if ( ctxt.mycol() < last_active_process_col && bufsize > 0)
        {
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(&sbuf_right[k*mloc],acol[k][bot[nploc-1]],
                   mloc*sizeof(double));
          }
          memcpy(&sbuf_right[a.size()*mloc], ucol[bot[nploc-1]],
                 mloc*sizeof(double) );
          ctxt.dsend(bufsize,1,&sbuf_right[0],bufsize,
                     ctxt.myrow(),ctxt.mycol()+1);
          ctxt.drecv(bufsize,1,&rbuf_right[0],bufsize,
                     ctxt.myrow(),ctxt.mycol()+1);
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(acol[k][bot[nploc-1]],&rbuf_right[k*mloc],
                   mloc*sizeof(double));
          }
          memcpy(ucol[bot[nploc-1]], &rbuf_right[a.size()*mloc],
                 mloc*sizeof(double) );
          //cout << ctxt.mype() << ": received vector jglobal="
          //     << jglobal[bot[nploc-1]] << " from right" << endl;
        }
        //ewdif ( ctxt.mycol() != 0 )
        if ( ctxt.mycol() != 0 && bufsize > 0 )
        {
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(&sbuf_left[k*mloc],acol[k][top[0]],mloc*sizeof(double));
          }
          memcpy(&sbuf_left[a.size()*mloc],ucol[top[0]],mloc*sizeof(double) );
          ctxt.dsend(bufsize,1,&sbuf_left[0],bufsize,
                     ctxt.myrow(),ctxt.mycol()-1);
          ctxt.drecv(bufsize,1,&rbuf_left[0],bufsize,
                     ctxt.myrow(),ctxt.mycol()-1);
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(acol[k][top[0]],&rbuf_left[k*mloc],mloc*sizeof(double) );
          }
          memcpy(ucol[top[0]],&rbuf_left[a.size()*mloc],mloc*sizeof(double) );
          //cout << ctxt.mype() << ": received vector jglobal="
          //     << jglobal[top[0]] << " from left" << endl;
        }
        tm_comm.stop();
      } // if nploc > 0

      // end of step

    } // for irot

    // sweep is complete
    tm_comm.start();
    ctxt.dsum('r',1,1,&diag_change,1);
    tm_comm.stop();

    // compute sum of squares of diagonal elements

    if ( debug_diag_sum )
    {
      // compute sum of squares of diagonal elements using current values
      // (after rotation)
      previous_diag_sum = diag_sum;
      diag_sum = 0.0;
      for ( int k = 0; k < a.size(); k++ )
      {
        for ( int ipair = 0; ipair < nploc; ipair++ )
        {
          double tmp[2] = { 0.0, 0.0 };
          // compute the diagonal elements
          // skip dummy vectors
          int one = 1;
          int mloc = a[k]->mloc();
          if ( jglobal[top[ipair]] >= 0 )
          {
            const double *ap = acol[k][top[ipair]];
            const double *up = ucol[top[ipair]];
            tmp[0] = ddot(&mloc,ap,&one,up,&one);
          }
          if ( jglobal[bot[ipair]] >= 0 )
          {
            const double *aq = acol[k][bot[ipair]];
            const double *uq = ucol[bot[ipair]];
            tmp[1] = ddot(&mloc,aq,&one,uq,&one);
          }
          // tmp now contains partial sums of app and aqq
          ctxt.dsum('c',2,1,tmp,2);
          // tmp now contains the diagonal elements app and aqq
          diag_sum += tmp[0]*tmp[0] + tmp[1]*tmp[1];
        }
      }
      ctxt.dsum('r',1,1,&diag_sum,1);
      const double diag_sum_increase = diag_sum - previous_diag_sum;
      if ( ctxt.onpe0() )
        cout << " jade: nsweep=" << nsweep
             << " dsum: "
             << setw(15) << setprecision(10) << diag_sum
             << " dsum_inc: "
             << setw(15) << setprecision(10) << diag_sum_increase << endl;
    }

    if ( ctxt.onpe0() )
      cout << " jade: nsweep=" << nsweep
           << " dchange: "
           << setw(15) << setprecision(10) << diag_change << endl;

    done = ( ( fabs(diag_change) < tol ) || ( nsweep >= maxsweep ) );

  } // while !done

  // if a dummy vector was used, (i.e. if nloc_odd), the dummy vector
  // may end up anywhere in the array after all rotations are completed.
  // The array a_aux may contain a (non-dummy) vector.

  if ( nloc_odd )
  {
    // find position of the dummy vector and copy a_aux onto it
    int idum = 0;
    while ( jglobal[idum] != -1 && idum < 2*nploc ) idum++;
    //cout << ctxt.mype() << ": idum=" << idum << endl;
    if ( idum != 2*nploc-1 )
    {
      for ( int k = 0; k < a.size(); k++ )
      {
        memcpy(acol[k][idum],&a_aux[k][0],mloc*sizeof(double));
      }
      memcpy(ucol[idum],&u_aux[0],mloc*sizeof(double));
    }
  }

  // compute diagonal values
  for ( int k = 0; k < a.size(); k++ )
  {
    for ( int i = 0; i < a[k]->n(); i++ )
      adiag[k][i] = 0.0;
    for ( int jblock = 0; jblock < a[k]->nblocks(); jblock++ )
      for ( int y = 0; y < a[k]->nbs(jblock); y++ )
      {
        // j is the global column index
        int j = a[k]->j(jblock,y);
        int jjj = y + jblock*a[k]->nb();
        const double *ap = a[k]->valptr(jjj*a[k]->mloc());
        const double *up = u.valptr(jjj*u.mloc());
        int mloc = a[k]->mloc();
        int one = 1;
        adiag[k][j] = ddot(&mloc,ap,&one,up,&one);
      }
    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    ctxt.dsum(a[k]->n(),1,&adiag[k][0],a[k]->n());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }

  if ( ctxt.onpe0() )
    cout << " jade: comm time: " << tm_comm.real() << endl;

  return nsweep;
}
