////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// jacobi.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <cmath>
#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_SCALAPACK
#include "blacs.h"
#endif

#include "Context.h"
#include <math/Matrix.h>
#include <math/blas.h>

int jacobi(int maxsweep, double threshold, DoubleMatrix& a, DoubleMatrix& u,
              vector<double>& e)
{
  assert(threshold>1.e-14);
  const Context& ctxt = a.context();
  // The input matrix is a
  // the orthogonal transformation is returned in u
  // on exit, a is diagonal, u contains eigenvectors, e contains eigenvalues
  assert(a.m()==a.n());
  assert(u.m()==u.n());
  assert(a.m()==u.m());
  assert(a.mb()==u.mb());
  assert(a.nb()==u.nb());

  const int mloc = a.mloc();
  const int nloc = a.nloc();
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
  e.resize(a.n());

  // check if the local number of rows is odd
  const bool nloc_odd = ( a.nloc()%2 != 0 );

  // if nloc is odd, an auxiliary array is created to host an extra column
  // for both a and u
  vector<double> a_aux, u_aux;
  if ( nloc_odd )
  {
    a_aux.resize(mloc);
    u_aux.resize(mloc);
  }

  // compute local number of pairs nploc
  const int nploc = (a.nloc()+1)/2;
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
  for ( int jblock = 0; jblock < a.nblocks(); jblock++ )
    for ( int y = 0; y < a.nbs(jblock); y++ )
      jglobal[y + jblock*a.nb()] = a.j(jblock,y);

  // store addresses of columns of a and of u in acol and ucol
  vector<double*> acol(2*nploc);
  vector<double*> ucol(2*nploc);
  for ( int i = 0; i < a.nloc(); i++ )
    acol[i] = a.valptr(i*a.mloc());
  // if nloc is odd, store the address of vector 2*nploc-1
  if ( nloc_odd )
    acol[2*nploc-1] = &a_aux[0];
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
  while ( !done )
  {
    int npairs_above_threshold = 0;
    // sweep: process local pairs and rotate 2*np-1 times
    nsweep++;
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

      // perform Jacobi rotation for all local pairs

      // compute off-diagonal matrix elements apq for all pairs
      // skip the pair if one or both of the vectors is a dummy vector
      // i.e. a vector having jglobal==-1
      vector<double> apq(nploc);
      for ( int ipair = 0; ipair < nploc; ipair++ )
      {
        //cout << ctxt.mype() << ": computing apq for global pair "
        //     << jglobal[top[ipair]] << " " << jglobal[bot[ipair]] << endl;
        apq[ipair] = 0.0;
        if ( jglobal[top[ipair]] >= 0 && jglobal[bot[ipair]] >= 0 )
        {
          const double *ap = acol[top[ipair]];
          const double *uq = ucol[bot[ipair]];
          int one = 1;
          int mloc = a.mloc();
          //cout << ctxt.mype() << ": apq ddot: ipair=" << ipair << endl;
          //cout << ctxt.mype() << ": apq ddot: top=" << top[ipair]
          //     << " bot=" << bot[ipair] << endl;
          //cout << ctxt.mype() << ": ap=" << ap << " uq=" << uq << endl;

          //for ( int i = 0; i < mloc; i++ )
          //  cout << ap[i] << " " << uq[i] << endl;
          apq[ipair] = ddot(&mloc,ap,&one,uq,&one);
        }
      }
      // apq now contains partial sums of apq
      ctxt.dsum('c',nploc,1,&apq[0],nploc);
      // apq now contains the off-diagonal elements apq

      // compute the diagonal elements and perform the rotation only on
      // pairs having abs(apq)>threshold
      // skip pairs having dummy vectors
      vector<double> appqq(2*nploc);
      for ( int ipair = 0; ipair < nploc; ipair++ )
      {
        appqq[2*ipair]   = 0.0;
        appqq[2*ipair+1] = 0.0;
        if ( jglobal[top[ipair]] >= 0 &&
             jglobal[bot[ipair]] >= 0 &&
             fabs(apq[ipair]) > threshold )
        {
          // compute diagonal matrix elements
          const double *ap = acol[top[ipair]];
          const double *aq = acol[bot[ipair]];
          const double *up = ucol[top[ipair]];
          const double *uq = ucol[bot[ipair]];
          // compute matrix elements app, apq, aqq
          int one = 1;
          int mloc = a.mloc();
          appqq[2*ipair]   = ddot(&mloc,ap,&one,up,&one);
          appqq[2*ipair+1] = ddot(&mloc,aq,&one,uq,&one);
        }
      }
      // appqq now contains partial sums of app and aqq
      ctxt.dsum('c',2*nploc,1,&appqq[0],2*nploc);
      // appqq now contains the off-diagonal elements app and aqq

      for ( int ipair = 0; ipair < nploc; ipair++ )
      {
        // keep count of the number of pairs above threshold
        if ( jglobal[top[ipair]] >= 0 &&
             jglobal[bot[ipair]] >= 0 &&
             fabs(apq[ipair]) > threshold )
        {
          npairs_above_threshold++;

          // compute rotation sine and cosine
          const double app = appqq[2*ipair];
          const double aqq = appqq[2*ipair+1];

          // fabs(apq[ipair]) is larger than machine epsilon
          // since it is larger than the threshold (which is > macheps)
          const double tau = ( aqq - app ) / ( 2.0 * apq[ipair]);
          const double sq = sqrt( 1.0 + tau*tau );
          double t = 1.0 / ( fabs(tau) + sq );
          if ( tau < 0.0 ) t = -t;
          double c = 1.0 / sqrt(1.0+t*t);
          double s = t * c;

          // the drot function computes
          // c*x + s*y -> x
          //-s*x + c*y -> y
          // call drot with args c, -s
          // change the sign of s before call to drot
          s = -s;

          //cout << " p=" << jglobal[top[ipair]] << " q=" << jglobal[bot[ipair]]
          //     << " app=" << app << " aqq=" << aqq
          //     << " apq=" << apq[ipair] << endl;
          //cout << " tau=" << tau << " c=" << c << " s=" << s << endl;

          // update columns of a and u
          double *ap = acol[top[ipair]];
          double *aq = acol[bot[ipair]];
          double *up = ucol[top[ipair]];
          double *uq = ucol[bot[ipair]];
          int one = 1;
          int mloc = a.mloc();
          drot(&mloc,ap,&one,aq,&one,&c,&s);
          drot(&mloc,up,&one,uq,&one,&c,&s);
          //cout << " apq_check=" << ddot(&mloc,ap,&one,uq,&one);
          //cout << " aqp_check=" << ddot(&mloc,aq,&one,up,&one) << endl;
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

        int rbufi_left, rbufi_right, sbufi_left, sbufi_right;
        // send buffers contain a column of a and of u
        vector<double> sbuf_left(2*a.mloc()), sbuf_right(2*a.mloc());
        vector<double> rbuf_left(2*a.mloc()), rbuf_right(2*a.mloc());

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

        if ( ctxt.mycol() < last_active_process_col )
        {
          memcpy(&sbuf_right[0],    acol[bot[nploc-1]], mloc*sizeof(double) );
          memcpy(&sbuf_right[mloc], ucol[bot[nploc-1]], mloc*sizeof(double) );
          ctxt.dsend(2*mloc,1,&sbuf_right[0],2*mloc,ctxt.myrow(),ctxt.mycol()+1);
          ctxt.drecv(2*mloc,1,&rbuf_right[0],2*mloc,ctxt.myrow(),ctxt.mycol()+1);
          memcpy(acol[bot[nploc-1]], &rbuf_right[0],    mloc*sizeof(double) );
          memcpy(ucol[bot[nploc-1]], &rbuf_right[mloc], mloc*sizeof(double) );
          //cout << ctxt.mype() << ": received vector jglobal="
          //     << jglobal[bot[nploc-1]] << " from right" << endl;
        }
        if ( ctxt.mycol() != 0 )
        {
          memcpy(&sbuf_left[0],       acol[top[0]],     mloc*sizeof(double) );
          memcpy(&sbuf_left[mloc],    ucol[top[0]],     mloc*sizeof(double) );
          ctxt.dsend(2*mloc,1,&sbuf_left[0],2*mloc,ctxt.myrow(),ctxt.mycol()-1);
          ctxt.drecv(2*mloc,1,&rbuf_left[0],2*mloc,ctxt.myrow(),ctxt.mycol()-1);
          memcpy(acol[top[0]],       &rbuf_left[0],     mloc*sizeof(double) );
          memcpy(ucol[top[0]],       &rbuf_left[mloc],  mloc*sizeof(double) );
          //cout << ctxt.mype() << ": received vector jglobal="
          //     << jglobal[top[0]] << " from left" << endl;
        }
      } // if nploc > 0
      ctxt.barrier();
      // end of step

    } // for irot

    // sweep is complete
    // accumulate number of pairs above threshold
    ctxt.isum('r',1,1,&npairs_above_threshold,1);

    done = ( ( npairs_above_threshold == 0 ) || ( nsweep >= maxsweep ) );

  } // while !done
  // cout << " nsweep=" << nsweep << endl;

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
      memcpy(acol[idum],&a_aux[0],mloc*sizeof(double));
      memcpy(ucol[idum],&u_aux[0],mloc*sizeof(double));
    }
  }

  // compute eigenvalues
  for ( int i = 0; i < a.n(); i++ )
    e[i] = 0.0;
  for ( int jblock = 0; jblock < a.nblocks(); jblock++ )
    for ( int y = 0; y < a.nbs(jblock); y++ )
    {
      // j is the global column index
      int j = a.j(jblock,y);
      int jjj = y + jblock*a.nb();
      const double *ap = a.valptr(jjj*a.mloc());
      const double *up = u.valptr(jjj*u.mloc());
      int mloc = a.mloc();
      int one = 1;
      e[j] = ddot(&mloc,ap,&one,up,&one);
    }
  // e now contains the partial sums of the diagonal elements of a
  ctxt.dsum(a.n(),1,&e[0],a.n());
  // e contains the eigenvalues of a
  // u contains the eigenvectors of a

  return nsweep;
}
