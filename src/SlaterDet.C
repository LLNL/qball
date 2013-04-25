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
// SlaterDet.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Context.h"
#include "blas.h" // daxpy
#include "Base64Transcoder.h"
#include "SharedFilePtr.h"
#include "AtomSet.h"
#include "Species.h"
#include "Timer.h"
#include "PrintMem.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#if USE_CSTDIO_LFS
#include <sstream>
#include <cstdio>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(Context& ctxt, const Context& my_col_ctxt, D3vector kpoint,
                     bool ultrasoft, bool force_complex) : ctxt_(ctxt), c_(ctxt),
                                                           ctxtsq_(ctxt), spsi_(ctxt)
{
  ultrasoft_ = ultrasoft;
  force_complex_ = (ultrasoft_ || force_complex);
  basis_ = new Basis(my_col_ctxt,kpoint,force_complex_);
  gram_reshape_ = false;
  highmem_ = false;
  // set seed for randomization
  srand48(ctxt_.myproc());
}
////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(const SlaterDet& rhs) : ctxt_(rhs.context()),
  basis_(new Basis(*(rhs.basis_))), c_(rhs.c_), gram_reshape_(rhs.gram_reshape_),
  spsi_(rhs.spsi_), highmem_(rhs.highmem_), ultrasoft_(rhs.ultrasoft_) {}
////////////////////////////////////////////////////////////////////////////////
SlaterDet::~SlaterDet()  {
  delete basis_;
  for (int i=0; i<betag_.size(); i++) {
    delete betag_[i];
    delete betapsi_[i];
    delete dbetapsi_[i];
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::init_usfns(AtomSet* atoms) {
  atoms_ = atoms;
  const int nsp = atoms_->nsp();
  const int npcol = ctxt_.npcol();
  const int mycol = ctxt_.mycol();
  const int myrow = ctxt_.myrow();
  const int mb = basis_->maxlocalsize();
  const int m = ctxt_.nprow() * mb;

  // if this gets called twice, need to delete previous matrices
  for (int i=0; i<betag_.size(); i++) {
    delete betag_[i];
    delete betapsi_[i];
    delete dbetapsi_[i];
  }
  betag_.resize(nsp);
  betapsi_.resize(nsp);
  dbetapsi_.resize(nsp);
  for ( int is = 0; is < nsp; is++ ) {
    betag_[is] = new ComplexMatrix(ctxt_);
    betapsi_[is] = new ComplexMatrix(ctxt_);
    dbetapsi_[is] = new ComplexMatrix(ctxt_);
  }
  
  atoms_->usloc_atind.resize(nsp);
  atoms_->usloc_atind_t.resize(nsp);
  atoms_->usloc_nat.resize(nsp);
  atoms_->usloc_nat_t.resize(nsp);
  atoms_->naloc_max.resize(nsp);
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_->species_list[is];
    if (s->ultrasoft()) {
      int nbetalm = s->nbetalm();
      int na = atoms_->na(is);
      
      // naloc = maximum number of local atoms
      int naloc = na/npcol;
      if (na%npcol != 0) naloc++;
      atoms_->naloc_max[is] = naloc;
      
      atoms_->usloc_atind[is].clear();
      atoms_->usloc_atind_t[is].clear();
      for (int ia = 0; ia < na ; ia++) {
        int thisat = ia/naloc;
        if (thisat == mycol)
          atoms_->usloc_atind[is].push_back(ia);
        if (thisat == myrow)
          atoms_->usloc_atind_t[is].push_back(ia);
      }
      atoms_->usloc_nat[is] = atoms_->usloc_atind[is].size();
      atoms_->usloc_nat_t[is] = atoms_->usloc_atind_t[is].size();

      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.init_usfns, species " << is << ", nbetalm = " << nbetalm << ", naloc = " << atoms_->usloc_nat[is] << ", naloc_t = " << atoms_->usloc_nat_t[is] << endl;
#endif


      int ntot, nloc;
      if (highmem_) {
        ntot = na*nbetalm;
        nloc = naloc*nbetalm;
      }
      else {
        //ewd: bit of a hack, this puts a copy of betag_ on all columns
        ntot = nbetalm*ctxt_.npcol();
        nloc = nbetalm;
      }
      betag_[is]->resize(m,ntot,mb,nloc);

      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.init_usfns, species " << is << " betag matrix size = " << m << " x " << ntot << ", local size = " << mb << " x " << nloc << endl;
#endif
  
    }
  }

  calc_betag();
  calc_betapsi();
  
  for (int is=0; is<nsp; is++) {
    Species *s = atoms_->species_list[is];
    if (s->ultrasoft()) { 
      vector<complex<double> > qnm;
      vector<double> qaug;
      s->calc_qnmg(basis_,qnm,qaug);
      set_qaug(is,qaug);
    }
  }
  calc_spsi();
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_usfns() {
  calc_betapsi();
  calc_spsi();      
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::resize(const UnitCell& cell, const UnitCell& refcell,
  double ecut, int nst) {
  // Test in next line should be replaced by test on basis min/max indices
  // to signal change in basis vectors
  //if ( basis_->refcell().volume() != 0.0 && !refcell.encloses(cell) )
  //{
    //cout << " SlaterDet::resize: cell=" << cell;
    //cout << " SlaterDet::resize: refcell=" << basis_->refcell();
    //throw SlaterDetException("could not resize: cell not in refcell");
  //}

  try {
    // create a temporary copy of the basis and of the coefficient matrix
    Basis btmp(basis_->context(),basis_->kpoint(),force_complex_);
    btmp.resize(basis_->cell(),basis_->refcell(),basis_->ecut());
    
    ComplexMatrix ctmp(c_);
    
    // perform normal resize operations, possibly resetting contents of c_
    basis_->resize(cell,refcell,ecut);
    occ_.resize(nst);
    eig_.resize(nst);

    int mb = basis_->maxlocalsize();
    const int m = ctxt_.nprow() * mb;
    int nb = nst/ctxt_.npcol() + (nst%ctxt_.npcol() > 0 ? 1 : 0);

    //ewd:  hacky, but works for now
#ifdef BGQ 
    if (basis_->real())
    {
       while (mb%8 != 0)
          mb++;
       while (nb%8 != 0)
          nb++;
    }
    else
    {
       while (mb%4 != 0)
          mb++;
       while (nb%8 != 0)
          nb++;
    }
#endif

    // Determine if plane wave coefficients must be reset after the resize
    // This is needed if the dimensions of the matrix c_ must be changed
    const bool needs_reset =
      m!=c_.m() || nst!=c_.n() || mb!=c_.mb() || nb!=c_.nb();

    c_.resize(m,nst,mb,nb);

    if ( needs_reset )
      reset();
    
    // check if data can be copied from temporary copy
    // It is assumed that nst and ecut are not changing at the same time
    // Only the cases where one change at a time occurs is covered

    // consider only cases where the dimensions are finite
    if ( c_.m() > 0 && c_.n() > 0 ) {

      // first case: only nst has changed
      if ( c_.m() == ctmp.m() && c_.n() != ctmp.n() ) {
        //cout << "SlaterDet::resize: c_m/n=   "
        //     << c_.m() << "/" << c_.n() << endl;
        //cout << "SlaterDet::resize: ctmp_m/n=" << ctmp.m()
        //     << "/" << ctmp.n() << endl;
        // nst has changed, basis is unchanged
        // copy coefficients up to min(n_old, n_new)
        if ( c_.n() < ctmp.n() ) {
          c_.getsub(ctmp,ctmp.m(),c_.n(),0,0);
        }
        else {
          c_.getsub(ctmp,ctmp.m(),ctmp.n(),0,0);
        }
        gram();
      }
      // second case: basis was resized, nst unchanged
      if ( btmp.ecut() > 0.0 && basis_->ecut() > 0.0 &&
           c_.m() != ctmp.m() && c_.n() == ctmp.n() ) {
        // transform all states to real space and interpolate
        int np0 = max(basis_->np(0),btmp.np(0));
        int np1 = max(basis_->np(1),btmp.np(1));
        int np2 = max(basis_->np(2),btmp.np(2));
        //cout << " SlaterDet::resize: grid: np0/1/2: "
        //     << np0 << " " << np1 << " " << np2 << endl;
        // FourierTransform tf1(oldbasis, new basis grid)
        // FourierTransform tf2(newbasis, new basis grid)
        FourierTransform ft1(btmp,np0,np1,np2);
        FourierTransform ft2(*basis_,np0,np1,np2);
        // allocate real-space grid
        valarray<complex<double> > tmpr(ft1.np012loc());
        // transform each state from old basis to grid to new basis
        for ( int n = 0; n < nstloc(); n++ ) {
          ft1.backward(ctmp.cvalptr(n*ctmp.mloc()),&tmpr[0]);
          ft2.forward(&tmpr[0], c_.valptr(n*c_.mloc()));
        }
      }
    }
  }
  catch ( bad_alloc ) {
    cout << "<ERROR> bad_alloc exception caught in SlaterDet::resize </ERROR>" << endl;
    throw;
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::reshape(const Context& newctxt, const Context& new_col_ctxt, bool setnewctxt) {
  try {
    bool is_finite_ = ( c_.m() > 0 && c_.n() > 0 );

    Basis btmp(basis_->context(),basis_->kpoint(),ultrasoft_);  // copy of basis_
    btmp.resize(basis_->cell(),basis_->refcell(),basis_->ecut());
    Context oldctxt = ctxt_;

    ComplexMatrix ctmp(c_);                          // copy of c_
    // if ctxt.nprow() has changed, reconstruct Basis 
    if (new_col_ctxt.nprow() != basis_->context().nprow()) {
      UnitCell tmpcell = basis_->cell();
      UnitCell tmprefcell = basis_->refcell();
      double tmpecut = basis_->ecut();
      D3vector tmpkpt = basis_->kpoint();
      delete basis_;
      basis_ = new Basis(new_col_ctxt,tmpkpt,force_complex_);
      basis_->resize(tmpcell,tmprefcell,tmpecut);
    }
    int mb = basis_->maxlocalsize();
    const int m = newctxt.nprow() * mb;
    int nb = ctmp.n()/newctxt.npcol() + (ctmp.n()%newctxt.npcol() > 0 ? 1 : 0);

    //ewd:  hacky, but works for now
#ifdef BGQ 
    if (basis_->real())
    {
       while (mb%8 != 0)
          mb++;
       while (nb%8 != 0)
          nb++;
    }
    else
    {
       while (mb%4 != 0)
          mb++;
       while (nb%8 != 0)
          nb++;
    }
#endif

    if (setnewctxt)
      ctxt_ = newctxt;
    c_.set_context(newctxt);
    c_.resize(m,ctmp.n(),mb,nb);
    if (ultrasoft_)
    {
       spsi_.set_context(newctxt);
       spsi_.resize(m,ctmp.n(),mb,nb);
    }
    
    // consider only cases where the dimensions are finite
    if (is_finite_) {
      // transform all states to real space and interpolate
      int np0 = max(basis_->np(0),btmp.np(0));
      int np1 = max(basis_->np(1),btmp.np(1));
      int np2 = max(basis_->np(2),btmp.np(2));
      FourierTransform ft1(btmp,np0,np1,np2);
      FourierTransform ft2(*basis_,np0,np1,np2);
      assert(ft1.np012() == ft2.np012());

      const bool real_basis = basis_->real();
      const int wftmpr_size = real_basis ? ft1.np012() : 2*ft1.np012();
      const int wftmpr_loc_size1 = real_basis ? ft1.np012loc() : 2*ft1.np012loc();
      const int wftmpr_loc_size2 = real_basis ? ft2.np012loc() : 2*ft2.np012loc();
      vector<complex<double> > wftmp1(ft1.np012loc());
      vector<complex<double> > wftmp2(ft2.np012loc());
      vector<double> wftmpr(wftmpr_size);

      for ( int n = 0; n < nst(); n++ ) {
        // For each state:
        //    1.  transform to real space
        //    2.  collect full state on first pe in current column
        //    3.  redistribute to appropriate column of new context
        //    4.  transform back to reciprocal space

        // Barrier to limit the number of messages sent to first column tasks 
        // that don't have a receive posted

        oldctxt.barrier();
    
        bool pe0_sending = false;
        // check if state n resides on current process in old context
        if ( ctmp.pc(n) == oldctxt.mycol() ) {
          int nloc = ctmp.y(n); // local index
          ft1.backward(ctmp.cvalptr(nloc*ctmp.mloc()),&wftmp1[0]);
          double *a = (double*) &wftmp1[0];
          if ( real_basis ) 
            for ( int i = 0; i < wftmpr_loc_size1; i++ )
              wftmpr[i] = a[2*i];
          else
            for ( int i = 0; i < wftmpr_loc_size1; i++ )
              wftmpr[i] = a[i];


            /*
          if ( real_basis ) {
            double *a = (double*) &wftmp1[0];
            for ( int i = 0; i < ft1.np012loc(); i++ )
              wftmpr[i] = a[2*i];
          }
          else { 
            double *a = (double*) &wftmp1[0];
            for ( int i = 0; i < 2*ft1.np012loc(); i++ )
              wftmpr[i] = a[i];
            //memcpy((void*)&wftmpr[0],(void*)&wftmp1[0],
            //       ft1.np012loc()*sizeof(complex<double>));
          }
            */

          for ( int i = 0; i < oldctxt.nprow(); i++ ) {
            if ( i == oldctxt.myrow() ) {
              int size = wftmpr_loc_size1;
              oldctxt.isend(1,1,&size,1,0,ctmp.pc(n));
              oldctxt.dsend(wftmpr_loc_size1,1,&wftmpr[0],1,0,ctmp.pc(n));
            }
          }
          if ( oldctxt.myrow() == 0 ) {
            for ( int i = 0; i < oldctxt.nprow(); i++ ) {
              int size = 0;
              oldctxt.irecv(1,1,&size,1,i,ctmp.pc(n));
              int istart = ft1.np0() * ft1.np1() * ft1.np2_first(i);
              if ( !real_basis )
                istart *= 2;
              oldctxt.drecv(size,1,&wftmpr[istart],1,i,ctmp.pc(n));
            }
            // wftmpr is now complete on first task in column
            // send row and column number of current process 
            // for new context to pe 0
            int srcpe[2];
            srcpe[0] = newctxt.myrow();
            srcpe[1] = newctxt.mycol();
            if (srcpe[0] != 0 || srcpe[1] != 0)
              oldctxt.isend(2,1,&srcpe[0],1,0,0);
            else
              pe0_sending = true;  // skip the isend/irecv 
          }
        }

        // pe 0 handles traffic routing between contexts
        if (newctxt.myrow() == 0 && newctxt.mycol() == 0) {
          // receive location of process containing state n, wftmpr
          int srcpe[2] = {0,0};
          if (!pe0_sending)
            oldctxt.irecv(2,1,&srcpe[0],1,0,ctmp.pc(n));
          // send it to new process
          newctxt.isend(2,1,&srcpe[0],1,0,c_.pc(n));
        }

        // send the data
        if ( ctmp.pc(n) == oldctxt.mycol() ) {
          if ( oldctxt.myrow() == 0 ) {
            // send data to first pe on new process column
            int newpecol = c_.pc(n);
            int size = wftmpr_size;
            newctxt.isend(1,1,&size,1,0,newpecol);
            newctxt.dsend(size,1,&wftmpr[0],1,0,newpecol);
          }
        }

        // receive state n data on first process in new context column
        if ( c_.pc(n) == newctxt.mycol()) {
          if (newctxt.myrow() == 0) {
            int srcpe[2] = {0,0};
            newctxt.irecv(2,1,&srcpe[0],1,0,0);
            int srcrow = srcpe[0];
            int srccol = srcpe[1];
            int size = 0;
            newctxt.irecv(1,1,&size,1,srcrow,srccol);
            newctxt.drecv(size,1,&wftmpr[0],1,srcrow,srccol);
            for ( int i = 1; i < newctxt.nprow(); i++ ) {
              int size = -1;
              newctxt.irecv(1,1,&size,1,i,c_.pc(n));
              int istart = ft2.np0() * ft2.np1() * ft2.np2_first(i);
              if ( !real_basis )
                istart *= 2;
              if (size > 0) 
                newctxt.dsend(size,1,&wftmpr[istart],1,i,c_.pc(n));
            }
          }
          int size = -1;
          for ( int i = 0; i < newctxt.nprow(); i++ ) {
            if ( i == newctxt.myrow() ) {
              size = wftmpr_loc_size2;
              newctxt.isend(1,1,&size,1,0,c_.pc(n));
              if (size > 0) 
                newctxt.drecv(size,1,&wftmpr[0],1,0,c_.pc(n));
            }
          }
          for ( int i = 0; i < newctxt.nprow(); i++ ) {
            if ( i == newctxt.myrow() ) {
              if (i==0)
                size = wftmpr_loc_size2;

              // copy to complex array
              if (real_basis) 
                for ( int j = 0; j < ft2.np012loc(); j++ )
                  wftmp2[j] = complex<double>(wftmpr[j],0.0);
              else 
                for ( int j = 0; j < ft2.np012loc(); j++ )
                  wftmp2[j] = complex<double>(wftmpr[2*j],wftmpr[2*j+1]);
              int nloc = c_.y(n); // local index
              ft2.forward(&wftmp2[0], c_.valptr(nloc*c_.mloc()));
            }
          }
        }
      }
    }
  }
  catch ( bad_alloc ) {
    cout << "<ERROR> bad_alloc exception caught in SlaterDet::reshape </ERROR>" << endl;
    throw;
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::copyTo(SlaterDet* newsd) {
  // transfer data between two differently sized and/or dimensioned SlaterDets
  // without running into context scope problems
  //
  // currently written under assumption that (*this) is active on all tasks, i.e.
  // use for reshaping large SlaterDet (this) to smaller one (newsd).  
  
  try {
    bool is_finite_ = ( c_.m() > 0 && c_.n() > 0 );
    if (is_finite_) {
      // transform all states to real space and interpolate
      int np0 = basis_->np(0);
      int np1 = basis_->np(1);
      int np2 = basis_->np(2);
      FourierTransform ft1(*basis_,np0,np1,np2);

      const bool real_basis = basis_->real();
      const int wftmpr_size = real_basis ? ft1.np012() : 2*ft1.np012();
      const int wftmpr_loc_size1 = real_basis ? ft1.np012loc() : 2*ft1.np012loc();
      vector<complex<double> > wftmp1(ft1.np012loc());
      vector<double> wftmpr(wftmpr_size);

      for ( int n = 0; n < nst(); n++ ) {
        // For each state:
        //    1.  transform to real space
        //    2.  collect full state on first pe in current column
        //    3.  redistribute to appropriate column of new context
        //    4.  transform back to reciprocal space

        // Barrier to limit the number of messages sent to first column tasks 
        // that don't have a receive posted
        ctxt_.barrier();
    
        // check if state n resides on current process in old context
        const int statecol = c_.pc(n);
        if ( statecol == ctxt_.mycol() ) {
          int nloc = c_.y(n); // local index
          ft1.backward(c_.cvalptr(nloc*c_.mloc()),&wftmp1[0]);
          double *a = (double*) &wftmp1[0];
          if ( real_basis ) 
            for ( int i = 0; i < wftmpr_loc_size1; i++ )
              wftmpr[i] = a[2*i];
          else
            for ( int i = 0; i < wftmpr_loc_size1; i++ )
              wftmpr[i] = a[i];

          for ( int i = 0; i < ctxt_.nprow(); i++ ) {
            if ( i == ctxt_.myrow() ) {
              int size = wftmpr_loc_size1;
              ctxt_.isend(1,1,&size,1,0,statecol);
              ctxt_.dsend(wftmpr_loc_size1,1,&wftmpr[0],1,0,statecol);
            }
          }
          if ( ctxt_.myrow() == 0 ) {
            for ( int i = 0; i < ctxt_.nprow(); i++ ) {
              int size = 0;
              ctxt_.irecv(1,1,&size,1,i,statecol);
              int istart = ft1.np0() * ft1.np1() * ft1.np2_first(i);
              if ( !real_basis )
                istart *= 2;
              ctxt_.drecv(size,1,&wftmpr[istart],1,i,statecol);
            }
          }
          // wftmpr is now complete on first task in column
        }

        // determine appropriate pe of newctxt to send data
        int newproc = -1;
        int oldproc = -1;
        if (newsd != 0) {
          if ( newsd->context().active()) {
            Context& newctxt = newsd->context();
            const int newcol = newsd->c().pc(n);
            if (newctxt.mycol() == newcol && newctxt.myrow() == 0)
              newproc = 1;
            else
              newproc = 0;
          }
        }
        if (ctxt_.mycol() == statecol && ctxt_.myrow() == 0)
          oldproc = 1;

        int mype = ctxt_.mype();
        MPI_Status status;
        MPI_Send(&newproc,1,MPI_INT,0,mype,ctxt_.comm());
        MPI_Send(&oldproc,1,MPI_INT,0,mype,ctxt_.comm());

        // pe 0 will tell source and destination pes about each other
        int srcpe = -1;
        int destpe = -1;
        if (ctxt_.oncoutpe()) {
          for (int i=0; i<ctxt_.size(); i++) {
            int tmpnew, tmpold;
            MPI_Recv(&tmpnew,1,MPI_INT,i,i,ctxt_.comm(),&status);
            MPI_Recv(&tmpold,1,MPI_INT,i,i,ctxt_.comm(),&status);
            if (tmpnew == 1) destpe = i;
            if (tmpold == 1) srcpe = i;
          }
          MPI_Send(&srcpe,1,MPI_INT,destpe,destpe,ctxt_.comm());
          MPI_Send(&destpe,1,MPI_INT,srcpe,srcpe,ctxt_.comm());
        }
        if (newproc == 1)
          MPI_Recv(&srcpe,1,MPI_INT,0,mype,ctxt_.comm(),&status);
        if (oldproc == 1)
          MPI_Recv(&destpe,1,MPI_INT,0,mype,ctxt_.comm(),&status);

        if (! (oldproc == 1 && newproc == 1)) {
          if (oldproc == 1)
            MPI_Send(&wftmpr[0],wftmpr_size,MPI_DOUBLE,destpe,destpe,ctxt_.comm());
          if (newproc == 1)
            MPI_Recv(&wftmpr[0],wftmpr_size,MPI_DOUBLE,srcpe,mype,ctxt_.comm(),&status);
        }

        // now newsd can grab data, transform it back to reciprocal space
        if (newsd != 0) {
          if ( newsd->context().active()) {
            FourierTransform ft2(newsd->basis(),np0,np1,np2);
            assert(ft1.np012() == ft2.np012());
            const int wftmpr_loc_size2 = real_basis ? ft2.np012loc() : 2*ft2.np012loc();
            vector<complex<double> > wftmp2(ft2.np012loc());
            Context& newctxt = newsd->context();
            const int newcol = newsd->c().pc(n);

            if (newctxt.mycol() == newcol) {
              if (newctxt.myrow() == 0) {
                for ( int i = 1; i < newctxt.nprow(); i++ ) {
                  int size = -1;
                  newctxt.irecv(1,1,&size,1,i,newcol);
                  int istart = ft2.np0() * ft2.np1() * ft2.np2_first(i);
                  if ( !real_basis )
                    istart *= 2;
                  if (size > 0) 
                    newctxt.dsend(size,1,&wftmpr[istart],1,i,newcol);
                }
              }
              int size = -1;
              for ( int i = 0; i < newctxt.nprow(); i++ ) {
                if ( i == newctxt.myrow() ) {
                  size = wftmpr_loc_size2;
                  newctxt.isend(1,1,&size,1,0,newcol);
                  if (size > 0) 
                    newctxt.drecv(size,1,&wftmpr[0],1,0,newcol);
                }
              }
              for ( int i = 0; i < newctxt.nprow(); i++ ) {
                if ( i == newctxt.myrow() ) {
                  if (i==0)
                    size = wftmpr_loc_size2;
                  
                  // copy to complex array
                  if (real_basis) 
                    for ( int j = 0; j < ft2.np012loc(); j++ )
                      wftmp2[j] = complex<double>(wftmpr[j],0.0);
                  else 
                    for ( int j = 0; j < ft2.np012loc(); j++ )
                      wftmp2[j] = complex<double>(wftmpr[2*j],wftmpr[2*j+1]);
                  int nloc = newsd->c().y(n); // local index
                  int mloc = newsd->c().mloc();
                  ft2.forward(&wftmp2[0], newsd->c().valptr(nloc*mloc));
                }
              }
            }
          }
        }
      }
    }
  }
  catch ( bad_alloc ) {
    cout << "<ERROR> bad_alloc exception caught in SlaterDet::reshape </ERROR>" << endl;
    throw;
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::reset(void) {

  // initialize coefficients with lowest plane waves
  if ( c_.n() <= basis_->size() ) {
    // initialize c_
    c_.clear();
    const double s2i = 1.0 / sqrt(2.0);
 
    // for each n, find the smallest g vector and initialize
    int ismallest = 0;
    // on each process, basis.isort(ismallest) is the index of the smallest
    // local g vector
    for ( int n = 0; n < c_.n(); n++ ) {
      double value = 1.0;
      if ( basis().real() && n != 0 )
        value = s2i;
 
      // find process row holding the smallest g vector
      double kpg2 = basis_->kpg2(basis_->isort(ismallest));
      // cout << "smallest vector on proc " << ctxt_.mype()
      //      << " has norm " << kpg2 << endl;
      int minrow, mincol;
      ctxt_.dmin('c',' ',1,1,&kpg2,1,&minrow,&mincol,1,-1,-1);
 
      // find column hosting state n
      int pc = c_.pc(n);
      int pr = minrow;
      if ( pr == ctxt_.myrow() ) {
        int iii = basis_->isort(ismallest);
        ismallest++; // increment on entire process row
        if ( pc == ctxt_.mycol() ) {
          // cout << " n=" << n << " on process "
          //      << pr << "," << pc
          //      << " vector " << basis_->idx(3*iii) << " "
          //      << basis_->idx(3*iii+1) << " "
          //      << basis_->idx(3*iii+2) << " norm="
          //      << basis_->kpg2(iii) << " "
          //      << value << endl;
          int jjj = c_.m(n) * c_.nb() + c_.y(n);
          int index = iii+c_.mloc()*jjj;
          c_[index] = complex<double> (value,0.0);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::compute_density(FourierTransform& ft, 
  double weight, double* rho) const {

  //Timer tm_ft, tm_rhosum;
  // compute density of the states residing on my column of ctxt_
  assert(occ_.size() == c_.n());
  vector<complex<double> > tmp(ft.np012loc());
  
  assert(basis_->cell().volume() > 0.0);
  const double prefac = weight / basis_->cell().volume();  // weight = kpoint weight/total weightsum
  const int np012loc = ft.np012loc();
  
  //ewd:  add result to rho instead of overwriting to allow for averaging over local kpoints
  // (rho zeroed before SlaterDet.compute_density is called, i.e. in ChargeDensity.update_density)
  //for ( int i = 0; i < np012loc; i++ )
  //  rho[i] = 0.0;
  
  //ewd DEBUG:  transform one state at a time for both cases
  //if ( basis_->real() ) {
  if ( 1==0 && basis_->real() ) {
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ ) {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = prefac * occ_[nn];
      const double fac2 = prefac * occ_[nn+1];
      
      if ( fac1 + fac2 > 0.0 ) {
        //tm_ft.start();
        ft.backward(c_.cvalptr(n*c_.mloc()),c_.cvalptr((n+1)*c_.mloc()),&tmp[0]);
        //tm_ft.stop();
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        //tm_rhosum.start();
        for ( int i = 0; i < np012loc; i++ ) {
          const double psi1 = psi[ii];
          const double psi2 = psi[ii+1];
          rho[i] += fac1 * psi1 * psi1 + fac2 * psi2 * psi2;
          ii++; ii++;
        }
        //tm_rhosum.start();
      }
    }
    if ( nstloc() % 2 != 0 ) {
      const int n = nstloc()-1;
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = prefac * occ_[nn];
      
      if ( fac1 > 0.0 ) {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        for ( int i = 0; i < np012loc; i++ ) {
          const double psi1 = psi[ii];
          rho[i] += fac1 * psi1 * psi1;
          ii++; ii++;
        }
      }
    }
  }
  else {
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ ) {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac = prefac * occ_[nn];
      if ( fac > 0.0 ) {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);

        //ewd DEBUG:  try threading this loop:
        #pragma omp parallel for
        for ( int i = 0; i < np012loc; i++ )
           rho[i] += fac * (real(tmp[i])*real(tmp[i]) + imag(tmp[i])*imag(tmp[i]));
        //rho[i] += fac * norm(tmp[i]);
      }
    }
  }

  // cout << "SlaterDet: compute_density: ft_bwd time: " 
  //      << tm_ft.real() << endl;
  // cout << "SlaterDet: compute_density: rhosum time: " 
  //      << tm_rhosum.real() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::rs_mul_add(FourierTransform& ft, 
 const double* v, SlaterDet& sdp) const {

  // transform states to real space, multiply states by v[r] in real space
  // transform back to reciprocal space and add to sdp
  // sdp[n] += v * sd[n]
  
  vector<complex<double> > tmp(ft.np012loc());
  vector<complex<double> > ctmp(2*c_.mloc());
  
  const int np012loc = ft.np012loc();
  const int mloc = c_.mloc();
  double* p = (double*) &tmp[0];
  double* dcp = (double*) sdp.c().valptr();
  complex<double>* zcp = sdp.c().valptr();

  if ( basis_->real() ) {
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ ) {
      ft.backward(c_.cvalptr(n*mloc),
                 c_.cvalptr((n+1)*mloc),&tmp[0]);
      int ii = 0;
      for ( int i = 0; i < np012loc; i++ ) {
        const double psi1 = p[ii];
        const double psi2 = p[ii+1];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = vii * psi2;
        ii++; ii++;
      }
      ft.forward(&tmp[0], &ctmp[0], &ctmp[mloc]);
      int len = 4 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
    if ( nstloc() % 2 != 0 ) {
      const int n = nstloc()-1;
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);
      int ii = 0;
      for ( int i = 0; i < np012loc; i++ ) {
        const double psi1 = p[ii];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = 0.0;
        ii++; ii++;
      }
      ft.forward(&tmp[0], &ctmp[0]);
      int len = 2 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
  }
  else {
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ ) {
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);
      #pragma omp parallel for
      for ( int i = 0; i < np012loc; i++ )
        tmp[i] *= v[i];
      ft.forward(&tmp[0], &ctmp[0]);
      int len = mloc;
      int inc1 = 1;
      complex<double> alpha = complex<double>(1.0,0.0);
      zaxpy(&len,&alpha,&ctmp[0],&inc1,&zcp[n*mloc],&inc1);
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::gram() {

   //ewd DEBUG:  disable gram reshaping for now
   gram_reshape_ = false;

   
   //if (ultrasoft_)
   //   update_usfns();   // calculate betapsi, spsi

   if ( basis_->real() ) {
    // k = 0 case
    // create a DoubleMatrix proxy for c_
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sc_proxy(spsi_);
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

    if (ultrasoft_) {   // orthogonalize psi and S*psi
      s.gemm('t','n',2.0,sc_proxy,c_proxy,0.0);
      s.ger(-1.0,sc_proxy,0,c_proxy,0);
    }
    else {
      s.syrk('l','t',2.0,c_proxy,0.0); 
      s.syr('l',-1.0,c_proxy,0,'r');
    }
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
    // solve triangular system X * L^T = C
    c_proxy.trsm('r','l','t','n',1.0,s);

    /*
    // create a square context for the Cholesky decomposition
    int nsq = (int) sqrt((double) ctxt_.size());
    Context csq(nsq,nsq);
    DoubleMatrix ssq(csq,c_.n(),c_.n(),c_.nb(),c_.nb());
    ssq.getsub(s,s.m(),s.n(),0,0);
    ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
    s.getsub(ssq,s.m(),s.n(),0,0);
    // solve triangular system X * L^T = C
    c_proxy.trsm('r','l','t','n',1.0,s);
    */
  }
  else {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    if (ultrasoft_)  // orthogonalize psi and S*psi
      s.gemm('c','n',1.0,c_,spsi_,0.0);
    else
    {
       // gemm may actually be faster on BG/Q, try it 
#ifdef BGQ
       s.herk('l','c',1.0,c_,0.0);
       //s.gemm('c','n',1.0,c_,c_,0.0);
#else       
       s.herk('l','c',1.0,c_,0.0);
#endif
    }
    if (gram_reshape_) {
      int mbsq = c_.n()/ctxtsq_.nprow() + (c_.n()%ctxtsq_.nprow() == 0 ? 0 : 1);
      int nbsq = c_.n()/ctxtsq_.npcol() + (c_.n()%ctxtsq_.npcol() == 0 ? 0 : 1);
      
      if (mbsq > nbsq) 
        nbsq = mbsq;
      else
        mbsq = nbsq;
                
      ComplexMatrix ssq(ctxtsq_,c_.n(),c_.n(),nbsq,nbsq);
      ssq.getsub(s,s.m(),s.n(),0,0);
      ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
      s.getsub(ssq,ssq.m(),ssq.n(),0,0);
    }
    else {
      s.potrf('l'); // Cholesky decomposition: S = L * L^T
    }
    // solve triangular system X * L^T = C
    c_.trsm('r','l','c','n',1.0,s);

    //ewd:  let's have functions that use betapsi call this themselves
    //if (ultrasoft_)
    //  calc_betapsi(); // update betapsi w. new wfs
  }

   if (ultrasoft_)
      update_usfns();   // calculate betapsi, spsi

}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::set_gram_reshape(bool reshape) {
  gram_reshape_ = reshape;
  if (gram_reshape_) {
    // create most square rectangular context for Cholesky decomposition
    int tmpcol = ctxt_.size();
    int tmprow = 1;
    while (tmpcol > tmprow) {
      tmprow *= 2;
      tmpcol /= 2;
    }

    if (ctxt_.oncoutpe()) 
      cout << "<!-- SlaterDet.set_gram_reshape:  reshaping context from " << ctxt_.nprow() << " x " << ctxt_.npcol() << " to " << tmprow << " x " << tmpcol << " -->" << endl;
    
    ctxtsq_ = Context(ctxt_,tmprow,tmpcol);
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::riccati(SlaterDet& sd) {
  //ewd DEBUG
  if (ultrasoft_) cout << "ERROR:  ultrasoft not implemented w. riccati yet!!" << endl;
  
  if ( basis_->real() ) {
    // k = 0 case
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();
 
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // DoubleMatrix proxy for c_ and sd.c()
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());
    
    // Factor -1.0 in next line: -0.5 from definition of s, 2.0 for G and -G
    s.syrk('l','t',-1.0,c_proxy,0.5); // s = 0.5 * ( I - A )
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',0.5,c_proxy,0,'r');
    // factor -2.0 in next line: G and -G
    r.gemm('t','n',-2.0,sdc_proxy,c_proxy,1.0); // r = ( I - B )
    // rank-1 update using first row of sdc_proxy() and c_proxy
    r.ger(1.0,sdc_proxy,0,c_proxy,0);
 
    xm = s;
    xm.symmetrize('l');
 
    s.syrk('l','t',0.5,r,1.0); // s = s + 0.5 * r^T * r
    s.symmetrize('l');
 
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // x = s - 0.5 * ( r - xm )^T * ( r - xm )
      // Note: t and r are not symmetric, x, xm, and s are symmetric

      for ( int i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];
 
      x = s;
      x.syrk('l','t',-0.5,t,1.0);
 
      // get full matrix x
      x.symmetrize('l');
 
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];
 
      diff = t.nrm2();
 
      xm = x;
      iter++;
    }
    c_proxy.symm('r','l',1.0,x,sdc_proxy,1.0);
  }
  else {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();
 
    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // s = 0.5 * ( I - A )
    s.herk('l','c',-0.5,c_,0.5);
    // r = ( I - B )
    r.gemm('c','n',-1.0,sd.c(),c_,1.0);
 
    xm = s;
    xm.symmetrize('l');
 
    // s = s + 0.5 * r^H * r
    s.herk('l','c',0.5,r,1.0);
    s.symmetrize('l');
 
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // x = s - 0.5 * ( r - xm )^H * ( r - xm )
      // Note: t and r are not hermitian, x, xm, and s are hermitian

      for ( int i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];
 
      x = s;
      x.herk('l','c',-0.5,t,1.0);
      x.symmetrize('l');
 
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];
 
      diff = t.nrm2();
 
      xm = x;
      iter++;
    }
    c_.hemm('r','l',1.0,x,sd.c(),1.0);
  }
}
  
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::lowdin() {
  // Higham algorithm for polar decomposition
  if ( basis_->real() ) {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    l.clear();
    l.syrk('l','t',2.0,c_proxy,0.0); 
    l.syr('l',-1.0,c_proxy,0,'r');
    
    //cout << "SlaterDet::lowdin: A=\n" << l << endl;
    
    // Cholesky decomposition of A=Y^T Y
    l.potrf('l');
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::lowdin: L=\n" << l << endl;
    
    // Compute the polar decomposition of R = L^T

    x.transpose(1.0,l,0.0);
    // x now contains R
    //cout << "SlaterDet::lowdin: R=\n" << x << endl;
    
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );
      
 
      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];
 
      diff = t.nrm2();
      
      //cout << " SlaterDet::lowdin: diff=" << diff << endl;
      
      x = xp;
      //cout << "SlaterDet::lowdin: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor U of the 
    // polar decomposition R = UH
    
    //cout << " SlaterDet::lowdin: orthogonal polar factor=\n" << x << endl;
    
    // Compute L^-1
    l.trtri('l','n');
    // l now contains L^-1
    
    // Form the product L^-T U
    t.gemm('t','n',1.0,l,x,0.0);
    
    // Multiply c by L^-T U
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,t,0.0);
    
  }
  else { 
    //ewd can we use constructor to copy this, or should we copy with = operator?
    ComplexMatrix c_tmp(c_);     
    ComplexMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    l.clear();
    if (ultrasoft_)  // orthogonalize psi and S*psi
       l.gemm('c','n',1.0,c_,spsi_,0.0);
    else
       l.herk('l','c',1.0,c_,0.0);

    //cout << "SlaterDet::lowdin: A=\n" << l << endl;
    
    // Cholesky decomposition of A=Y^T Y
    l.potrf('l');
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::lowdin: L=\n" << l << endl;
    
    // Compute the polar decomposition of R = L^T

    x.transpose(1.0,l,0.0);
    // x now contains R
    //cout << "SlaterDet::lowdin: R=\n" << x << endl;
    
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );
      
 
      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];
 
      diff = t.nrm2();
      
      //cout << " SlaterDet::lowdin: diff=" << diff << endl;
      
      x = xp;
      //cout << "SlaterDet::lowdin: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor U of the 
    // polar decomposition R = UH
    
    //cout << " SlaterDet::lowdin: orthogonal polar factor=\n" << x << endl;
    
    // Compute L^-1
    l.trtri('l','n');
    // l now contains L^-1
    
    // Form the product L^-T U
    t.gemm('c','n',1.0,l,x,0.0);
    
    // Multiply c by L^-T U
    c_.gemm('n','n',1.0,c_tmp,t,0.0);
    
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::ortho_align(const SlaterDet& sd) {
  // Orthogonalize *this and align with sd
  // Higham algorithm for polar decomposition
  if ( basis_->real() ) {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sc_proxy(spsi_);
    DoubleMatrix sdc_proxy(sd.c());
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    l.clear();
    if (ultrasoft_) {   // orthogonalize psi and S*psi
      l.gemm('t','n',2.0,sc_proxy,c_proxy,0.0);
      l.ger(-1.0,sc_proxy,0,c_proxy,0);
    }
    else {
      l.syrk('l','t',2.0,c_proxy,0.0); 
      l.syr('l',-1.0,c_proxy,0,'r');
    }    
    //cout << "SlaterDet::ortho_align: A=\n" << l << endl;
    
    // Cholesky decomposition of A=Y^T Y
    l.potrf('l');
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::ortho_align: L=\n" << l << endl;
    
    // Compute the polar decomposition of L^-1 B
    // where B = C^T sd.C
    
    // Compute B: store result in x
    // factor -2.0 in next line: G and -G
    x.gemm('t','n',2.0,c_proxy,sdc_proxy,0.0);
    // rank-1 update using first row of sdc_proxy() and c_proxy
    x.ger(-1.0,c_proxy,0,sdc_proxy,0);
    
    // Form the product L^-1 B, store result in x
    // triangular solve: L X = B
    // trtrs: solve op(*this) * X = Z, output in Z
    l.trtrs('l','n','n',x);
    // x now contains L^-1 B

    //cout << "SlaterDet::ortho_align: L^-1 B=\n" << x << endl;
    
    // compute the polar decomposition of L^-1 B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );
      
 
      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];
 
      diff = t.nrm2();
      
      //cout << " SlaterDet::ortho_align: diff=" << diff << endl;
      
      x = xp;
      //cout << "SlaterDet::ortho_align: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor X of the 
    // polar decomposition L^-1 B = XH
    
    //cout << " SlaterDet::ortho_align: orthogonal polar factor=\n" 
    //     << x << endl;
    
    // Form the product L^-T Q
    // Solve trans(L) Z = X
    l.trtrs('l','t','n',x);
    
    // x now contains L^-T Q
    
    // Multiply c by L^-T Q
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,x,0.0);
    
  }
  else { 
    ComplexMatrix c_tmp(c_);
    ComplexMatrix sdc_proxy(sd.c());
    ComplexMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    l.clear();
    if (ultrasoft_)  // orthogonalize psi and S*psi
      l.gemm('c','n',1.0,spsi_,c_,0.0);
    else 
      l.herk('l','c',1.0,c_,0.0); 
    
    //cout << "SlaterDet::ortho_align: A=\n" << l << endl;
    
    // Cholesky decomposition of A=Y^T Y
    l.potrf('l');
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::ortho_align: L=\n" << l << endl;
    
    // Compute the polar decomposition of L^-1 B
    // where B = C^T sd.C
    
    // Compute B: store result in x
    if (ultrasoft_)
      x.gemm('c','n',1.0,spsi_,sdc_proxy,0.0);
    else
      x.gemm('c','n',1.0,c_,sdc_proxy,0.0);
    
    // Form the product L^-1 B, store result in x
    // triangular solve: L X = B
    // trtrs: solve op(*this) * X = Z, output in Z
    l.trtrs('l','n','n',x);
    // x now contains L^-1 B

    //cout << "SlaterDet::ortho_align: L^-1 B=\n" << x << endl;
    
    // compute the polar decomposition of L^-1 B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );
      
 
      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];
 
      diff = t.nrm2();
      
      //cout << " SlaterDet::ortho_align: diff=" << diff << endl;
      
      x = xp;
      //cout << "SlaterDet::ortho_align: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor X of the 
    // polar decomposition L^-1 B = XH
    
    //cout << " SlaterDet::ortho_align: orthogonal polar factor=\n" 
    //     << x << endl;
    
    // Form the product L^-T Q
    // Solve trans(L) Z = X
    l.trtrs('l','c','n',x);
    
    // x now contains L^-T Q
    
    // Multiply c by L^-T Q
    c_.gemm('n','n',1.0,c_tmp,x,0.0);
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::align(const SlaterDet& sd) {
  // Align *this with sd
  // Higham algorithm for polar decomposition
  if ( basis_->real() ) {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    
    // Compute the polar decomposition of B
    // where B = C^T sd.C
    
    // Compute B: store result in x
    // factor -2.0 in next line: G and -G
    x.gemm('t','n',2.0,c_proxy,sdc_proxy,0.0);
    // rank-1 update using first row of sdc_proxy() and c_proxy
    x.ger(-1.0,c_proxy,0,sdc_proxy,0);
    
    // x now contains B

    //cout << "SlaterDet::align: B=\n" << x << endl;
    
    // Compute the distance | c - sdc | before alignment
    //for ( int i = 0; i < c_proxy.size(); i++ )
    //  c_tmp_proxy[i] = c_proxy[i] - sdc_proxy[i];
    //cout << " SlaterDet::align: distance before: "
    //     << c_tmp_proxy.nrm2() << endl;
    
    // compute the polar decomposition of B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );
      
 
      // Next lines: use t as temporary to compute || x - xp ||_F
      //for ( int i = 0; i < t.size(); i++ )
      //  t[i] = x[i] - xp[i];
 
      //diff = t.nrm2();
      
      //cout << " SlaterDet::align: diff=" << diff << endl;
      
      x = xp;
      //cout << "SlaterDet::align: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor X of the 
    // polar decomposition B = XH
    
    //cout << " SlaterDet::align: orthogonal polar factor=\n" << x << endl;
        
    // Multiply c by X
    c_tmp_proxy = c_proxy;
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,x,0.0);
    
    // Compute the distance | c - sdc | after alignment
    //for ( int i = 0; i < c_proxy.size(); i++ )
    //  c_tmp_proxy[i] = c_proxy[i] - sdc_proxy[i];
    //cout << " SlaterDet::align: distance after:  "
    //     << c_tmp_proxy.nrm2() << endl;
    
  }
  else { 
    ComplexMatrix c_tmp(c_);
    ComplexMatrix sdc_proxy(sd.c());
    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // Compute the polar decomposition of B
    // where B = C^T sd.C
    
    // Compute B: store result in x
    x.gemm('t','n',1.0,c_,sdc_proxy,0.0);
    // x now contains B
    
    // compute the polar decomposition of B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon ) {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      // t now contains X^-T
      
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );

      x = xp;
      //cout << "SlaterDet::align: X=\n" << x << endl;
 
      iter++;
    }
    
    // x now contains the orthogonal polar factor X of the 
    // polar decomposition B = XH
    
    //cout << " SlaterDet::align: orthogonal polar factor=\n" << x << endl;
        
    // Multiply c by X
    c_tmp = c_;
    c_.gemm('n','n',1.0,c_tmp,x,0.0);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::dot(const SlaterDet& sd) const {
  // dot product of Slater determinants: dot = tr (V^T W)
  if ( basis_->real() ) {
    // DoubleMatrix proxy for c_ and sd.c()
    const DoubleMatrix c_proxy(c_);
    const DoubleMatrix sdc_proxy(sd.c());
    // factor 2.0: G and -G
    double d = 2.0 * c_proxy.dot(sdc_proxy);
    
    // correct double counting of first element
    double sum = 0.0;
    if ( ctxt_.myrow() == 0 ) {
      // compute the scalar product of the first rows of c_ and sd.c_
      const double *c = c_proxy.cvalptr(0);
      const double *sdc = sdc_proxy.cvalptr(0);
      int len = c_proxy.nloc();
      // stride of scalar product is mloc
      int stride = c_proxy.mloc();
      sum = ddot(&len,c,&stride,sdc,&stride);
    }
    ctxt_.dsum(1,1,&sum,1);
    return d - sum;
  }
  else {
    //ewdconst ComplexMatrix &c_proxy = c_;
    const ComplexMatrix &sdc_proxy = sd.c();
    complex<double> d = c_.dot(sdc_proxy);
    return d.real();
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::sdot(const SlaterDet& sd) const {
  // dot product of Slater determinant w. overlap matrix: dot = tr (V^T SW)
  if ( basis_->real() ) {
    // DoubleMatrix proxy for S*c_ and sd.c()
    const DoubleMatrix spsi_proxy(spsi_);
    const DoubleMatrix sdc_proxy(sd.c());
    // factor 2.0: G and -G
    double d = 2.0 * spsi_proxy.dot(sdc_proxy);
    
    // correct double counting of first element
    double sum = 0.0;
    if ( ctxt_.myrow() == 0 ) {
      // compute the scalar product of the first rows of c_ and sd.c_
      const double *sp = spsi_proxy.cvalptr(0);
      const double *sdc = sdc_proxy.cvalptr(0);
      int len = spsi_proxy.nloc();
      // stride of scalar product is mloc
      int stride = spsi_proxy.mloc();
      sum = ddot(&len,sp,&stride,sdc,&stride);
    }
    ctxt_.dsum(1,1,&sum,1);
    return d - sum;
  }
  else {
    //ewdconst ComplexMatrix &c_proxy = c_;
    const ComplexMatrix &sdc_proxy = sd.c();
    complex<double> d = spsi_.dot(sdc_proxy);
    return d.real();
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(int nel, int nspin) {
  // compute occupation numbers as 0.0, 1.0 or 2.0
  // if nspin = 1: use 0, 1 or 2
  // if nspin = 2: use 0 or 1;
  assert (nel >= 0);
  assert (occ_.size() == c_.n());
  if ( nspin == 1 ) {
    assert (nel <= 2*c_.n());
    int ndouble = nel/2;
    for ( int n = 0; n < ndouble; n++ )
      occ_[n] = 2.0;
    for ( int n = ndouble; n < ndouble+nel%2; n++ )
      occ_[n] = 1.0;
    for ( int n = ndouble+nel%2; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else if ( nspin == 2 ) {
    assert (nel <= c_.n());
    for ( int n = 0; n < nel; n++ )
      occ_[n] = 1.0;
    for ( int n = nel; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::total_charge(void) {
  // compute total charge from occ_[i]
  double sum = 0.0;
  for ( int n = 0; n < occ_.size(); n++ ) {
    sum += occ_[n];
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(int nspin, double mu, double temp, int ngauss) {
  // compute occupation numbers using either a Fermi distribution f(mu,temp)
  // or a Gaussian smearing (Methfessel-Paxton) scheme, and the eigenvalues in eig_[i]

  assert(nspin==1 || nspin==2);
  assert (occ_.size() == c_.n());
  assert (eig_.size() == c_.n());

  if (ngauss < 0) {  // Fermi distribution
    if ( nspin == 1 ) {
      for ( int n = 0; n < eig_.size(); n++ ) {
        occ_[n] = 2.0 * fermi(eig_[n],mu,temp);
      }
    }
    else if ( nspin == 2 ) {
      for ( int n = 0; n < eig_.size(); n++ ) {
        occ_[n] = fermi(eig_[n],mu,temp);
      }
    }
    else {
      // incorrect value of nspin_
      assert(false);
    }
  }
  else {
    if ( nspin == 1 ) {
      for ( int n = 0; n < eig_.size(); n++ ) {
        occ_[n] = 2.0 * methfessel(eig_[n],mu,temp,ngauss);
      }
    }
    else if ( nspin == 2 ) {
      for ( int n = 0; n < eig_.size(); n++ ) {
        occ_[n] = methfessel(eig_[n],mu,temp,ngauss);
      }
    }
    else {
      // incorrect value of nspin_
      assert(false);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::promote_occ(double occ_change, int origin_level, int destination_level)
{
   assert (occ_.size() == c_.n());
   assert (eig_.size() == c_.n());

   occ_[origin_level] -= occ_change;
   occ_[destination_level] += occ_change;

   assert ( (occ_[origin_level]>=0.0) && (occ_[origin_level]<=2.0) );
   assert ( (occ_[destination_level]>=0.0) && (occ_[destination_level]<=2.0) );
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::fermi(double e, double mu, double fermitemp) {
  // e, mu in Hartree, fermitemp in Rydbergs

  if ( fermitemp == 0.0 ) {
    if ( e < mu ) return 1.0;
    else if ( e == mu ) return 0.5;
    else return 0.0;
  }

  //const double kb = 3.1667907e-6; // Hartree/Kelvin
  //const double kt = kb * fermitemp;
  const double kt = 0.5 * fermitemp;  // convert to Hartree
  double arg = ( e - mu ) / kt;

  if ( arg < -30.0 ) return 1.0;
  if ( arg >  30.0 ) return 0.0;

  return 1.0 / ( 1.0 + exp ( arg ) );
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::methfessel(double e, double mu, double width, int ngauss) {
  // e, mu in Hartree, width in Rydbergs

  const double kt = 0.5 * width;  // convert to Hartree
  const double arg = (e - mu) / kt;

  // if n=0, this is equivalent to simple Gaussian smearing
  double sum = 0.5 - 0.5*erf(arg);

  // add additional terms in Methfessel-Paxton expansion   
  //   sum += A_n * H_(2n-1)(arg) * exp(-arg^2)
  //     where A_n = (-1)^n/(n! 4^n sqrt(pi))
  //           H_n(x) = Hermite polynomial of degree n
  //                   = 2x H_(n-1)(x) - 2n H_(n-2)(x)

  double an = 1.0/sqrt(3.141592653589793);
  double hodd = 0.0;     // H_(2n-1)(x)
  double heven = 1.0;    // H_2n(x)
  for (int n=1; n<=ngauss; n++) {
    double nd = (double)n;

    hodd = 2.0*arg*heven - 2.0*(nd-2.0)*hodd;  // H_(2n-1)(arg)
    heven = 2.0*arg*hodd - 2.0*(nd-1.0)*heven; // H_(2n)(arg)
    an *= -0.25/nd;
    sum += an * hodd * exp(-arg*arg);
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::entropy(int nspin) {
  // return dimensionless entropy
  // the contribution to the free energy is - t_kelvin * k_boltz * wf.entropy()

  assert(nspin==1 || nspin==2);
  const double fac = ( nspin > 1 ) ? 1.0 : 2.0;
  double sum = 0.0;
  for ( int n = 0; n < occ_.size(); n++ ) {
    const double f = occ_[n] / fac;
    if ( f > 0.0  &&  f < 1.0 ) {
      sum -= fac * ( f * log(f) + (1.0-f) * log(1.0-f) );
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::ortho_error(void) {
  // deviation from orthogonality of c_
  double error;
  if ( basis_->real() ) {
    // k = 0 case
    // declare a proxy DoubleMatrix for c_
    DoubleMatrix c_proxy(c_);
    
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // real symmetric rank-k update
    // factor 2.0 in next line: G and -G
    s.syrk('l','t',2.0,c_proxy,0.0); // compute real overlap matrix
 
    // correct for double counting of G=0
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',-1.0,c_proxy,0,'r');
 
    DoubleMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();
    
    s -= id; // subtract identity matrix from S
    
    error = s.nrm2();
  }
  else {
    // k != 0 case
    
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.herk('l','c',1.0,c_,0.0);
    
    ComplexMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();
 
    s -= id; // subtract identity matrix from S
    
    error = s.nrm2();
  }
  return error;
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::calc_betag() {
  // calculate betag on current basis, multiplied by structure factor and (-i)^l
  // call this function whenever atom positions or cell changes
  const int nsp = atoms_->nsp();
  const int ngwl = basis_->localsize();
  const int mloc = basis_->maxlocalsize();
  const double *kpg = basis_->kpg_ptr();
  const double *const kpgx = basis_->kpgx_ptr(0);
  const double *const kpgy = basis_->kpgx_ptr(1);
  const double *const kpgz = basis_->kpgx_ptr(2);
  const double omega = basis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  // calculate betag on wavefunction basis
  if (highmem_) {
    vector<vector<double> > tau;
    atoms_->get_positions(tau,true);
    for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_->species_list[is];
      if (s->ultrasoft()) {
        int nbeta = s->nbeta();
        int nbetalm = s->nbetalm();
        for (int ibl = 0; ibl < atoms_->usloc_atind[is].size(); ibl++) {
          int ia = atoms_->usloc_atind[is][ibl];

          // calculate structure factor
          vector<double> ckpgr(ngwl); 
          vector<double> skpgr(ngwl); 
          for ( int ig = 0; ig < ngwl; ig++ ) {
            const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
            skpgr[ig] = sin(arg);
            ckpgr[ig] = cos(arg);
          }
          double* bg = (double*) betag_[is]->valptr();
          int lm = 0;
          for (int b=0; b < nbeta; b++) {
            int nbetam = 2*s->betal(b)+1;
            for (int m=0; m < nbetam ; m++) {
              const int betal = s->betalm_l(lm);
              const int betam = s->betalm_m(lm);
              assert(betal == s->betal(b));
              const int ind0 = 2*ibl*nbetalm*mloc + 2*lm*mloc;
            
              if ( betal == 0 ) {
                for ( int ig = 0; ig < ngwl; ig++ ) {
                  double bval;
                  s->betag(b,kpg[ig],bval); 
                  const double ylmg = s->ylm(betal,betam,kpgx[ig],kpgy[ig],kpgz[ig]);
                  bg[ind0+2*ig]   = bval * ylmg * ckpgr[ig];
                  bg[ind0+2*ig+1] = -bval * ylmg * skpgr[ig];
                }
              }
              else if ( betal == 1 ) {
                for ( int ig = 0; ig < ngwl; ig++ ) {
                  double bval;
                  s->betag(b,kpg[ig],bval);
                  const double ylmg = s->ylm(betal,betam,kpgx[ig],kpgy[ig],kpgz[ig]);
                  bg[ind0+2*ig]   = -bval * ylmg * skpgr[ig];
                  bg[ind0+2*ig+1] = -bval * ylmg * ckpgr[ig];
                }
              }
              else if ( betal == 2 ) {
                for ( int ig = 0; ig < ngwl; ig++ ) {
                  double bval;
                  s->betag(b,kpg[ig],bval);
                  const double ylmg = s->ylm(betal,betam,kpgx[ig],kpgy[ig],kpgz[ig]);
                  bg[ind0+2*ig]   = -bval * ylmg * ckpgr[ig];
                  bg[ind0+2*ig+1] = bval * ylmg * skpgr[ig];
                }
              }
              else if ( betal == 3 ) {
                for ( int ig = 0; ig < ngwl; ig++ ) {
                  double bval;
                  s->betag(b,kpg[ig],bval);
                  const double ylmg = s->ylm(betal,betam,kpgx[ig],kpgy[ig],kpgz[ig]);
                  bg[ind0+2*ig]   = bval * ylmg * skpgr[ig];
                  bg[ind0+2*ig+1] = bval * ylmg * ckpgr[ig];
                }
              }
              lm++;
            }
          }
        }
      }
    }
  }
  else {  // !highmem, leave structure factor term out
    for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_->species_list[is];
      if (s->ultrasoft()) {
        complex<double>* bg = betag_[is]->valptr();
        int lm = 0;
        int nbeta = s->nbeta();
        for (int b=0; b < nbeta; b++) {
          int nbetam = 2*s->betal(b)+1;
          for (int m=0; m < nbetam ; m++) {
            const int betal = s->betalm_l(lm);
            const int betam = s->betalm_m(lm);
            assert(betal == s->betal(b));
            complex<double> il;
            if (betal == 0) il = complex<double>(1.0,0.0);
            else if (betal == 1) il = complex<double>(0.0,-1.0);
            else if (betal == 2) il = complex<double>(-1.0,0.0);
            else if (betal == 3) il = complex<double>(0.0,1.0);
            for ( int ig = 0; ig < ngwl; ig++ ) {
              double bval;
              s->betag(b,kpg[ig],bval); 
              const double ylmg = s->ylm(betal,betam,kpgx[ig],kpgy[ig],kpgz[ig]);
              bg[lm*mloc+ig] = bval*ylmg*il;
            }
            lm++;
          }
        }
      }
    }
  }

  // resize qaug using AtomSet
  if (qaug_.size() != nsp)
    qaug_.resize(nsp);

}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::calc_betapsi(void) {
  if (highmem_) {
    for (int is=0; is<betapsi_.size(); is++) {
      ComplexMatrix* bg = betag_[is];
      betapsi_[is]->resize(bg->n(),c_.n(),bg->nb(),c_.nb());

      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << " betapsi matrix size = " << betapsi_[is]->m() << " x " << betapsi_[is]->n() << ", local size = " << betapsi_[is]->mloc() << " x " << betapsi_[is]->nloc() << endl;
#endif
      
      if ( basis_->real() ) {
        // correct for G=0 term
        double* bgp = (double*) betag_[is]->valptr();
        int bg_nloc = betag_[is]->nloc();
        int bg_mloc = betag_[is]->mloc();
        // correct for G=0:  replace w. blas function (dscal/zscal)?
        //#pragma omp parallel for
        for (int ib=0; ib < bg_nloc; ib++) {
          bgp[2*ib*bg_mloc] *= 0.5;
          bgp[2*ib*bg_mloc+1] *= 0.5;
        }

      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", calling gemm of betag (" << bg->m() << " x " << bg->n() << ", " << bg->mloc() << " x " << bg->nloc() << ") and c (" << c_.m() << " x " << c_.n() << ", " << c_.mloc() << " x " << c_.nloc() << ")" << endl;
#endif


        betapsi_[is]->gemm('c','n',2.0,*bg,c_,0.0);
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", betapsi gemm finished." << endl;
#endif

        // restore betag:  replace w. blas function?
      //#pragma omp parallel for
        for (int ib=0; ib < bg_nloc; ib++) {
          bgp[2*ib*bg_mloc] *= 2.;
          bgp[2*ib*bg_mloc+1] *= 2.;
        }
      }
      else {
      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", calling gemm of betag (" << bg->m() << " x " << bg->n() << ", " << bg->mloc() << " x " << bg->nloc() << ") and c (" << c_.m() << " x " << c_.n() << ", " << c_.mloc() << " x " << c_.nloc() << ")" << endl;
#endif
      
        betapsi_[is]->gemm('c','n',1.0,*bg,c_,0.0);
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", betapsi gemm finished." << endl;
#endif
      }
    }
  }
  else {
    vector<vector<double> > tau;
    atoms_->get_positions(tau,true);
    int nsp = atoms_->nsp();
    for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_->species_list[is];
      if (s->ultrasoft()) {
        int nbeta = s->nbeta();
        int nbetalm = s->nbetalm();
        const int ngwl = basis_->localsize();
        ComplexMatrix* bg = betag_[is];
        const complex<double>* bgp = bg->cvalptr();

        int na = atoms_->na(is);
        int npcol = ctxt_.npcol();
        int bp_m = na*nbetalm;
        int bp_mloc = atoms_->naloc_max[is]*nbetalm;
        betapsi_[is]->resize(bp_m,c_.n(),bp_mloc,c_.nb());
        complex<double>* bpp = betapsi_[is]->valptr();
        int bg_mloc = betag_[is]->mloc();
        int naloc = atoms_->usloc_nat[is];
        int naloc_t = atoms_->usloc_nat_t[is];
        int naloc_max = atoms_->naloc_max[is];

        // create temporary matrices to do parallel gemms for each local atom
        int bpt_m = nbetalm*npcol;
        ComplexMatrix bptmp(ctxt_,bpt_m,c_.n(),nbetalm,c_.nb());
        complex<double>* bptmpp = bptmp.valptr();
        ComplexMatrix bgsf(ctxt_,betag_[is]->m(),betag_[is]->n(),betag_[is]->mb(),betag_[is]->nb());
        complex<double>* bgsfp = bgsf.valptr();
        for (int ibl = 0; ibl < naloc_max; ibl++) {
          if (ibl < naloc) { 
            int ia = atoms_->usloc_atind[is][ibl]; 
            // calculate structure factor
            vector<double> ckpgr(ngwl); 
            vector<double> skpgr(ngwl); 
            const double *const kpgx = basis_->kpgx_ptr(0);
            const double *const kpgy = basis_->kpgx_ptr(1);
            const double *const kpgz = basis_->kpgx_ptr(2);
            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
              skpgr[ig] = sin(arg);
              ckpgr[ig] = cos(arg);
            }

            for (int lm=0; lm<nbetalm; lm++)
              for ( int ig = 0; ig < ngwl; ig++ )
                bgsfp[lm*bg_mloc+ig] = bgp[lm*bg_mloc+ig]*complex<double>(ckpgr[ig],-skpgr[ig]);
          }
      //ewd DEBUG
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", ibl = " << ibl << ", calling gemm of bgsf (" << bgsf.m() << " x " << bgsf.n() << ", " << bgsf.mloc() << " x " << bgsf.nloc() << ") and c (" << c_.m() << " x " << c_.n() << ", " << c_.mloc() << " x " << c_.nloc() << ")" << endl;
#endif
      
          bptmp.gemm('c','n',1.0,bgsf,c_,0.0);
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_betapsi, species " << is << ", bptmp gemm finished." << endl;
#endif

          if (ibl < naloc_t) { 
            for ( int n = 0; n < nstloc(); n++ ) {
              for (int lm = 0; lm < nbetalm; lm++) {
                int bpind = bp_mloc*n+ibl*nbetalm+lm;
                bpp[bp_mloc*n+ibl*nbetalm+lm] = bptmpp[nbetalm*n+lm];
              }
            }
          }
        }
      }
    }
  }

  /*
  //ewd DEBUG
  vector<vector<double> > tau;
  atoms_->get_positions(tau,true);
  int nsp = atoms_->nsp();
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_->species_list[is];
    if (s->ultrasoft()) {
      int nbetalm = s->nbetalm();
      int naloc_t = atoms_->usloc_nat_t[is];
      int naloc_max = atoms_->naloc_max[is];
      int bp_mloc = atoms_->naloc_max[is]*nbetalm;
      complex<double>* bpp = betapsi_[is]->valptr();
      for (int ibl = 0; ibl < naloc_max; ibl++) {
        int ia = atoms_->usloc_atind_t[is][ibl]; 
        if (ibl < naloc_t) { 
          for ( int n = 0; n < nstloc(); n++ ) {
            for (int lm = 0; lm < nbetalm; lm++) {
              cout << "SD.BETAPSI, mype = " << ctxt_.mype() << ", is = " << is << ", ibl = " << ibl << ", ia = " << ia << ", n = " << n << ", lm = " << lm << ", bpsi = " << bpp[bp_mloc*n+ibl*nbetalm+lm] << endl;
            }
          }
        }
      }
    }
  }
  //ewd DEBUG
  */

}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::calc_dbetapsi(int j) {
  const int ngwl = basis_->localsize();
  const int nsp = betapsi_.size();
  if (highmem_) {
    for (int is=0; is<nsp; is++) {
      ComplexMatrix* bg = betag_[is];
      complex<double>* bgp = betag_[is]->valptr();
      ComplexMatrix dbgj(ctxt_,betag_[is]->m(),betag_[is]->n(),betag_[is]->mb(),betag_[is]->nb());
      dbetapsi_[is]->resize(bg->n(),c_.n(),bg->nb(),c_.nb());
      dbetapsi_[is]->clear();
      int bg_nloc = betag_[is]->nloc();
      int bg_mloc = betag_[is]->mloc();
      complex<double>* dbgpj = dbgj.valptr();
      const double *const kpgj = basis_->kpgx_ptr(j);
      for (int ib=0; ib < bg_nloc; ib++)
        for (int ig=0; ig<ngwl; ig++)
          dbgpj[ib*bg_mloc+ig] = bgp[ib*bg_mloc+ig]*complex<double>(0.0,-kpgj[ig]);

      if ( basis_->real() ) {
        // correct for G=0 term:  replace w. blas function?
         //#pragma omp parallel for
        for (int ib=0; ib < bg_nloc; ib++)
          dbgpj[ib*bg_mloc] *= 0.5;
        dbetapsi_[is]->gemm('c','n',2.0,dbgj,c_,0.0);
      }
      else {
        dbetapsi_[is]->gemm('c','n',1.0,dbgj,c_,0.0);
      }
    }
  }
  else {
    vector<vector<double> > tau;
    atoms_->get_positions(tau,true);
    int nsp = atoms_->nsp();
    for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_->species_list[is];
      if (s->ultrasoft()) {
        int nbeta = s->nbeta();
        int nbetalm = s->nbetalm();
        ComplexMatrix* bg = betag_[is];
        const complex<double>* bgp = bg->cvalptr();

        int na = atoms_->na(is);
        int npcol = ctxt_.npcol();
        int bp_m = na*nbetalm;
        int bp_mloc = atoms_->naloc_max[is]*nbetalm;
        dbetapsi_[is]->resize(bp_m,c_.n(),bp_mloc,c_.nb());
        dbetapsi_[is]->clear();
        complex<double>* dbp = dbetapsi_[is]->valptr();
        int bg_mloc = betag_[is]->mloc();
        int naloc = atoms_->usloc_nat[is];
        int naloc_t = atoms_->usloc_nat_t[is];
        int naloc_max = atoms_->naloc_max[is];

        // create temporary matrices to do parallel gemms for each local atom
        int bpt_m = nbetalm*npcol;
        ComplexMatrix bptmp(ctxt_,bpt_m,c_.n(),nbetalm,c_.nb());
        complex<double>* bptmpp = bptmp.valptr();
        ComplexMatrix dbgsf(ctxt_,betag_[is]->m(),betag_[is]->n(),betag_[is]->mb(),betag_[is]->nb());
        complex<double>* dbgsfp = dbgsf.valptr();
        
        for (int ibl = 0; ibl < naloc_max; ibl++) {
          if (ibl < naloc) { 
            int ia = atoms_->usloc_atind[is][ibl]; 
            // calculate structure factor
            vector<double> ckpgr(ngwl); 
            vector<double> skpgr(ngwl); 
            const double *const kpgj = basis_->kpgx_ptr(j);
            const double *const kpgx = basis_->kpgx_ptr(0);
            const double *const kpgy = basis_->kpgx_ptr(1);
            const double *const kpgz = basis_->kpgx_ptr(2);
            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
              skpgr[ig] = sin(arg);
              ckpgr[ig] = cos(arg);
            }

            for (int lm=0; lm<nbetalm; lm++)
              for ( int ig = 0; ig < ngwl; ig++ ) 
                dbgsfp[lm*bg_mloc+ig] = bgp[lm*bg_mloc+ig]*complex<double>(-kpgj[ig]*skpgr[ig],-kpgj[ig]*ckpgr[ig]);

          }
          bptmp.gemm('c','n',1.0,dbgsf,c_,0.0);

          if (ibl < naloc_t) { 
            for ( int n = 0; n < nstloc(); n++ ) {
              for (int lm = 0; lm < nbetalm; lm++) {
                int bpind = bp_mloc*n+ibl*nbetalm+lm;
                dbp[bp_mloc*n+ibl*nbetalm+lm] = bptmpp[nbetalm*n+lm];
              }
            }
          }
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::calc_wfphase() {
  // this is a debugging function for the case of gamma-pt calcs w. a complex wf
  for ( int n = 0; n < c_.nloc(); n++ ) {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    int phasecnt = 0;
    double phaseavg = 0.0;
    double resum = 0.0;
    double imsum = 0.0;
    for ( int i = 0; i < basis_->localsize(); i++ ) {
      double re = real(p[i]);
      double im = imag(p[i]);
      resum += re;
      imsum += im;
      if (re != 0.0) {
        phaseavg += im/re;
        phasecnt++;
      }
    }
    if (ctxt_.oncoutpe())
      cout << "SD.CALC_WFPHASE, local state " << n << ", phaseavg = " << phaseavg/phasecnt << ", phasecnt = " << phasecnt << ", resum = " << resum << ", imsum = " << imsum << endl;
  }  
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::set_qaug(int is, vector<double>& qaug) {
  // store qaug when it's calculated in NonLocalPotential to avoid recalculation

  assert(is < qaug_.size());
  const int nbetalmsq = qaug.size();
  if (qaug_[is].size() != nbetalmsq) {
    qaug_[is].resize(nbetalmsq);
  }
  for (int in=0; in<nbetalmsq; in++)
    qaug_[is][in] = qaug[in];
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::calc_spsi() {

  // ultrasoft:  calculate product of overlap matrix S w. wavefunction
  // S*|psi> = |psi> + SUM qaug*|beta_n><beta_m|psi>
  const double omega = basis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;
  vector<vector<double> > tau;
  if (!highmem_)
    atoms_->get_positions(tau,true);

  spsi_.resize(c_.m(),c_.n(),c_.mb(),c_.nb());
  spsi_ = c_;
  
  int nsp = atoms_->nsp();
  for (int is=0; is<nsp; is++) {
    Species *s = atoms_->species_list[is];
    if (s->ultrasoft()) { 
      int naloc_t = atoms_->usloc_nat_t[is];
      int naloc = atoms_->usloc_nat[is];
      int naloc_max = atoms_->naloc_max[is];
      int nqtot = s->nqtot();
      int nbetalm = s->nbetalm();

      ComplexMatrix* betag = betag_[is];
      ComplexMatrix* bpsi = betapsi_[is];
      const int bp_mloc = bpsi->mloc();
      const int bg_mloc = betag_[is]->mloc();
      const int sp_mloc = spsi_.mloc();
      const int nstloc = bpsi->nloc();
      complex<double>* bp = bpsi->valptr();
      if (highmem_) {
        // multiply qaug_nm and <beta|psi>
        ComplexMatrix bpsisum(betapsi_[is]->context(),betapsi_[is]->m(),betapsi_[is]->n(),betapsi_[is]->mb(),betapsi_[is]->nb());
        complex<double>* bpsp = bpsisum.valptr();
        if (bpsi->size() > 0) { 
          vector<complex<double> > bpsum(nstloc*bp_mloc);
          int bpsize = bpsum.size();
          for (int i=0; i<bpsize; i++)
            bpsum[i] = complex<double>(0.0,0.0);
          for (int n=0; n<nstloc; n++) {
            for (int ibl = 0; ibl < naloc_t; ibl++) {
              int ia = atoms_->usloc_atind_t[is][ibl];
              for (int qind=0; qind<nqtot; qind++) {
                int lm1 = s->qnm_lm1(qind);
                int lm2 = s->qnm_lm2(qind);
                int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                bpsum[ind1] += qaug_[is][qind]*bp[ind2];
                if (lm1 != lm2) 
                  bpsum[ind2] += qaug_[is][qind]*bp[ind1];
                
              }
            }
          }
          for (int i=0; i<bpsize; i++)
            bpsp[i] = omega_inv*bpsum[i];
        }
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_spsi, species " << is << ", calling gemm of betag (" << betag->m() << " x " << betag->n() << ", " << betag->mloc() << " x " << betag->nloc() << ") and bpsisum (" << bpsisum.m() << " x " << bpsisum.n() << ", " << bpsisum.mloc() << " x " << bpsisum.nloc() << ")" << endl;
#endif
        spsi_.gemm('n','n',1.0,*betag,bpsisum,1.0);
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_spsi, species " << is << ", spsi gemm finished." << endl;
#endif
      }
      else { // !highmem
        // need to multiply betag by structure factor first
        const double *const kpgx = basis_->kpgx_ptr(0);
        const double *const kpgy = basis_->kpgx_ptr(1);
        const double *const kpgz = basis_->kpgx_ptr(2);
        const int ngwl = basis_->localsize();
        const complex<double>* bgp = betag_[is]->cvalptr();
        complex<double>* sp = spsi_.valptr();

        // create temporary matrices to do parallel gemms for each local atom
        ComplexMatrix bgsf(ctxt_,betag_[is]->m(),betag_[is]->n(),betag_[is]->mb(),betag_[is]->nb());
        complex<double>* bgsfp = bgsf.valptr();
        int tbp_m = nbetalm*ctxt_.npcol();
        ComplexMatrix tbpsum(ctxt_,tbp_m,c_.n(),nbetalm,c_.nb());
        complex<double>* tbpp = tbpsum.valptr();
        for (int ibl = 0; ibl < naloc_max; ibl++) {
          vector<complex<double> > bpsum(nstloc*nbetalm);
          for (int i=0; i<bpsum.size(); i++)
            bpsum[i] = complex<double>(0.0,0.0);

          if (ibl < naloc) {
            int ia = atoms_->usloc_atind[is][ibl];
            vector<double> ckpgr(ngwl); 
            vector<double> skpgr(ngwl); 
            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
              skpgr[ig] = sin(arg);
              ckpgr[ig] = cos(arg);
            }
          
            for (int lm=0; lm<nbetalm; lm++)
              for ( int ig = 0; ig < ngwl; ig++ )
                bgsfp[lm*bg_mloc+ig] = bgp[lm*bg_mloc+ig]*complex<double>(ckpgr[ig],-skpgr[ig]);
          }
          //else {
          //  for (int lm=0; lm<nbetalm; lm++)
          //    for ( int ig = 0; ig < ngwl; ig++ )
          //      bgsfp[lm*bg_mloc+ig] = complex<double>(0.0,0.0);
          //}

          if (ibl < naloc_t) { 
            // copy local atom's betapsi data to tbpsum
            if (bpsi->size() > 0) { 
              for (int n=0; n<nstloc; n++) {
                for (int qind=0; qind<nqtot; qind++) {
                  int lm1 = s->qnm_lm1(qind);
                  int lm2 = s->qnm_lm2(qind);
                  int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                  int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                  int tind1 = nbetalm*n + lm1;
                  int tind2 = nbetalm*n + lm2;
                  bpsum[tind1] += qaug_[is][qind]*bp[ind2];
                  if (lm1 != lm2) 
                    bpsum[tind2] += qaug_[is][qind]*bp[ind1];
                  
                }
              }
            }
          }
          if (tbpsum.size() >  0) {
            for (int i=0; i<nstloc*nbetalm; i++)
              tbpp[i] = omega_inv*bpsum[i];
          }
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_spsi, species " << is << ", calling gemm of bgsf (" << bgsf.m() << " x " << bgsf.n() << ", " << bgsf.mloc() << " x " << bgsf.nloc() << ") and tbpsum (" << tbpsum.m() << " x " << tbpsum.n() << ", " << tbpsum.mloc() << " x " << tbpsum.nloc() << ")" << endl;
#endif
          spsi_.gemm('n','n',1.0,bgsf,tbpsum,1.0);
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "SD.calc_spsi, species " << is << ", spsi gemm finished." << endl;
#endif
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::randomize(double amplitude) {
  // Note: randomization results depend on the process grid size and shape

  srand48(ctxt_.myproc());
  for ( int n = 0; n < c_.nloc(); n++ ) {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    for ( int i = 0; i < basis_->localsize(); i++ ) {
      double re = drand48();
      double im = drand48();
      p[i] += amplitude * complex<double>(re,im);
      //p[i] = amplitude * complex<double>(re,im);
    }
  }
  cleanup();
  gram();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::randomize_us(double amplitude, AtomSet& as) {
  // Note: randomization results depend on the process grid size and shape

  //ewd DEBUG
  if (ctxt_.mype() == 0) {
    if (highmem_)
      cout << "<!-- Randomize_wf_us:  highmem ON -->" << endl;
    else
      cout << "<!-- Randomize_wf_us:  highmem OFF -->" << endl;
  }
  //ewd DEBUG

  
  srand48(ctxt_.myproc());
  for ( int n = 0; n < c_.nloc(); n++ ) {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    for ( int i = 0; i < basis_->localsize(); i++ ) {
      double re = drand48();
      double im = drand48();
      //p[i] += amplitude * complex<double>(re,im);
      p[i] = amplitude * complex<double>(re,im);
    }
  }
  if (ctxt_.mype() == 0)
     cout << "Calling cleanup()" << endl;

  cleanup();

  if (ctxt_.mype() == 0)
     cout << "Calling init_usfns()" << endl;

  // calculate ultrasoft functions
  init_usfns(&as);

  if (ctxt_.mype() == 0)
     cout << "Calling gram()" << endl;

  // orthogonalize
  gram();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::randomize_real(double amplitude)
{
  if ( basis_->size() == 0 )
    return;
  // Note: randomization results depend on the process grid size and shape
  srand48(ctxt_.myproc());
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    for ( int i = 0; i < basis_->localsize(); i++ )
    {
      double re = drand48();
      double im = drand48();
      p[i] = amplitude * complex<double>(re,im);
      p[basis_->index_of_minus_g(i)] = amplitude * complex<double>(re,-1.0*im);
    }

    p[basis_->zero_index()] = complex<double>(real(p[basis_->zero_index()]),0.0);
  }
  // gram does an initial cleanup
  gram();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::rescale(double factor) {
  c_ *= factor;
}

////////////////////////////////////////////////////////////////////////////////
// AS: shift state n_state by the vector (shift_x, shift_y, shift_z)
void SlaterDet::shift_wf(double shift_x,double shift_y,double shift_z,int n_state)
{
  if ( basis_->size() == 0 )
    return;

  // AS: remark: nstloc() is the same as c_.nloc()

  // AS: first version: shifts all the states
  // for ( int n = 0; n < nstloc(); n++ )
  //   {
  //     for ( int i = 0; i < basis_->localsize(); i++ )
  //       {
  //         double kpg_x = basis_->kpgx(i);
  //         double kpg_y = basis_->kpgx(i+basis_->localsize()  );
  //         double kpg_z = basis_->kpgx(i+basis_->localsize()*2);
  //         c_[i+n*c_.mloc()] *= exp( complex<double>(0,1) * (kpg_x*shift_x+kpg_y*shift_y+kpg_z*shift_z) );
  //       }
  //   }

  // AS: second version: shifts only n_state, worked
  // for ( int i = 0; i < basis_->localsize(); i++ )
  //   {
  //     double kpg_x = basis_->kpgx(i);
  //     double kpg_y = basis_->kpgx(i+basis_->localsize()  );
  //     double kpg_z = basis_->kpgx(i+basis_->localsize()*2);
  //     c_[i+(n_state-1)*c_.mloc()] *= exp( complex<double>(0,1) * (kpg_x*shift_x+kpg_y*shift_y+kpg_z*shift_z) );
  //   }

  // AS: final version: shifts only n_state, works as well

  if ( ctxt_.mycol() != c_.pc(n_state-1) )
    return;

  assert( c_.y(n_state-1) <= c_.nloc() );
  complex<double>* p = c_.valptr(c_.mloc()*c_.y(n_state-1));

  for ( int i = 0; i < basis_->localsize(); i++ )
  {
    double kpg_x = basis_->kpgx(i);
    double kpg_y = basis_->kpgx(i+basis_->localsize());
    double kpg_z = basis_->kpgx(i+basis_->localsize()*2);
    p[i] *= exp( complex<double>(0,1) * (kpg_x*shift_x+kpg_y*shift_y+kpg_z*shift_z) );
  }

  // AS: no orthogonalization should be necessary here
  // gram();
}

////////////////////////////////////////////////////////////////////////////////
// AS: change phase of the wave function to make it real for Gamma only
void SlaterDet::phase_wf_real(void)
{
  if ( basis_->size() == 0 )
    return;

  for ( int n = 0; n < nstloc(); n++ )
  {
    // AS: pick the plane-wave coefficient
    complex<double>* p = c_.valptr( c_.mloc() * n );
    double as_renorm = arg( p[basis_->zero_index()] );

    // AS: This would be an alternative way if doing the phase change, but it does give the same result
    // complex<double> as_renorm2 = abs( p[basis_->zero_index()] ) / p[basis_->zero_index()];

    // AS: DEBUG output
    // cout << "AS: ZERO_INDEX: " << basis_->zero_index() << endl;
    // cout.precision(20);
    // cout << "AS: STATE: " << n << " p[basis_->zero_index()]        : " << p[basis_->zero_index()] << endl;
    // cout << "AS: STATE: " << n << " arg( p[basis_->zero_index()] ) : " << arg(p[basis_->zero_index()]) << endl;
    // cout << "AS: STATE: " << n << " as_renorm                      : " << as_renorm << endl;
    // cout << "AS: STATE: " << n << " as_renorm2                     : " << as_renorm2 << endl;
    // cout << "AS: STATE: " << n << " abs( p[basis_->zero_index()])  : " << abs(p[basis_->zero_index()]) << endl;
    // cout << "AS: STATE: " << n << " exp( -1.0 * i * as_renorm )    : " << exp( -1.0 * complex<double>(0.0,1.0) * as_renorm ) << endl;

    // AS: that did not work for very small G=0
    // for ( int i = 0; i < basis_->localsize(); i++ )
    // {
    //   p[i] *= exp( -1.0 * complex<double>(0.0,1.0) * as_renorm );
    // }

    for ( int i = 0; i < basis_->localsize(); i++ )
    {
      complex<double> c_g = p[i];
      complex<double> c_minusg = p[basis_->index_of_minus_g(i)];

      p[i] = 0.5*(c_g + conj(c_minusg));
      p[basis_->index_of_minus_g(i)] = 0.5*(c_minusg + conj(c_g));
    }

    p[basis_->zero_index()] = complex<double>(real(p[basis_->zero_index()]),0.0);
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::cleanup(void) {
  // set Im( c(G=0) ) to zero and 
  // set the empty rows of the matrix c_ to zero
  // The empty rows are located between i = basis_->localsize() and 
  // c_.mloc(). Empty rows are necessary to insure that the 
  // local size c_.mloc() is the same on all processes, while the 
  // local basis size is not.
  for ( int n = 0; n < c_.nloc(); n++ ) {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    // reset imaginary part of G=0 component to zero
    if ( ctxt_.myrow() == 0 ) {
      // index of G=0 element
      int izero;
      if ( basis_->real() )
        izero = 0;
      else
        izero = basis_->rod_size(0)/2;
      //cout << " izero = " << izero << " G = " << basis_->kv(3*izero) << " " 
      //     << basis_->kv(3*izero+1) << " " << basis_->kv(3*izero+2) << endl;
      p[izero] = complex<double> ( p[izero].real(), 0.0);
    }
    // reset values of empty rows of c_ to zero
    for ( int i = basis_->localsize(); i < c_.mloc(); i++ ) {
      p[i] = 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
SlaterDet& SlaterDet::operator=(SlaterDet& rhs) {
  if ( this == &rhs ) return *this;
  assert(ctxt_.ictxt() == rhs.context().ictxt());
  c_ = rhs.c_;
  if (rhs.ultrasoft_)
     spsi_ = rhs.spsi_;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::memsize(void) const {
  return basis_->memsize() + c_.memsize();
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::localmemsize(void) const {
  return basis_->localmemsize() + c_.localmemsize();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::print(ostream& os, string encoding) {
  FourierTransform ft(*basis_,basis_->np(0),basis_->np(1),basis_->np(2));
  vector<complex<double> > wftmp(ft.np012loc());
  vector<double> wftmpr(ft.np012());
  Base64Transcoder xcdr;
  
  if ( ctxt_.oncoutpe() ) {

    //ewd: SlaterDet does not have weight of each k-point
    const double weight = 1.0; //!! fixed determinant weight to 1.0
    //!! no spin attribute written
    os << "<slater_determinant kpoint=\"" << basis_->kpoint() << "\""
       << "  weight=\"" << weight << "\""
       << " size=\"" << nst() << "\">" << endl;
 
    os << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">" 
       << endl;
    os.setf(ios::fixed,ios::floatfield);
    os.setf(ios::right,ios::adjustfield);
    for ( int i = 0; i < nst(); i++ ) {
      os << " " << setprecision(8) << occ_[i];
      if ( i%10 == 9 )
        os << endl;
    }
    if ( nst()%10 != 0 )
      os << endl;
    os << "</density_matrix>" << endl;
  }
  
  for ( int n = 0; n < nst(); n++ ) {
    // Barrier to limit the number of messages sent to task 0 
    // that don't have a receive posted
    ctxt_.barrier();
    
    // check if state n resides on mype
    if ( c_.pc(n) == ctxt_.mycol() ) {
      //cout << " state " << n << " is stored on column " 
      //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
      int nloc = c_.y(n); // local index
      ft.backward(c_.cvalptr(c_.mloc()*nloc),&wftmp[0]);
      
      double *a = (double*) &wftmp[0];
      for ( int i = 0; i < ft.np012loc(); i++ )
        wftmpr[i] = a[2*i];
        
      for ( int i = 0; i < ctxt_.nprow(); i++ ) {
        if ( i == ctxt_.myrow() ) {
          int size = ft.np012loc();
          //cout << " process " << ctxt_.mype() << " sending block " << i
          //     << " of state "
          //     << n << " to task 0, size = " << size << endl;
          ctxt_.isend(1,1,&size,1,0,0);
          ctxt_.dsend(size,1,&wftmpr[0],1,0,0);
        }
      }
    }
    if ( ctxt_.oncoutpe() ) {
      for ( int i = 0; i < ctxt_.nprow(); i++ ) {
        int size = 0;
        ctxt_.irecv(1,1,&size,1,i,c_.pc(n));
        int istart = ft.np0() * ft.np1() * ft.np2_first(i);
        //cout << " task 0 (" << ctxt_.mype() << ") receiving block " << i
        //     << " of state "
        //     << n << " size=" << size << " istart=" << istart << endl;
        ctxt_.drecv(size,1,&wftmpr[istart],1,i,c_.pc(n));
      }
      
      // wftmpr is now complete on task 0

      if ( encoding == "base64" ) {
        #if AIX
        xcdr.byteswap_double(ft.np012(),&wftmpr[0]);
        #endif
        int nbytes = ft.np012()*sizeof(double);
        int outlen = xcdr.nchars(nbytes);
        char* b = new char[outlen];
        assert(b!=0);
        xcdr.encode(nbytes,(byte*) &wftmpr[0],b);
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"base64\">" << endl;
        xcdr.print(outlen,(char*) b, os);
        os << "</grid_function>\n";
        delete [] b;
      }
      else {
        // encoding == "text" or unknown encoding
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"text\">" << endl;
        int count = 0;
        for ( int k = 0; k < ft.np2(); k++ )
          for ( int j = 0; j < ft.np1(); j++ )
            for ( int i = 0; i < ft.np0(); i++ ) {
              os << " " << wftmpr[ft.index(i,j,k)];
              if ( count++%4 == 3)
                os << "\n";
            }
        if ( count%4 != 0 )
          os << "\n";
        os << "</grid_function>\n";
      }

    }
  }
  if ( ctxt_.oncoutpe() )
    os << "</slater_determinant>" << endl;
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::write(SharedFilePtr& sfp, string encoding, double weight, int ispin,
                      int nspin) const {
  FourierTransform ft(*basis_,basis_->np(0),basis_->np(1),basis_->np(2));
  vector<complex<double> > wftmp(ft.np012loc());
  const bool real_basis = basis_->real();
  const int wftmpr_loc_size = real_basis ? ft.np012loc() : 2*ft.np012loc();
  vector<double> wftmpr(wftmpr_loc_size);
  Base64Transcoder xcdr;

  char* wbuf = 0;
  size_t wbufsize = 0;

  // Segment n on process iprow is sent to row (n*nprow+iprow)/(nprow)
  const Context& colctxt = basis_->context();
  const int nprow = ctxt_.nprow();
  vector<int> scounts(nprow), sdispl(nprow), rcounts(nprow), rdispl(nprow);
  string header;

  //if ( ctxt_.oncoutpe() )
  if ( ctxt_.myproc() == 0 )
    {
      ostringstream ostr_hdr;
      string spin = (ispin > 0) ? "down" : "up";
      ostr_hdr << "<slater_determinant";
      if ( nspin == 2 )
        ostr_hdr << " spin=\"" << spin << "\"";
      ostr_hdr << " kpoint=\"" << basis_->kpoint() << "\"\n"
               << "  weight=\"" << setprecision(12) << weight << "\""
               << " size=\"" << nst() << "\">" << endl;

      ostr_hdr << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">"
               << endl;
      ostr_hdr.setf(ios::fixed,ios::floatfield);
      ostr_hdr.setf(ios::right,ios::adjustfield);
      for ( int i = 0; i < nst(); i++ )
        {
          ostr_hdr << " " << setprecision(8) << occ_[i];
          if ( i%10 == 9 )
            ostr_hdr << endl;
        }
      if ( nst()%10 != 0 )
        ostr_hdr << endl;
      ostr_hdr << "</density_matrix>" << endl;
      header = ostr_hdr.str();
    }

  // serialize all local columns of c and store in segments seg[n]
  string seg;
  for ( int n = 0; n < nstloc(); n++ )
    {
      seg.clear();
      if ( n == 0 && ctxt_.myrow() == 0 )
        seg = header;

      ostringstream ostr;
      //cout << " state " << n << " is stored on column "
      //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
      ft.backward(c_.cvalptr(c_.mloc()*n),&wftmp[0]);

      if ( real_basis )
        {
          double *a = (double*) &wftmp[0];
          for ( int i = 0; i < ft.np012loc(); i++ )
            wftmpr[i] = a[2*i];
        }
      else
        {
          memcpy((void*)&wftmpr[0],(void*)&wftmp[0],
                 ft.np012loc()*sizeof(complex<double>));
        }

      // find index of last process holding some data
      int lastproc = ctxt_.nprow()-1;
      while ( lastproc >= 0 && ft.np2_loc(lastproc) == 0 ) lastproc--;
      assert(lastproc>=0);

      // Adjust number of values on each task to have a number of values
      // divisible by three. This is necessary in order to have base64
      // encoding without trailing '=' characters.
      // The last node in the process column may have a number of values
      // not divisible by 3.

      // data now resides in wftmpr, distributed on ctxt_.mycol()
      // All nodes in the process column except the last have the
      // same wftmpr_loc_size
      // Use group-of-three redistribution algorithm to make all sizes
      // multiples of 3. In the group-of-three algorithm, nodes are divided
      // into groups of three nodes. In each group, the left and right members
      // send 1 or 2 values to the center member so that all three members
      // end up with a number of values divisible by three.

      // Determine how many values must be sent to the center-of-three node
      int ndiff;
      const int myrow = ctxt_.myrow();
      if ( myrow == 0 )
        {
          ndiff = wftmpr_loc_size % 3;
          ctxt_.ibcast_send('c',1,1,&ndiff,1);
        }
      else
        {
          ctxt_.ibcast_recv('c',1,1,&ndiff,1,0,ctxt_.mycol());
        }
      // assume that all nodes have at least ndiff values
      if ( myrow <= lastproc ) assert(wftmpr_loc_size >= ndiff);

      // Compute number of values to be sent to neighbors
      int nsend_left=0, nsend_right=0, nrecv_left=0, nrecv_right=0;
      if ( myrow % 3 == 0 )
        {
          // mype is the left member of a group of three
          // send ndiff values to the right if not on the last node
          if ( myrow < lastproc )
            nsend_right = ndiff;
        }
      else if ( myrow % 3 == 1 )
        {
          // mype is the center member of a group of three
          if ( myrow <= lastproc )
            nrecv_left = ndiff;
          if ( myrow <= lastproc-1 )
            nrecv_right = ndiff;
        }
      else if ( myrow % 3 == 2 )
        {
          // mype is the right member of a group of three
          // send ndiff values to the left if not on the first or last node
          if ( myrow <= lastproc && myrow > 0 )
            nsend_left = ndiff;
        }

      double rbuf_left[2], rbuf_right[2], sbuf_left[2], sbuf_right[2];
      int tmpr_size = wftmpr_loc_size;
      if ( nsend_left > 0 )
        {
          for ( int i = 0; i < ndiff; i++ )
            sbuf_left[i] = wftmpr[i];
          ctxt_.dsend(ndiff,1,sbuf_left,ndiff,ctxt_.myrow()-1,ctxt_.mycol());
          tmpr_size -= ndiff;
        }
      if ( nsend_right > 0 )
        {
          for ( int i = 0; i < ndiff; i++)
            sbuf_right[i] = wftmpr[wftmpr_loc_size-ndiff+i];
          ctxt_.dsend(ndiff,1,sbuf_right,ndiff,ctxt_.myrow()+1,ctxt_.mycol());
          tmpr_size -= ndiff;
        }
      if ( nrecv_left > 0 )
        {
          ctxt_.drecv(ndiff,1,rbuf_left,ndiff,ctxt_.myrow()-1,ctxt_.mycol());
          tmpr_size += ndiff;
        }
      if ( nrecv_right > 0 )
        {
          ctxt_.drecv(ndiff,1,rbuf_right,ndiff,ctxt_.myrow()+1,ctxt_.mycol());
          tmpr_size += ndiff;
        }

      // check that size is a multiple of 3 (except on last node)
      // cout << ctxt_.mype() << ": tmpr_size: " << tmpr_size << endl;
      if ( ctxt_.myrow() != lastproc )
        assert(tmpr_size%3 == 0);
      vector<double> tmpr(tmpr_size);

      // Note: all nodes either receive data or send data, not both
      if ( nrecv_left > 0 || nrecv_right > 0 )
        {
          // this node is a receiver
          int index = 0;
          if ( nrecv_left > 0 )
            {
              for ( int i = 0; i < ndiff; i++ )
                tmpr[index++] = rbuf_left[i];
            }
          for ( int i = 0; i < wftmpr_loc_size; i++ )
            tmpr[index++] = wftmpr[i];
          if ( nrecv_right > 0 )
            {
              for ( int i = 0; i < ndiff; i++ )
                tmpr[index++] = rbuf_right[i];
            }
          assert(index==tmpr_size);
        }
      else if ( nsend_left > 0 || nsend_right > 0 )
        {
          // this node is a sender
          int index = 0;
          int istart=0, iend=wftmpr_loc_size;
          if ( nsend_left > 0 )
            istart = ndiff;
          if ( nsend_right > 0 )
            iend = wftmpr_loc_size - ndiff;
          for ( int i = istart; i < iend; i++ )
            tmpr[index++] = wftmpr[i];
          assert(index==tmpr_size);
        }
      else
        {
          // no send and no recv
          for ( int i = 0; i < wftmpr_loc_size; i++ )
            tmpr[i] = wftmpr[i];
          assert(tmpr_size==wftmpr_loc_size);
        }

      // All nodes (except the last) now have a number of values
      // divisible by 3 in tmpr[]

      if ( ctxt_.myrow()!=lastproc ) assert(tmpr_size%3==0);

      // convert local data to base64 and write to outfile

      // tmpr contains either a real or a complex array

      const string element_type = real_basis ? "double" : "complex";

      if ( encoding == "base64" )
        {
      #if PLT_BIG_ENDIAN
          xcdr.byteswap_double(tmpr_size,&tmpr[0]);
      #endif
          int nbytes = tmpr_size*sizeof(double);
          int outlen = xcdr.nchars(nbytes);
          char* b = new char[outlen];
          assert(b!=0);
          xcdr.encode(nbytes,(byte*) &tmpr[0],b);
          // Note: optional x0,y0,z0 attributes not used, default is zero
          if ( ctxt_.myrow() == 0 )
            {
              // if on first row, write grid function header
              ostr << "<grid_function type=\"" << element_type << "\""
                   << " nx=\"" << ft.np0()
                   << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
                   << " encoding=\"base64\">" << endl;
            }
          xcdr.print(outlen,(char*) b, ostr);
          if ( ctxt_.myrow() == lastproc )
            ostr << "</grid_function>\n";
          delete [] b;
        }
      else
        {
          // encoding == "text" or unknown encoding
          // Note: optional x0,y0,z0 attributes not used, default is zero
          if ( ctxt_.myrow() == 0 )
            {
              // if on first row, write grid function header
              ostr << "<grid_function type=\"" << element_type << "\""
                   << " nx=\"" << ft.np0()
                   << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
                   << " encoding=\"text\">" << endl;
            }
          int count = 0;
          for ( int k = 0; k < ft.np2(); k++ )
            for ( int j = 0; j < ft.np1(); j++ )
              for ( int i = 0; i < ft.np0(); i++ )
                {
                  int index = ft.index(i,j,k);
                  if ( real_basis )
                    ostr << " " << tmpr[index];
                  else
                    ostr << " " << tmpr[2*index] << " " << tmpr[2*index+1];
                  if ( count++%4 == 3)
                    ostr << "\n";
                }
          if ( count%4 != 0 )
            ostr << "\n";
          if ( ctxt_.myrow() == lastproc )
            ostr << "</grid_function>\n";
        }
      // copy contents of ostr stringstream to segment
      seg += ostr.str();
      // cout << ctxt_.mype() << ": segment " << n << " size: " << seg.size()
      //      << endl;

      // seg is defined

      // redistribute segments to tasks within each process column

      for ( int i = 0; i < nprow; i++ )
        {
          scounts[i] = 0;
          sdispl[i] = 0;
          rcounts[i] = 0;
          rdispl[i] = 0;
        }

      int idest = (n*nprow+ctxt_.myrow())/nstloc();
      scounts[idest] = seg.size();

      // send sendcounts to all procs
      MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,colctxt.comm());

      // dimension receive buffer
      int rbufsize = rcounts[0];
      rdispl[0] = 0;
      for ( int i = 1; i < ctxt_.nprow(); i++ )
        {
          rbufsize += rcounts[i];
          rdispl[i] = rdispl[i-1] + rcounts[i-1];
        }
      char* rbuf = new char[rbufsize];

      int err = MPI_Alltoallv((void*)seg.data(),&scounts[0],&sdispl[0],
                              MPI_CHAR,rbuf,&rcounts[0],&rdispl[0],MPI_CHAR,colctxt.comm());

      if ( err != 0 )
        cout << ctxt_.mype()
             << " SlaterDet::write: error in MPI_Alltoallv" << endl;

      if ( rbufsize > 0 )
        {
          // append rbuf to wbuf
          char* tmp;
      try
        {
          tmp = new char[wbufsize+rbufsize];
        }
      catch ( bad_alloc )
        {
          cout << ctxt_.mype() << " bad_alloc in wbuf append "
               << " n=" << n
               << " rbufsize=" << rbufsize
               << " wbufsize=" << wbufsize << endl;
        }
      memcpy(tmp,wbuf,wbufsize);
      memcpy(tmp+wbufsize,rbuf,rbufsize);
      delete [] wbuf;
      wbuf = tmp;
      wbufsize += rbufsize;
        }
      delete [] rbuf;
    }
  // wbuf now contains the data to be written in the correct order

  ctxt_.barrier();

  // compute offsets
  //ewd DEBUG
  //sfp.sync();
  // do a sync only on pes in this SlaterDet
  long long int s_off,tmpoff_;
  s_off = sfp.offset();
  MPI_Allreduce(&s_off,&tmpoff_,1,MPI_LONG_LONG,MPI_MAX,ctxt_.comm());
  sfp.set_offset(tmpoff_);

  MPI_Offset off;
  long long int local_offset,current_offset;
  current_offset = sfp.offset();

  // compute local offset of next write
  long long int local_size = wbufsize;
  MPI_Scan(&local_size, &local_offset, 1,
           MPI_LONG_LONG, MPI_SUM, ctxt_.comm());
  // add base and correct for inclusive scan by subtracting local_size
  local_offset += current_offset - local_size;
  off = local_offset;

  MPI_Status status;

  // write wbuf from all tasks using computed offset
  int len = wbufsize;
  //ewdint err = MPI_File_write_at_all(sfp.file(),off,(void*)wbuf,len,
  //ewd                                MPI_CHAR,&status);
  //ewd:  comment out collective blocking write_at_all since SlaterDet
  //ewd:  can exist on subset of all pes when nparallelkpts > 1
  int err = MPI_File_write_at(sfp.file(),off,(void*)wbuf,len,
                                  MPI_CHAR,&status);
  if ( err != 0 )
    cout << ctxt_.mype()
         << " error in MPI_File_write_at_all" << endl;
  sfp.set_offset(local_offset+len);

  //ewd DEBUG
  //sfp.sync();
  s_off = sfp.offset();
  tmpoff_ = sfp.offset();
  MPI_Allreduce(&s_off,&tmpoff_,1,MPI_LONG_LONG,MPI_MAX,ctxt_.comm());
  sfp.set_offset(tmpoff_);

  delete [] wbuf;

  //if ( ctxt_.oncoutpe() )
  if ( ctxt_.myproc() == 0 )
    {
      string s("</slater_determinant>\n");
      int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*) s.data(),
                                  s.size(),MPI_CHAR,&status);
      if ( err != 0 )
        cout << ctxt_.mype()
             << " error in MPI_File_write, slater_determinant trailer"
             << endl;
      sfp.advance(s.size());
    }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::info(ostream& os) {  

  //if ( ctxt_.myproc()==0 ) {

    //ewd: SlaterDet does not have weight of each k-point
    const double weight = 1.0; //!! fixed determinant weight to 1.0
    //!! no spin attribute written
    os << "<slater_determinant kpoint=\"" << basis_->kpoint() << "\""
       << "  weight=\"" << weight << "\""
       << " size=\"" << nst() << "\">" << endl;
    os << " <!-- sdcontext: " << ctxt_.nprow() << "x" << ctxt_.npcol() << " -->"
       << endl;
    //os << " <!-- sdcontext: " << ctxt_ << " -->" << endl;
    os << "<grid nx=\"" << basis_->np(0) << "\""
       <<      " ny=\"" << basis_->np(1) << "\""
       <<      " nz=\"" << basis_->np(2) << "\"/>" << endl;
    os << " <!-- basis size: " << basis_->size() << " -->" << endl;
    os << " <!-- c dimensions: "
       << c_.m() << "x" << c_.n()
       << "   (" << c_.mb() << "x" << c_.nb() << " blocks)" << " -->" << endl;
    os << " <density_matrix form=\"diagonal\" size=\"" << nst() << "\">" 
       << endl;
    os << " </density_matrix>" << endl;
    os << "</slater_determinant>" << endl;

    //}
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, SlaterDet& sd) {
  sd.print(os,"text");
  return os;
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::print_memory(ostream& os, int kmult, int kmultloc, double& totsum, double& locsum) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setprecision(3);
  PrintMem pm;
    
  const int ngw = basis_->size();
  const int ngwl = basis_->localsize();
  
  double zmult = (double)kmult*sizeof(complex<double>);
  double zmultloc = (double)kmultloc*sizeof(complex<double>);
  double psi_size = (double)ngw*(double)nst()*zmult;
  double psi_locsize = (double)ngwl*(double)nstloc()*zmultloc;

  double betag_size = 0.0;
  double betag_locsize = 0.0;
  double betapsi_size = 0.0;
  double betapsi_locsize = 0.0;
  if (ultrasoft_) {
    int nsp = betag_.size();
    for (int is=0; is<nsp; is++) {
      betag_size += kmult*betag_[is]->memsize();
      betag_locsize += kmultloc*betag_[is]->localmemsize();
      betapsi_size += kmult*betapsi_[is]->memsize();
      betapsi_locsize += kmultloc*betapsi_[is]->localmemsize();
    }
  }
  
  if (ultrasoft_) {
    totsum += 3.*psi_size + betag_size + betapsi_size;
    locsum += 3.*psi_locsize + betag_locsize + betapsi_locsize;
  }
  else {
    totsum += 2.*psi_size;
    locsum += 2.*psi_locsize;
  }
  
  string psi_unit = pm.memunit(psi_size);
  string psi_locunit = pm.memunit(psi_locsize);
  os << "<!-- memory sd.psi      :  " << setw(7) << psi_size << psi_unit << "  (" << psi_locsize << psi_locunit << " local) -->" << endl;
  os << "<!-- memory sd.hpsi     :  " << setw(7) << psi_size << psi_unit << "  (" << psi_locsize << psi_locunit << " local) -->" << endl;
  if (ultrasoft_) {
    os << "<!-- memory sd.spsi     :  " << setw(7) << psi_size << psi_unit << "  (" << psi_locsize << psi_locunit << " local) -->" << endl;
    string betag_unit = pm.memunit(betag_size);
    string betag_locunit = pm.memunit(betag_locsize);
    os << "<!-- memory sd.betag    :  " << setw(7) << betag_size << betag_unit << "  (" << betag_locsize << betag_locunit << " local) -->" << endl;
    string betapsi_unit = pm.memunit(betapsi_size);
    string betapsi_locunit = pm.memunit(betapsi_locsize);
    os << "<!-- memory sd.betapsi  :  " << setw(7) << betapsi_size << betapsi_unit << "  (" << betapsi_locsize << betapsi_locunit << " local) -->" << endl;
  }
}
