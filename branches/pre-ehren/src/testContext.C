////////////////////////////////////////////////////////////////////////////////
//
// testContext.c
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testContext.C,v 1.4 2007/04/30 16:15:10 draeger1 Exp $

#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#ifdef USE_MPI  
#include <mpi.h>
#endif

#include "Context.h"

int main(int argc, char **argv)
{
  int mype;
  int npes;
#ifdef USE_MPI  

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  int nr = atoi(argv[1]);
  int nc = atoi(argv[2]);
  
  { // start Context scope
  
    MPI_Comm tcomm_;
    cout << "mype = " << mype << ", tcomm_ = " << tcomm_ << endl;
    cout << "MPI_COMM_NULL = " << MPI_COMM_NULL << endl;

    Context ctxt;
    //cout << mype << ":" << ctxt.mype() << ":" << ctxt.myproc() 
    //     << " base: " << ctxt;


    if (mype == 0)
      cout << "Creating spincontext" << endl;
    Context* spincontext = new Context(nr,nc);


    //Context* subctxt_ = new Context(*spincontext,nr,nc/2,0,nc/2);
    //Context* subctxt_ = new Context(*spincontext,nr,nc/2,0,0);
    Context* subctxt_ = new Context(*spincontext,nr,nc,0,0);
    cout << "Created subctxt on pe = " << mype << endl;

    if (subctxt_->active()) {
      cout << "subctxt active on pe = " << mype << endl;

      Context* my_col_ctxt = 0;
      for ( int icol = 0; icol < subctxt_->npcol(); icol++ ) {

        cout << "creating single-column context for column " << icol << " on pe = " << mype << ", subctxt->nprow = " << subctxt_->nprow() << endl;

        Context* col_ctxt = new Context(*subctxt_,subctxt_->nprow(),1,0,icol);

        cout << "returned from single-column context constructor for column " << icol << " on pe = " << mype << endl;

        subctxt_->barrier();
        if ( icol == subctxt_->mycol() )
          my_col_ctxt = col_ctxt;
        else
          delete col_ctxt;
      }
      cout << "finished creating single-column contexts" << endl;

    }

    cout << "All done, exiting." << endl;

#if 0
    // create subcontext from spincontext
    if (mype == 0)
      cout << "Creating subcontext" << endl;
    Context* subcontext = new Context(*spincontext,nr,nc/2,0,nc/2);

    //if (spincontext->mycol() >= nc/2) {

    if (spincontext->mype() == npes/2) 
      cout << "Start creating col_ctxt" << endl;

    if (subcontext->active()) {

      if (spincontext->mype() <= npes/2) 
        cout << "WARNING:  subcontext->active on pe <= " << npes/2 << endl;


      // try creating one single-column context
      int icol = 2;
      Context* my_col_ctxt_ = 0;
      Context* col_ctxt = new Context(*subcontext,subcontext->nprow(),1,0,icol);
      if (spincontext->mype() == npes/2) 
        cout << "Finished creating col_ctxt" << endl;
      if ( icol == subcontext->mycol() )
        my_col_ctxt_ = col_ctxt;
      else
        delete col_ctxt;
      if (my_col_ctxt_ != 0) {
        cout << "subcontext pe " << subcontext->mype() << " on column " << subcontext->mycol() << " is col_ctxt pe " << my_col_ctxt_->mype() << endl;
      }
    }
      
    //}
    //else {
    //delete subcontext;
    //}

#endif    
#if 0
    vector<Context*> subcontext;

    int nkpar = 2;
    subcontext.resize(nkpar);

    for (int k=0; k<nkpar; k++)
      subcontext[k] = 0;

    if (spincontext->active()) {
      for (int k=0; k<nkpar; k++) {
        int nsubcol = nc/nkpar;
        if (mype == 0)
          cout << "Creating subcontext" << k << endl;
        subcontext[k] = new Context(*spincontext,nr,nsubcol,0,k*nsubcol);
        if (mype == 0)
          cout << "Done." << endl;

        if (subcontext[k]->active()) {
          // mimic SlaterDet constructor by creating single-column subcontexts
          for ( int icol = 0; icol < subcontext[k]->npcol(); icol++ ) {
            Context* my_col_ctxt_ = 0;
            Context* col_ctxt = new Context(*subcontext[k],subcontext[k]->nprow(),1,0,icol);
            subcontext[k]->barrier();
            if ( icol == subcontext[k]->mycol() )
              my_col_ctxt_ = col_ctxt;
            else
              delete col_ctxt;
          }
        }
      }    
    }
#endif

#if 0
    vector<Context*> c;
    
    c.push_back(new Context(nr,nc));
    c.push_back(new Context(*c[0],2,2,1,1));
    for ( int icol = 0; icol < c[0]->npcol(); icol++ )
    {
      ctxt.barrier();
      c.push_back(new Context(*c[0],c[0]->nprow(),1,0,icol));
    }
#if 0
#endif
    
    for ( int i = 0; i < c.size(); i++ )
    {
      Context* pc = c[i];
      cout << mype << ":" << pc->mype() << ":" << pc->myproc()
           << " at (" << pc->myrow() << "," << pc->mycol() << ")" 
           << " in c" << i << ": " << *pc;
    }
    
#if 0
    MPI_Comm comm = c[1]->comm();
    int mype_in_c1,size_of_c1;
    MPI_Comm_rank(comm,&mype_in_c1);
    MPI_Comm_size(comm,&size_of_c1);
    cout << mype << ": mype_in_c1: " << mype_in_c1 
         << " size_of_c1=" << size_of_c1
         << " comm[c1]=" << comm << endl;
    
    // test dgsum2d function
    double a = c[1]->mype();
    cout << c[1]->mype() << ": a     = " << a << endl;
    c[1]->dsum('R',1,1,&a,1);
    cout << c[1]->mype() << ": a_sum_row = " << a << endl;
    c[1]->dsum('C',1,1,&a,1);
    cout << c[1]->mype() << ": a_sum_all = " << a << endl;
    
#endif
    for ( int i = 0; i < c.size(); i++ )
    {
      delete c[i];
    }
    
//   // test reference counting
//   if ( npes%2 == 0 && npes >= 4 )
//   {
//     Context *c1 = new Context(npes/2,2);
//     cout << "c1: " << *c1 << endl;
//     Context *c2;
//     if ( c1->active() )
//       c2 = new Context(*c1,npes/2,npes/2,1);
//     else
//       c2 = 0;
//     // this line causes crash: Context *c2 = new Context(*c1,1,1,1);
//     delete c1;
//     if ( c2 != 0 ) cout << c2->mype() << " c2: " << *c2;
//     delete c2;
//   }
  
#if 0
  }
#endif
#endif
  } // end Context scope

  MPI_Finalize();

#else

  mype=0;
  npes=1;
  {
    BlacsContext c1(1,1);
    cout << " c1.ictxt = " << c1.ictxt() << endl;
  }

#endif
}
