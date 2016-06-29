////////////////////////////////////////////////////////////////////////////////
//
// testContext.c
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testContext.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

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
  
    Context ctxt;
    cout << mype << ":" << ctxt.mype() << ":" << ctxt.myproc() 
         << " base: " << ctxt;

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
