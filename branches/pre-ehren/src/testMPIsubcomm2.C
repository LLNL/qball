////////////////////////////////////////////////////////////////////////////////
//
// testMPIsubcomm2
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testMPIsubcomm2.C,v 1.1 2007/04/30 16:15:10 draeger1 Exp $

#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#ifdef USE_MPI  
#include <mpi.h>
#endif

int main(int argc, char **argv) {
  int mype;
  int npes;
#ifdef USE_MPI  

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  int nr = atoi(argv[1]);
  int nc = atoi(argv[2]);
  
  //ewd create subcommunicator manually to see if we get same error on uP
  {
    MPI_Comm comm_;
    vector<int> pmap_;
    pmap_.resize(npes);
    for ( int i = 0; i < npes; i++ )
      pmap_[i] = i;

    MPI_Group group_world, subgroup;
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);
    MPI_Group_incl(group_world,npes,&pmap_[0],&subgroup);
    MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm_);
    MPI_Group_free(&group_world);
    MPI_Group_free(&subgroup);

    // now create subcommunicator of this

    vector<int> pmap2_;
    pmap2_.resize(npes);
    // build pmap2
    int i = 0;
    int icstart = 0;
    int irstart = 0;
    for ( int ic = icstart; ic < icstart+0.5*nc; ic++ )
      for ( int ir = irstart; ir < irstart+nr; ir++ ) {
        pmap2_[i] = pmap_[i]-pmap_[0];
        i++;
      }

    MPI_Comm comm2_;
    MPI_Group c_group, subgroup2;
    MPI_Comm_group(comm_,&c_group);
    MPI_Group_incl(c_group,0.5*npes,&pmap2_[0],&subgroup2);
    MPI_Comm_create(comm_,subgroup2,&comm2_);
    MPI_Group_free(&c_group);
    MPI_Group_free(&subgroup2);

    cout << "mype = " << mype << ", communicator handle = " << comm_ << ", subcommunicator handle = " << comm2_ << " (MPI_COMM_NULL = " << MPI_COMM_NULL << ")" << endl;
    
    MPI_Barrier(MPI_COMM_WORLD);

    //if (mype == 0)
    //  cout << "now call MPI_Comm_free..." << endl;
    MPI_Comm_free(&comm_);

    //comment out this if to cause a crash on uP
    //if (mype < 0.5*npes)
      MPI_Comm_free(&comm2_);

  }
  MPI_Finalize();

#endif

}

