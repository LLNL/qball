////////////////////////////////////////////////////////////////////////////////
//
// LoaddenCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoaddenCmd.C,v 1.1 2006/12/22 01:17:11 draeger1 Exp $


#include "LoaddenCmd.h"
#include "fstream"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int LoaddenCmd::action(int argc, char **argv) {

  if ( !(argc == 2) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: loadden filename -->" << endl;
    return 1;
  }

  char* filename = argv[1];

  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  bool loadtext = false;

  ofstream os;
  if (loadtext) {
    os.open(filename,ofstream::out);    // text output
    os.setf(ios::fixed,ios::floatfield);
    os << setprecision(8);
  }
  else {
    os.open(filename,ofstream::binary); // binary output
  }

  ChargeDensity cd_(s);
  cd_.update_density();
  const Context* ctxt_ = s->wf.spincontext(0);
  FourierTransform* ft_ = cd_.vft();

  if (ctxt_->mycol() == 0) {

    vector<double> rhortmp(ft_->np012loc());
    for (int j = 0; j < ft_->np012loc(); j++)
      rhortmp[j] = cd_.rhor[0][j];

    for ( int i = 0; i < ctxt_->nprow(); i++ ) {
      if ( i == ctxt_->myrow() ) {
        int size = ft_->np012loc();
        //cout << " process " << ctxt_->mype() << " sending block " << i
        //     << " of density to task 0, size = " << size << endl;
        ctxt_->isend(1,1,&size,1,0,0);
        ctxt_->dsend(size,1,&rhortmp[0],1,0,0);
      }
    }
    if ( ctxt_->oncoutpe() ) {
      for ( int i = 0; i < ctxt_->nprow(); i++ ) {
        int size = 0;
        ctxt_->irecv(1,1,&size,1,i,0);
        //int istart = cd_.vft.np0() * cd_.vft.np1() * cd_.vft.np2_first(i);
        //cout << " process " << ctxt_->mype() << " receiving block " << i
        //     << " of density on task 0, size = " << size << endl;
        ctxt_->drecv(size,1,&rhortmp[0],1,i,0);

        // write this portion of the density to file
        if (loadtext) {
          if (i==0)
            os << "  " << ft_->np0() << "  " << ft_->np1() << "  " << ft_->np2() << endl;
          for (int j=0; j<size; j++) 
            os << rhortmp[j] << endl;
        }
        else {
          if (i==0) {
            int np0 = ft_->np0();
            os.write((char*)&np0,sizeof(int));
            int np1 = ft_->np1();
            os.write((char*)&np1,sizeof(int));
            int np2 = ft_->np2();
            os.write((char*)&np2,sizeof(int));
          }
          os.write((char*)&rhortmp[0],sizeof(double)*size);
        }

      }
    }
  }
  os.close();


  /*
  // print out real-space density
  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
  for ( int i = 0; i < np0v; i++ ) {
    for ( int j = 0; j < np1v; j++ ) {
      for ( int k = 0; k < np2v; k++ ) {
        int pt = k*np0v*np1v + j*np0v + i;
        cout << "DENSITY:  " << i << "  "  << j << "  "  << k << "  " << setprecision(8) << prhor[pt] << endl;
      }
    }  
  }
  */

  return 0;
}
