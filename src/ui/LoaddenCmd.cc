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
// LoaddenCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>


#include "LoaddenCmd.h"
#include "fstream"
#include <qball/isodate.h>
#include <qball/release.h>
#include <qball/qbox_xmlns.h>
#include <qball/ChargeDensity.h>
#include <qball/FourierTransform.h>

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

  ChargeDensity cd_(*s);
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
