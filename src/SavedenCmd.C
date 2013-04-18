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
// SavedenCmd.C:
//
////////////////////////////////////////////////////////////////////////////////


#include "SavedenCmd.h"
#include "fstream"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SavedenCmd::action(int argc, char **argv) {

  if ( !(argc == 2 || argc == 3) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: saveden filename -->" << endl;
    return 1;
  }

  char* filename = argv[1];
  string format = "gopenmol";
  for ( int i = 1; i < argc; i++ ) {
    string arg(argv[i]);
    if ( arg == "-text" )
      format = "binary";
    else if ( arg=="-molmol" ) 
      format = "molmol";
    else if ( arg=="-gopenmol" )
      format = "gopenmol";
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else {
      if ( ui->oncoutpe() )
        cout << "  <!-- use: saveden [-binary|-molmol|-gopenmol] filename -->" 
             << endl;
      return 1;
    }
  }



  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  ofstream os;
  if (format == "binary") 
    os.open(filename,ofstream::binary); // binary output
  else {
    os.open(filename,ofstream::out);    // text output
    os.setf(ios::fixed,ios::floatfield);
    os << setprecision(8);
  }

  ChargeDensity cd_(*s);
  cd_.update_density();
  const Context* ctxt_ = s->wf.spincontext(0);
  FourierTransform* ft_ = cd_.vft();

  if (s->wf.nspin() == 2 && ctxt_->oncoutpe() )
    cout << "<WARNING> saveden command only prints spin = 0 density </WARNING>" << endl;
  
  //ewd DEBUG:  calculate maximum density
  double maxden = 0.0;

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

        //ewd DEBUG
        for (int j=0; j<size; j++) 
          if (rhortmp[j] > maxden) 
            maxden = rhortmp[j];
        
        // write this portion of the density to file
        if (format == "binary") {
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
        else {
          // headers for visualization formats

          if (i==0) {
            if (format == "molmol") {
              D3vector a0 = s->wf.cell().a(0);
              D3vector a1 = s->wf.cell().a(1);
              D3vector a2 = s->wf.cell().a(2);
              if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                  a2.x != 0.0 || a2.y != 0.0 )
                os << "Error writing header:  Molmol isosurface requires rectangular box!" << endl;
              double dx = a0.x/(double)ft_->np0();
              double dy = a1.y/(double)ft_->np1();
              double dz = a2.z/(double)ft_->np2();
              //    origin    npoints       grid spacing
              os << "0.0 " << ft_->np0() << " " << dx << endl;
              os << "0.0 " << ft_->np1() << " " << dy << endl;
              os << "0.0 " << ft_->np2() << " " << dz << endl;
            }
            else if (format == "gopenmol") {
              D3vector a0 = s->wf.cell().a(0);
              D3vector a1 = s->wf.cell().a(1);
              D3vector a2 = s->wf.cell().a(2);
              if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                  a2.x != 0.0 || a2.y != 0.0 )
                os << "Error writing header:  gOpenMol isosurface requires rectangular box!" << endl;
              os << "3 200" << endl;  //ewd copying this from another utility, not sure what it means
              os << ft_->np0() << " " << ft_->np1() << " " << ft_->np2() << endl;
              os << "0.0 " << a2.z << endl;
              os << "0.0 " << a1.y << endl;
              os << "0.0 " << a0.x << endl;
            }
          }
          ostringstream oss;
          for (int j=0; j<size; j++) 
            oss << rhortmp[j] << endl;
          string tos = oss.str();
          os.write(tos.c_str(),tos.length());
        }
      }
    }
  }
  os.close();

  //ewd DEBUG
  if ( ctxt_->oncoutpe() ) 
    cout << "<!-- Saveden: maximum density = " << maxden << " -->" << endl;

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
