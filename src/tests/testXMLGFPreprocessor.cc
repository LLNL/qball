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
// testXMLGFPreprocessor.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>


#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
using namespace std;

#include <qball/Context.h>
#include <math/Matrix.h>
#include "XMLGFPreprocessor.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  {
    assert(argc==4);
    const int nr = atoi(argv[1]);
    const int nc = atoi(argv[2]);
    const char* const filename = argv[3];
  
    Context ctxt(nr,nc,'c'); // context on which gfdata is defined
    DoubleMatrix gfdata(ctxt);
    string xmlcontent;
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << gfdata.context().mype() 
         << " on " << processor_name << endl;

    XMLGFPreprocessor xmlgfp;
    xmlgfp.process(filename,gfdata,xmlcontent);
    
#if 0
    // write all gfdata on file
    ofstream gf_file("gf.dat");
    gf_file << gfdata;
#endif
    
    if ( ctxt.oncoutpe() )
    {
      ofstream xmlfile("xml.dat");
      xmlfile << xmlcontent;
    }
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
