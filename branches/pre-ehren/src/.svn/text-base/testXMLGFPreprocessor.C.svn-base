////////////////////////////////////////////////////////////////////////////////
//
// testXMLGFPreprocessor.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testXMLGFPreprocessor.C,v 1.2 2008/04/07 22:00:37 draeger1 Exp $


#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
using namespace std;

#include "Context.h"
#include "Matrix.h"
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
