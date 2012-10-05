////////////////////////////////////////////////////////////////////////////////
//
// SampleWriter.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleWriter.C,v 1.2 2009/04/05 15:27:36 draeger1 Exp $


#include "SampleWriter.h"
#include "Sample.h"
#include "fstream"
#include "qbox_xmlns.h"
#include "Timer.h"
#include "SharedFilePtr.h"
#include <sstream>
#include <cstring>
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleWriter::SampleWriter(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SampleWriter::writeSample(const Sample& s, const string filename,
                              string description,
                              bool base64, bool atomsonly, bool serial)
{
  Timer tm;
  tm.start();
  // set default encoding
  string encoding =  base64 ? "base64" : "text";
  const char* filename_cstr = filename.c_str();

  long long file_size;

  if ( serial )
  {
    ofstream os;
    if ( ctxt_.oncoutpe() )
    {
      os.open(filename_cstr);
      cout << "  SaveCmd: saving to file " << filename
           << ", encoding=" << encoding << endl;

      os <<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         <<"<fpmd:sample xmlns:fpmd=\""
         << qbox_xmlns()
         << "\"\n"
         <<" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         <<" xsi:schemaLocation=\""
         << qbox_xmlns() << " sample.xsd\">"
         << endl;

      os << "<description> " << description
         << " </description>" << endl;
      os << s.atoms;
    }

    if ( !atomsonly )
    {
      s.wf.print(os,encoding,"wavefunction");
      if ( s.wfv != 0 )
        s.wfv->print(os,encoding,"wavefunction_velocity");
    }

    if ( ctxt_.oncoutpe() )
      os << "</fpmd:sample>" << endl;

    os.close();
  }
  else
  {
    MPI_File fh;
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Offset fsize;
    SharedFilePtr sfp(ctxt_.comm(),fh);

    int err;
    err = MPI_File_open(ctxt_.comm(),(char*) filename_cstr,
                        MPI_MODE_WRONLY|MPI_MODE_CREATE,info,&fh);
    if ( err != 0 )
      cout << s.ctxt_.mype() << ": error in MPI_File_open: " << err << endl;

    MPI_File_set_size(fh,0);
    ctxt_.barrier();

    MPI_Status status;
    if ( ctxt_.oncoutpe() )
    {
      string header("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      "<fpmd:sample xmlns:fpmd=\"");
      header += qbox_xmlns();
      header += string("\"\n");
      header += string(
      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
      header += string(" xsi:schemaLocation=\"");
      header += qbox_xmlns();
      header += string(" sample.xsd\">\n");

      string desc = string("<description> ") +
        description +
        string(" </description>\n");

      header += desc;

      ostringstream ss("");
      ss << s.atoms;
      header += ss.str();
      int len = header.size();
      err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)header.c_str(),
            len,MPI_CHAR,&status);
      if ( err != 0 )
        cout << ctxt_.mype() << ": error in MPI_File_write: header "
             << err << endl;
      sfp.advance(len);
    }
    sfp.sync();

    if ( !atomsonly )
    {
      s.wf.write(sfp,encoding,"wavefunction");
      if ( s.wfv != 0 )
        s.wfv->write(sfp,encoding,"wavefunction_velocity");
    }

    sfp.sync();

    if ( ctxt_.oncoutpe() )
    {
      char *trailer = "</fpmd:sample>\n";
      int len = strlen(trailer);
      err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)trailer,
              len,MPI_CHAR,&status);
      if ( err != 0 )
        cout << ctxt_.mype() << ": error in MPI_File_write: trailer "
             << err << endl;
      sfp.advance(len);
    }

    sfp.sync();

    file_size = sfp.offset();

    err = MPI_File_close(&fh);
    if ( err != 0 )
      cout << ctxt_.mype() << ": error in MPI_File_close: " << err << endl;
  }

  tm.stop();
  if ( ctxt_.oncoutpe() )
  {
    cout << " SampleWriter: write time: "
         << setprecision(3) << tm.real() << " s"
         << endl;
    if ( !serial )
    {
      cout << " SampleWriter: file size: " << file_size << endl;
      cout << " SampleWriter: aggregate write rate: "
           << setprecision(2) << file_size/(tm.real()*1024*1024)
           << " MB/s" << endl;
    }
  }
}
