////////////////////////////////////////////////////////////////////////////////
//
// tenv.C
//
////////////////////////////////////////////////////////////////////////////////

const char* const xmlns_url = "http://www.llnl.gov/casc/fpmd/qbox/1.0";

#include <iostream>
#include <string>
using namespace std;

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;

#include <sys/utsname.h>
#include <unistd.h>
#include <ctime>
#include <cstdlib>
#if AIX || OSF1
#include<filehdr.h>
#endif

string isodate(void)
{
  const time_t t = time(NULL);
  struct tm* tms = gmtime(&t);
  char s[32];
  const char* fmt = "%Y-%m-%dT%TZ";
  strftime(s,32,fmt,tms);
  string st(s);
  return st;
}

int main(int argc, char **argv, char **envp)
{
  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  // Identify executable name, checksum, size and link date
  if ( getlogin() != 0 ) 
    cout << "<user> " << getlogin() << " </user>" << endl;
  cout << "<effective_uid> " << getuid() << "</effective_uid>" << endl;
#if AIX || OSF1
  // read filehdr for link time
  filehdr hdr;
  FILE *fx = fopen(argv[0],"r");
  if ( fx != 0 )
  {
    size_t sz = fread((void*)&hdr,sizeof(filehdr),1,fx);
    fclose(fx);
    string s = ctime((time_t*)&hdr.f_timdat);
    cout << "<linktime> " << s << " </linktime>" << endl;
  }
#endif
  
  // Identify platform
  {
    struct utsname un;
    uname (&un);
    cout << "<sysname> " << un.sysname << " </sysname>" << endl;
    cout << "<nodename> " << un.nodename << " </nodename>" << endl;
  }

  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Never;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  
  // XMLPlatformUtils initialization on task 0 only
  //try
  {
    XMLPlatformUtils::Initialize();
  }

//  catch (const XMLException& toCatch)
//  {
//    cout << "  <!-- Sample::readSample: Error during XML initialization :\n"
//         << StrX(toCatch.getMessage()) << " -->" << endl;
//    ierr = 1;
//  }
  return 0;
}
