////////////////////////////////////////////////////////////////////////////////
//
// tb64.C: test Xerces Base64 encoding/decoding
//
// icc tb64.C -o tb64 -I$(XMLDIR)/include -L$(XMLDIR)/lib -lxercesc
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <cassert>
using namespace std;
#include <xercesc/util/Base64.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;

int main()
{
  double a[10];
  unsigned int outlen;
  const int n = 32;
  
  for ( int i = 0; i < n; i++ )
    a[i] = 3*i;
    
  //////////////////////////////////////////////////////////////////////////////
  // encode array of double into base64 format
  
  XMLByte* b = Base64::encode((XMLByte*)a, n*sizeof(double), &outlen);
  assert(b!=0);
  cout << " out length: " << outlen << endl;
  cout << "<vector type=\"double\" size=\"" << n << "\" encoding=\"base64\">\n"; 
  cout.write((char*) b, outlen);
  cout << "</vector>\n";
  XMLString::release(&b);

  //////////////////////////////////////////////////////////////////////////////
  // decode base64-encoded character string to an array of double

  string cs = 
  "AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQAAAAAAA\n"
  "ACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkAAAAAAAAA0QAAA\n"
  "AAAAADZAAAAAAAAAOEAAAAAAAAA6QAAAAAAAADxAAAAAAAAAPkAAAAAAAABA\n"
  "QAAAAAAAAEFAAAAAAAAAQkAAAAAAAABDQAAAAAAAAERAAAAAAAAARUAAAAAA\n"
  "AABGQAAAAAAAAEdAAAAAAAAASEAAAAAAAABJQAAAAAAAAEpAAAAAAAAAS0AA\n"
  "AAAAAABMQAAAAAAAAE1AAAAAAAAATkAAAAAAAABPQA==\n";

  b = Base64::decode((XMLByte*)cs.c_str(), &outlen);
  assert(b!=0);
  cout << " out length: " << outlen << endl;
  
  double* d = (double*) b;
  for ( int i = 0; i < outlen/sizeof(double); i++ )
    cout << d[i] << endl;

  XMLString::release(&b);
  
  //////////////////////////////////////////////////////////////////////////////
  // decode text-encoded character string to an array of double
  
  const int m = 7;
  string s = "0.1  0.2  0.3 \n  0.4 .5 0.6  \n 0.7";
  istringstream is(s); 

  double e[m];
  for ( int i = 0; i < m; i++ )
    is >> e[i];

  for ( int i = 0; i < m; i++ )
    cout << e[i] << endl;

  //////////////////////////////////////////////////////////////////////////////
  // test stringstream functionality
  stringstream ss;
  
  // following lines occur in the "character" function of the Handler
  ss << "0.123";
  ss << "0.456";
  
  // following lines occur in the "endElement" function of the Handler
  double f;
  ss >> f;
  cout << f << endl;
  ss >> f;
  cout << f << endl;
  
}
