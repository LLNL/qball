////////////////////////////////////////////////////////////////////////////////
//
// testBase64Transcoder.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Base64Transcoder.h"
#include <iostream>
using namespace std;

int main()
{
  const int n = 7;
  double a[n] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  double c[n] = { 0, 0, 0, 0, 0, 0, 0 };
  
  Base64Transcoder xcdr;
  
  int nbytes = n * sizeof(double);
  
  int nchars = xcdr.nchars(nbytes);
  
  cout << " nbytes=" << nbytes << endl;
  cout << " nchars=" << nchars << endl;
  
  char* b = new char[nchars];
  
  xcdr.encode(nbytes,(unsigned char*) &a[0],b);
  
  cout << " b=" << b << endl;
  
  xcdr.decode(nchars,b,(unsigned char*) &c[0]);
  
  for ( int i = 0; i < n; i++ )
    assert(a[i]==c[i]);
    
  cout << " done" << endl;
  
  return 0;
}
