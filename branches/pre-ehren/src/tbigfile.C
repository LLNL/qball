// 
// Write large files ( > 2 GB ) on Linux
//

#define _LARGEFILE_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <cstdio>
#include <cstdlib>

#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
  FILE *outfile;
  outfile = fopen( "t.dat", "w" );
  
  const off_t n = atoi(argv[1]);
  const off_t m = 1024*1024/8;
  
  double a[m];
  
  cout << " sizeof(int): " << sizeof(int) << endl;
  cout << " sizeof(size_t): " << sizeof(size_t) << endl;
  cout << " sizeof(long): " << sizeof(long) << endl;
  cout << " sizeof(off_t): " << sizeof(off_t) << endl;
  
  for ( int i = 0; i < n; i++ )
    fwrite(&a[0],sizeof(double),m,outfile);
  
  fclose(outfile);
  
  outfile = fopen("t.dat","r");
  off_t offset = 8*m*n-1000;
  fseeko(outfile,offset,SEEK_SET);
  off_t pos = ftello(outfile);
  
  cout << offset << endl;
  cout << pos << endl;
  fclose(outfile);
}

