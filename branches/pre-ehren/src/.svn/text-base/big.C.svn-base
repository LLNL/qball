// $Id: big.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
using namespace std;

int main()
{
   unsigned long long kb = 1024;
   unsigned long long mb = kb * kb;
   unsigned long long gb = kb * mb;
   unsigned long long four_gb = gb * 4;
   long long ten_gb = gb * 10;
   int write_count = four_gb / kb;

   cout << "sizeof(unsigned long long) " << sizeof(unsigned long long) << endl;
   cout << "sizeof(std::streamoff)     " << sizeof(std::streamoff) << endl;
   cout << "sizeof(std::streampos)     " << sizeof(std::streampos) << endl;
   cout << "sizeof(std::streamsize)    " << sizeof(std::streamsize) << endl;

   cout << "kb:       " << kb
        << "\nmb:       " << mb
        << "\ngb:       " << gb
        << "\nfour_gb:  " << four_gb
        << "\nten_gb:   " << ten_gb
        << endl;

   unsigned char* buf = new unsigned char[kb];
   memset(buf, 0, kb);

   double total_writes = write_count;
   double writes = 0;
   ofstream os("big_file");

   cout << setiosflags(ios::fixed)
        << setprecision(0);

   for (int i = 0; i < write_count; ++i)
   {
      os.write((char*)buf, kb);

      ++writes;

      cout << "\r"  << setw(3)
           << (writes / total_writes * 100.0) << "%"
           << flush;
   }

   os.close();

   cout << "\r100%\nFinished..." << endl;

   // Open for reading and test seek.
   ifstream is("big_file");
   std::streampos pos = four_gb-1000; // Arbitary position.
   is.seekg(pos);
   std::streampos new_pos = is.tellg();

   if (pos == new_pos)
   {
      cout << "Seek to " << pos << " worked!" << endl;
   }
   else
   {
      cout << "Seek to " << pos << " failed!" << endl;
   }

   is.close();

   return 0;
}


