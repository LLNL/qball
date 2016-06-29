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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
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


