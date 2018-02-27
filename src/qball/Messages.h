////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2016, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Xavier Andrade <xavier@llnl.gov>
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
// Messages.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef MESSAGES_H
#define MESSAGES_H

#include <math/d3vector.h>
#include <string>
using namespace std;

class Messages
{

  public:
  
  static void fatal(const string & info){
    cerr << endl;
    cerr << "**********************************************************************************" << endl;
    cerr << endl;
    cerr << "    Error: " << info << "." << endl;
    cerr << endl;
    cerr << "**********************************************************************************" << endl;
    cerr << endl;
    exit(1);
  }
  
};

#endif

