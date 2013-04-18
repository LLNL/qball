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
// compile using icc -o testdl testdl.C -ldl
#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;
#include <dlfcn.h>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
using namespace xercesc;

void (*XMLPlatformUtils_Initialize_ptr)(char const*);

int main(int argc, char **argv) {
    void *handle;
    double (*cosine)(double);
    char *error;
    char* libname = "libxerces-c.so.22";

    //handle = dlopen ("/lib/i686/libm.so.6", RTLD_LAZY);
    handle = dlopen (libname, RTLD_LAZY);
    if (!handle) {
        fputs (dlerror(), stderr);
        exit(1);
    }
    char* symbol="ZN11xercesc_2_216XMLPlatformUtils10InitializeEPKc";
    XMLPlatformUtils_Initialize_ptr = 
    (void (*)(char const*)) dlsym(handle,symbol);
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }
    (*XMLPlatformUtils_Initialize_ptr)("");
    
#if 0
    cosine = (double (*)(double)) dlsym(handle, "cos");
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }

    printf ("%f\n", (*cosine)(2.0));
#endif
  
    dlclose(handle);
}
