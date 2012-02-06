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
