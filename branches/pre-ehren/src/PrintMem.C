#include "PrintMem.h"
#include <iostream>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PrintMem::PrintMem() { }
////////////////////////////////////////////////////////////////////////////////
PrintMem::~PrintMem() { }
////////////////////////////////////////////////////////////////////////////////
string PrintMem::memunit(double& v) {
  const double megabyte = 1024.*1024.;
  const double gigabyte = megabyte*1024.;
  const double terabyte = gigabyte*1024.;

  string unit = " MB";
  if (v > terabyte) {
    unit = " TB";
    v /= terabyte;
  }
  else if (v > gigabyte) {
    unit = " GB";
    v /= gigabyte;
  }
  else {
    v /= megabyte;
  }
  
  return unit;
}
