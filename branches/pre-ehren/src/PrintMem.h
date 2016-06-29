////////////////////////////////////////////////////////////////////////////////
//
// PrintMem.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PrintMem.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef PRINTMEM_H
#define PRINTMEM_H

#include <string>
using namespace std;

class PrintMem
{
  public:

  string memunit(double& v);
  PrintMem(void);
  ~PrintMem(void);
};
#endif
