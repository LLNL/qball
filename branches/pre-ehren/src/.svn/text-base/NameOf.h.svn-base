////////////////////////////////////////////////////////////////////////////////
//
// NameOf.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NameOf.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef NAMEOF_H
#define NAMEOF_H

#include <string>
using namespace std;

// predicate class for searching T* containers by name
// T must be a pointer type to something that has a name() member
template <class T>
class NameOf
{
  public:
  string name;

  NameOf<T>(string s) : name(s) {};
  bool operator() (T t) const
  {
    return t->name() == name;
  }
};
#endif
