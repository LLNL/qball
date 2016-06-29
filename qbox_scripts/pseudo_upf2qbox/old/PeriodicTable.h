////////////////////////////////////////////////////////////////////////////////
//
// PeriodicTable.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PeriodicTable.h,v 1.1.1.1 2003/11/11 18:59:35 fgygi Exp $

#ifndef PERIODICTABLE_H
#define PERIODICTABLE_H

#include <map>
#include <vector>
#include <string>
using namespace std;

struct Element
{
  int z;
  string symbol;
  string config;
  double mass;
  Element(int zz, string s, string c, double m) : z(zz), symbol(s), config(c),
    mass(m) {}
};

class PeriodicTable
{
  private:

  vector<Element> ptable;
  map<string,int> zmap;

  public:

  PeriodicTable(void);
  int z(string symbol) const;
  string symbol(int zval) const;
  string configuration(int zval) const;
  string configuration(string symbol) const;
  double mass(int zval) const;
  double mass(string symbol) const;
  int size(void) const;

};
#endif
