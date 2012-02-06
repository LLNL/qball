#include "ClebschGordan.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ClebschGordan::ClebschGordan(void) {
}
////////////////////////////////////////////////////////////////////////////////
ClebschGordan::~ClebschGordan(void) {
}
////////////////////////////////////////////////////////////////////////////////
double ClebschGordan::coefficient(int j1, int m1, int j2, int m2, int j3, int m3) {

  if ((m1+m2) != m3) return 0.0;
  if ((j1-m1) < 0) return 0.0;
  if ((j2-m2) < 0) return 0.0;
  if ((j3-m3) < 0) return 0.0;
  if ((j1-j2+j3) < 0) return 0.0;
  if ((j2-j1+j3) < 0) return 0.0;
  if ((j1+j2-j3) < 0) return 0.0;
  
  int rmax = j1+j2-j3;
  if ((j1+m1) > rmax) rmax = j1+m1;
  if ((j1-m2-j3) > rmax) rmax = j1-m2-j3;
  if ((j2+m1-j3) > rmax) rmax = j2+m1-j3;
  if ((j2-m2) > rmax) rmax = j2-m2;

  //cout << "rmax = " << rmax << endl;
  
  double pref = (double)((2.*j3+1.)*dfactorial(j1-j2+j3)*dfactorial(j2-j1+j3)*dfactorial(j1+j2-j3))/(double)dfactorial(j1+j2+j3+1);

  double num = (double)(dfactorial(j3+m3)*dfactorial(j3-m3)*dfactorial(j2+m2)*dfactorial(j2-m2)*dfactorial(j1+m1)*dfactorial(j1-m1));

  double sum = 0.0;
  for (int r=0; r<=rmax; r++) {
    double denom = (double)(dfactorial(r)*dfactorial(j1+j2-j3-r)*dfactorial(j1+m1-r)*dfactorial(j3-j1+m2+r)*dfactorial(j3-m1-j2+r)*dfactorial(j2-m2-r));
    if (denom > 0.0) {
      double neg = pow(-1.0,r+j1+j2-j3);
      sum += (double)neg*sqrt(pref)*sqrt(num)/denom;
    }
  }
  return sum;
}
////////////////////////////////////////////////////////////////////////////////
int ClebschGordan::ifactorial(int v) {
  if (v < 0)
    return 0;
  else if (v == 0)
    return 1;
  else {
    int f = 1;
    for (int i=1; i<=v; i++)
      f*=i;
    return f;
  }
}
////////////////////////////////////////////////////////////////////////////////
double ClebschGordan::dfactorial(int v) {
  int f = ifactorial(v);
  return (double)f;
}
