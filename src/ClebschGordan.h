#ifndef CLEBSCHGORDAN_H
#define CLEBSCHGORDAN_H
#include <vector>

class ClebschGordan {

  private:

  public:

  ClebschGordan(void);
  ~ClebschGordan(void);

  double coefficient(int j1, int m1, int j2, int m2, int j3, int m3);
  int ifactorial(int v);
  double dfactorial(int v);

};
#endif
