#ifndef SPHERICALINTEGRATION_H
#define SPHERICALINTEGRATION_H
#include <vector>
using namespace std;

class SphericalIntegration {

  private:

  public:

  SphericalIntegration(void);
  ~SphericalIntegration(void);

  double clebsch_gordan(int j1, int m1, int j2, int m2, int j3, int m3);
  double clebsch_gordan_real(int j1, int m1, int j2, int m2, int j3, int m3);
  int ifactorial(int v);
  double dfactorial(int v);
  double ylm_real(int l, int m, double gx, double gy, double gz);
  void init_ll74grid(vector<double> &llsph_r, vector<double> &llsph_wt);
  void llsph_gen_oh(int code, double a, double b, double v, vector<double>& r, vector<double>& wt);

};
#endif
