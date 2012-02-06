////////////////////////////////////////////////////////////////////////////////
//
// Hugoniostat.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Hugoniostat.h,v 1.1 2008/05/05 19:37:43 draeger1 Exp $

#ifndef HUGONIOSTAT_H
#define HUGONIOSTAT_H

#include <vector>
using namespace std;

class Hugoniostat {
  private:

  double ref_etot_;
  double ref_pressure_;
  double ref_volume_;
  double th_temp_;
  double deltatemp_;
  double lasthugavg_;
  bool firststep_;

  vector<double> temp_history;
  vector<double> etot_history;
  vector<double> press_history;
  vector<double> vol_history;
  vector<double> hug_history;

  int updatefreq_;
  bool oncoutpe_;
  bool sameSign(double a, double b);

  public:

  bool updatenow;

  Hugoniostat(double etotref, double vref, double pref, double temp, bool oncoutpe);
  ~Hugoniostat();
  void set_updatefreq(int updatefreq);
  void set_deltatemp(double deltatemp);
  void addValues(double etot, double volume, double pressure, double temp);
  void updateTemp(double& temperature);
};

#endif
