////////////////////////////////////////////////////////////////////////////////
//
// ConvergenceDetector.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConvergenceDetector.h,v 1.1 2007/07/24 16:40:14 draeger1 Exp $

#ifndef CONVERGENCEDETECTOR_H
#define CONVERGENCEDETECTOR_H

using namespace std;

class ConvergenceDetector {
  protected:

  double threshold_;

  public:

  virtual void addValue(double val) = 0;
  virtual bool isConverged(void) = 0;
  virtual ~ConvergenceDetector() {}
    
};
#endif
