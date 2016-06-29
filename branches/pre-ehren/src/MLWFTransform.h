////////////////////////////////////////////////////////////////////////////////
//
// MLWFTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MLWFTransform.h,v 1.1 2008/05/13 18:28:54 draeger1 Exp $

#ifndef MLWFTRANSFORM_H
#define MLWFTRANSFORM_H

#include <vector>
#include <complex>
class SlaterDet;
class UnitCell;
class DoubleMatrix;
#include "D3vector.h"
#include "BasisMapping.h"

class MLWFTransform
{
  private:

  const SlaterDet& sd_;
  const UnitCell& cell_;
  const Context& ctxt_;

  BasisMapping bm_;
  std::vector<DoubleMatrix*> a_;  // cosine and sine matrices
  DoubleMatrix* u_;               // orthogonal transformation
  std::vector<std::vector<double> > adiag_; // diagonal elements

  public:

  void compute_transform(void);
  void compute_sincos(const int n, const std::complex<double>* f,
    std::complex<double>* fc, std::complex<double>* fs);
  void apply_transform(SlaterDet& sd);

  double spread2(int i, int j);
  double spread2(int i);
  double spread2(void);
  double spread(int i);
  double spread(void);
  D3vector center(int i);
  D3vector dipole(void);

  MLWFTransform(const SlaterDet& sd);
  ~MLWFTransform(void);
};
#endif
