////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionHandler.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#ifndef WavefunctionHANDLER_H
#define WavefunctionHANDLER_H

#include "StructureHandler.h"
#include "UnitCell.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class FourierTransform;

class WavefunctionHandler : public StructureHandler
{
  private:

  Wavefunction& wf_;
  DoubleMatrix& gfdata_;
  UnitCell uc;
  UnitCell ruc;
  double ecut;
  // dmat[ispin][ikp][i]
  std::vector<std::vector<std::vector<double> > > &dmat_;
  int& nx_;
  int& ny_;
  int& nz_;
  int current_gf_nx,current_gf_ny,current_gf_nz;
  std::string current_gf_encoding;
  int current_ispin,current_ikp,current_n,current_igf;
  std::vector<double> dmat_tmp;
  double current_kx, current_ky, current_kz, current_weight;
  int current_size;
  int read_from_gfdata;
  FourierTransform* ft;
  std::vector<std::complex<double> > wftmp;

  void byteswap_double(size_t n, double* x);

  public:

  // Start of the root element in the structure being handled
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of the root element in the structure being handled
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  WavefunctionHandler(Wavefunction& wf, DoubleMatrix& gfdata,
    int& nx, int& ny, int& nz,
    std::vector<std::vector<std::vector<double> > > &dmat);
  ~WavefunctionHandler();
};
#endif
