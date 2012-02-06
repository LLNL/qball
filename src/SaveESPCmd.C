////////////////////////////////////////////////////////////////////////////////
//
// SaveESPCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveESPCmd.C,v 1.3 2010/01/07 18:01:48 draeger1 Exp $

#include "SaveESPCmd.h"
#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Context.h"
#include "AtomSet.h"
#include "Species.h"
#include<fstream>
#include<iostream>
using namespace std;

int SaveESPCmd::action(int argc, char **argv)
{

  if ( !(argc == 2) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: saveesp filename -->" << endl;
    return 1;
  }
  char* filename = argv[1];

  StructureFactor sf;
  ChargeDensity cd_(*s);
  cd_.update_density();
  EnergyFunctional ef_(*s,s->wf,cd_);
  const Context* ctxt_ = s->wf.spincontext(0);
  FourierTransform* ft_ = cd_.vft();
  Basis* vbasis_ = cd_.vbasis();
  const UnitCell& cell = s->wf.cell();
  const AtomSet& atoms = s->atoms;

  vector<vector<double> > tau0, rhops;
  vector<complex<double> > rhopst, rhog, v_g, v_r;
  const int nsp_  = atoms.nsp();
  const int ngloc = vbasis_->localsize();
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const double *const g2i = vbasis_->g2i_ptr();
  const double fpi = 4.0 * M_PI;

  ofstream os;
  os.open(filename,ofstream::out);    // text output
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);

  tau0.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ ) {
    const int na = atoms.na(is);
    tau0[is].resize(3*na);
  }
  atoms.get_positions(tau0,true);
  sf.init(tau0,*vbasis_);
  sf.update(tau0,*vbasis_);
  
  rhops.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ )
    rhops[is].resize(ngloc);

  for ( int is = 0; is < nsp_; is++ ) {
    Species *s = atoms.species_list[is];
    const double * const g = vbasis_->g_ptr();
    complex<double> *sfg = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ ) {
      const complex<double> sg = sfg[ig];
      rhops[is][ig] = s->rhopsg(g[ig]) * omega_inv;
    }    
  }
  rhopst.resize(ngloc);
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < atoms.nsp(); is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
    }
  }

  if (s->wf.nspin() == 2 && ctxt_->oncoutpe() )
    cout << "<WARNING> saveesp command only prints spin = 0 potential </WARNING>" << endl;

  // compute total charge density:  electronic + ionic cores
  rhog.resize(ngloc);
  if ( s->wf.nspin() == 1 )
    for ( int ig = 0; ig < ngloc; ig++ ) 
      rhog[ig] = omega_inv * cd_.rhog[0][ig] + rhopst[ig];
  else 
    for ( int ig = 0; ig < ngloc; ig++ ) 
      rhog[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] ) + rhopst[ig];

  // vhart_g = 4 * pi * (rhoel + rhops) * g2i
  v_g.resize(ft_->np012loc());
  for ( int ig = 0; ig < ngloc; ig++ ) {
    v_g[ig] = fpi * rhog[ig] * g2i[ig];
  }
 
  // FT to v_r
  v_r.resize(ft_->np012loc());
  ft_->backward(&v_g[0],&v_r[0]);

  //ewd DEBUG
  if ( ctxt_->oncoutpe() ) 
    cout << "Preparing to print out potential..." << endl;

  // send data in first column to pe 0 for output
  if (ctxt_->mycol() == 0) {
    vector<double> vrtmp(ft_->np012loc());
    for (int j = 0; j < ft_->np012loc(); j++)
      vrtmp[j] = real(v_r[j]);
    
    for ( int i = 0; i < ctxt_->nprow(); i++ ) {
      if ( i == ctxt_->myrow() ) {
        int size = ft_->np012loc();
        ctxt_->isend(1,1,&size,1,0,0);
        ctxt_->dsend(size,1,&vrtmp[0],1,0,0);
      }
    }
    if ( ctxt_->oncoutpe() ) {
      for ( int i = 0; i < ctxt_->nprow(); i++ ) {
        int size = 0;
        ctxt_->irecv(1,1,&size,1,i,0);
        ctxt_->drecv(size,1,&vrtmp[0],1,i,0);

        // print grid info in header (use Molmol format)
        if (i==0) {
          D3vector a0 = s->wf.cell().a(0);
          D3vector a1 = s->wf.cell().a(1);
          D3vector a2 = s->wf.cell().a(2);
          double dx = a0.x/(double)ft_->np0();
          double dy = a1.y/(double)ft_->np1();
          double dz = a2.z/(double)ft_->np2();
          //    origin    npoints       grid spacing
          os << "0.0 " << ft_->np0() << " " << dx << endl;
          os << "0.0 " << ft_->np1() << " " << dy << endl;
          os << "0.0 " << ft_->np2() << " " << dz << endl;
        }
        for (int j=0; j<size; j++) 
          os << vrtmp[j] << endl;
      }
    }
  }
  os.close();
  return 0;
}
