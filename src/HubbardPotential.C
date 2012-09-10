////////////////////////////////////////////////////////////////////////////////
//
// HubbardPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HubbardPotential.C,v 1.1 2010/08/26 17:44:42 draeger1 Exp $

#include "HubbardPotential.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Species.h"
#include "Basis.h"
#include "blas.h"
#include "D3vector.h"
#include <iomanip>
#include <vector>
using namespace std;

#if AIX || BGL
extern "C" void vsincos(double *x, double *y, double *z, int *n);
#endif


////////////////////////////////////////////////////////////////////////////////
HubbardPotential::~HubbardPotential(void) {
}

////////////////////////////////////////////////////////////////////////////////
void HubbardPotential::init(void) {

  nsp = atoms_.nsp();
  na.resize(nsp);
  hub_l_.resize(nsp);
  hub_u_.resize(nsp);
  hub_alpha_.resize(nsp);
  lmsize_.resize(nsp);
  int nsym = symset_.nsym()+1; // include identity, in last array index [nsym]
  
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_.species_list[is];
    na[is] = atoms_.na(is);

    if (s->hubbard_l() < 0)
      lmsize_[is] = 0;
    else {
      const int lm = 2*s->hubbard_l()+1;
      lmsize_[is] = lm;
    }
  }
  phiylm.resize(wf_.nspin());
  for (int ispin=0; ispin<wf_.nspin(); ispin++)
    phiylm[ispin].resize(wf_.nkptloc());
  for (int ispin=0; ispin<wf_.nspin(); ispin++)
    for (int ikp=0; ikp<wf_.nkptloc(); ikp++) 
      phiylm[ispin][ikp].resize(nsp);
  for (int ispin=0; ispin<wf_.nspin(); ispin++)
    if (wf_.spinactive(ispin)) {
      for (int ikp=0; ikp<wf_.nkptloc(); ikp++) {
        for ( int is = 0; is < nsp; is++ ) {
          const SlaterDet* sd_ = wf_.sdloc(ispin,ikp);
          const Basis& basis_ = sd_->basis();
          const int ngwl = basis_.localsize();
          phiylm[ispin][ikp][is].resize(lmsize_[is]*ngwl);
        }
      }
    }

  if (nsym > 0) {
    int lmax = 3;
    ylmsym.resize(lmax);
    for (int l=0; l<lmax; l++)
      ylmsym[l].resize(nsym);

    for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_.species_list[is];
      int hubl = s->hubbard_l();
      if (hubl >= 0) 
        for (int isym=0; isym<nsym; isym++)
          ylmsym[hubl][isym].resize(lmsize_[is]*lmsize_[is]);
    }

    rmrand.resize(2*lmax+1);
    rmrandlen.resize(2*lmax+1);
    rmsym.resize(2*lmax+1);
    rmsymlen.resize(2*lmax+1);
    for (int i=0; i<2*lmax+1; i++) {
      rmsym[i].resize(nsym);
      rmsymlen[i].resize(nsym);
    }
    
    int seedval = 185832792;
    srand48(seedval);
    for (int i=0; i<2*lmax+1; i++) {
      rmrand[i].x = drand48()-0.5;
      rmrand[i].y = drand48()-0.5;
      rmrand[i].z = drand48()-0.5;
      rmrandlen[i] = length(rmrand[i]);
    }

    // apply symmetry ops to rmrand coordinates
    for (int i=0; i<2*lmax+1; i++) 
      for (int isym=0; isym<nsym; isym++)
        rmsym[i][isym].resize(3);

    for (int i=0; i<2*lmax+1; i++) {
      for (int isym=0; isym<symset_.nsym(); isym++) {
        D3vector r_sym = symset_.symlist[isym]->applyToVector(rmrand[i],true);
        rmsym[i][isym][0] = r_sym.x;
        rmsym[i][isym][1] = r_sym.y;
        rmsym[i][isym][2] = r_sym.z;
        rmsymlen[i][isym] = length(r_sym);
      }
      // put identity in last index
      rmsym[i][nsym-1][0] = rmrand[i].x;
      rmsym[i][nsym-1][1] = rmrand[i].y;
      rmsym[i][nsym-1][2] = rmrand[i].z;
      rmsymlen[i][nsym-1] = rmrandlen[i];
    }    
  }
}

////////////////////////////////////////////////////////////////////////////////
void HubbardPotential::update_phiylm(void) {
  // update array phiylm[ispin][ikp][is][lm][ig] following a change of cell dimensions
  
  const double pi = M_PI;
  const double fpi = 4.0 * pi;
  const double s14pi = sqrt(1.0/fpi);
  const double s34pi = sqrt(3.0/fpi);  
  const double s54pi = sqrt(5.0/fpi);
  const double s20pi = sqrt(20.0*pi);
  const double s20pi3 = sqrt(20.0*pi/3.0);
  const double s3 = sqrt(3.0);
  const double s74pi = sqrt(7.0/fpi);
  const double s2132pi = sqrt(21.0/(32.*pi));
  const double s3532pi = sqrt(35.0/(32.*pi));
  const double s1054pi = sqrt(105.0/fpi);
  
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_.species_list[is];
    hub_l_[is] = s->hubbard_l();    
    hub_u_[is] = s->hubbard_u();    
    hub_alpha_[is] = s->hubbard_alpha();    
  }

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for (int kloc=0; kloc<wf_.nkptloc(); kloc++) {
        const SlaterDet* sd_ = wf_.sdloc(ispin,kloc);
        const Basis& basis_ = sd_->basis();
        const int ngwl = basis_.localsize();

        const double omega = basis_.cell().volume();
        assert(omega != 0.0);
        const double omega_inv = 1.0 / omega;
        const double fpisqrtvol = fpi*sqrt(omega_inv);  // normalization factor for phi(G)
        
        const double *kpg   = basis_.kpg_ptr();
        const double *kpg2  = basis_.kpg2_ptr();
        const double *kpgi  = basis_.kpgi_ptr();
        const double *kpg_x = basis_.kpgx_ptr(0);
        const double *kpg_y = basis_.kpgx_ptr(1);
        const double *kpg_z = basis_.kpgx_ptr(2);

        // compute phiylm
        for ( int is = 0; is < nsp; is++ ) {
           Species *s = atoms_.species_list[is];

          if ( hub_l_[is] == 0 ) {
            // fill phiylm with phi_g*Y_lm at basis set vectors

            double *p0 = &phiylm[ispin][kloc][is][0];
            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double tg = kpg[ig];
              double phig_;
              s->phig(tg,phig_);
              phig_ *= fpisqrtvol;
              p0[ig] = s14pi * phig_;
            }
          }
          else if ( hub_l_[is] == 1 ) {
            double *p1 = &phiylm[ispin][kloc][is][0*ngwl];
            double *p2 = &phiylm[ispin][kloc][is][1*ngwl];
            double *p3 = &phiylm[ispin][kloc][is][2*ngwl];

            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double tg = kpg[ig];
              double phig_;
              s->phig(tg,phig_);
              phig_ *= fpisqrtvol;
              
              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgi = kpgi[ig];
              
              const double y1 = s34pi * tgx * tgi;
              const double y2 = s34pi * tgy * tgi;
              const double y3 = s34pi * tgz * tgi;
              
              p1[ig]  = y1 * phig_;
              p2[ig]  = y2 * phig_;
              p3[ig]  = y3 * phig_;
            }
          }
          else if ( hub_l_[is] == 2 ) {
            double *p4 = &phiylm[ispin][kloc][is][0*ngwl];
            double *p5 = &phiylm[ispin][kloc][is][1*ngwl];
            double *p6 = &phiylm[ispin][kloc][is][2*ngwl];
            double *p7 = &phiylm[ispin][kloc][is][3*ngwl];
            double *p8 = &phiylm[ispin][kloc][is][4*ngwl];

            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double tg = kpg[ig];
              double phig_;
              s->phig(tg,phig_);
              phig_ *= fpisqrtvol;
              
              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tgi = kpgi[ig];
              const double tg2 = tg * tg;
              const double tgi2 = tgi * tgi;
              
              const double tgxx = tgx2 * tgi2;
              const double tgyy = tgy2 * tgi2;
              const double tgzz = tgz2 * tgi2;
              const double tgxy = tgx * tgy * tgi2;
              const double tgyz = tgy * tgz * tgi2;
              const double tgxz = tgx * tgz * tgi2;
        
              const double y4 = s54pi * 0.5 * (3.0 * tgzz - 1.0 );
              const double y5 = s54pi * 0.5 * s3 * ( tgxx - tgyy );
              const double y6 = s54pi * s3 * tgxy;
              const double y7 = s54pi * s3 * tgyz;
              const double y8 = s54pi * s3 * tgxz;
              
              p4[ig]  = y4 * phig_;
              p5[ig]  = y5 * phig_;
              p6[ig]  = y6 * phig_;
              p7[ig]  = y7 * phig_;
              p8[ig]  = y8 * phig_;
            }
          }
          else if ( hub_l_[is] == 3 ) {
            double *p9 = &phiylm[ispin][kloc][is][0*ngwl];
            double *p10 = &phiylm[ispin][kloc][is][1*ngwl];
            double *p11 = &phiylm[ispin][kloc][is][2*ngwl];
            double *p12 = &phiylm[ispin][kloc][is][3*ngwl];
            double *p13 = &phiylm[ispin][kloc][is][4*ngwl];
            double *p14 = &phiylm[ispin][kloc][is][5*ngwl];
            double *p15 = &phiylm[ispin][kloc][is][6*ngwl];
            
            for ( int ig = 0; ig < ngwl; ig++ ) {
              const double tg = kpg[ig];
              double phig_;
              s->phig(tg,phig_);
              phig_ *= fpisqrtvol;

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tgi = kpgi[ig];
              const double tg2 = tg * tg;
              const double tgi2 = tgi * tgi;
              const double tgi3 = tgi2 * tgi;
              
              const double tgxx = tgx2 * tgi2;
              const double tgyy = tgy2 * tgi2;
              const double tgzz = tgz2 * tgi2;
              
              // real spherical harmonics:
              //   sqrt[ ((2l+1)/4pi) * (l-|m|)!/(l+|m|)! ] * P_l^|m|(cos(theta)) * F(phi)
              //     where F(phi) = sqrt(2)*cos(m*phi)   m > 0
              //                  = 1                    m = 0
              //                  = sqrt(2)*sin(|m|*phi) m < 0
              //
              const double y9 = s74pi * 0.5 * tgz * tgi * (5.0 * tgzz - 3.0 );
              const double y10 = s2132pi * tgx * tgi * ( 5.0 * tgzz - 1.0 );
              const double y11 = s2132pi * tgy * tgi * ( 5.0 * tgzz - 1.0 );
              const double y12 = s1054pi * tgx * tgy * tgz * tgi3;
              const double y13 = s1054pi * 0.5 * tgz * tgi * (tgxx - tgyy);
              const double y14 = s3532pi * tgx * tgi * ( tgxx - 3.0*tgyy);
              const double y15 = s3532pi * tgy * tgi * ( 3.0*tgxx - tgyy);
              
              p9[ig]  = y9 * phig_;
              p10[ig]  = y10 * phig_;
              p11[ig]  = y11 * phig_;
              p12[ig]  = y12 * phig_;
              p13[ig]  = y13 * phig_;
              p14[ig]  = y14 * phig_;
              p15[ig]  = y15 * phig_;
            }
          }
          else if ( hub_l_[is] > 3 ) {
            assert(false);
          }
        }
      }
    }
  }


  // calculate Ylm transform matrices to apply symmetry operations to hub_occ

  if (symset_.nsym() > 0) {
    int nsym = symset_.nsym()+1; // include identity, in last array index [nsym]
    for ( int is = 0; is < nsp; is++ ) {
      int hubl = hub_l_[is];
      int lmsize = lmsize_[is];

      if (hubl >= 0) {
        vector<double> ylm(lmsize*lmsize);
        vector<double> ylminv(lmsize*lmsize);
        vector<vector<double> > ylmstar;
        ylmstar.resize(nsym);
        for (int i=0; i<nsym; i++)
          ylmstar[i].resize(lmsize*lmsize);
      
        if (hubl == 0) {
          for (int m1=0; m1<lmsize; m1++) 
            ylm[m1*lmsize + 0] = s14pi;
          for (int isym=0; isym<nsym; isym++) 
            for (int m1=0; m1<lmsize; m1++) 
              ylmstar[isym][m1*lmsize+0] = s14pi;
        }
        else if (hubl == 1) {
          for (int m1=0; m1<lmsize; m1++) {
            ylm[m1*lmsize + 0] = s34pi*rmrand[m1].x/rmrandlen[m1];
            ylm[m1*lmsize + 1] = s34pi*rmrand[m1].y/rmrandlen[m1];
            ylm[m1*lmsize + 2] = s34pi*rmrand[m1].z/rmrandlen[m1];
          }
          for (int isym=0; isym<nsym; isym++) {
            for (int m1=0; m1<lmsize; m1++) {
              ylmstar[isym][m1*lmsize + 0] = s34pi*rmsym[m1][isym][0]/rmsymlen[m1][isym];
              ylmstar[isym][m1*lmsize + 1] = s34pi*rmsym[m1][isym][1]/rmsymlen[m1][isym];
              ylmstar[isym][m1*lmsize + 2] = s34pi*rmsym[m1][isym][2]/rmsymlen[m1][isym];
            }
          }
        }
        else if (hubl == 2) {
          for (int m1=0; m1<lmsize; m1++) {
            double x = rmrand[m1].x;
            double y = rmrand[m1].y;
            double z = rmrand[m1].z;
            double r = rmrandlen[m1];
            double r2inv = 1./(r*r);
            ylm[m1*lmsize + 0] = s54pi * 0.5 * (3.0 * z*z*r2inv - 1.0 );
            ylm[m1*lmsize + 1] = s54pi * 0.5 * s3 * ( x*x*r2inv - y*y*r2inv );
            ylm[m1*lmsize + 2] = s54pi * s3 * x*y*r2inv;
            ylm[m1*lmsize + 3] = s54pi * s3 * y*z*r2inv;
            ylm[m1*lmsize + 4] = s54pi * s3 * x*z*r2inv;
            for (int isym=0; isym<nsym; isym++) {
              for (int m1=0; m1<lmsize; m1++) {
                x = rmsym[m1][isym][0];
                y = rmsym[m1][isym][1];
                z = rmsym[m1][isym][2];
                r = rmsymlen[m1][isym];
                r2inv = 1./(r*r);
                ylmstar[isym][m1*lmsize + 0] = s54pi * 0.5 * (3.0 * z*z*r2inv - 1.0 );
                ylmstar[isym][m1*lmsize + 1] = s54pi * 0.5 * s3 * ( x*x*r2inv - y*y*r2inv );
                ylmstar[isym][m1*lmsize + 2] = s54pi * s3 * x*y*r2inv;
                ylmstar[isym][m1*lmsize + 3] = s54pi * s3 * y*z*r2inv;
                ylmstar[isym][m1*lmsize + 4] = s54pi * s3 * x*z*r2inv;
              }
            }
          }
        }
        else if (hubl == 3) {
          for (int m1=0; m1<lmsize; m1++) {
            double x = rmrand[m1].x;
            double y = rmrand[m1].y;
            double z = rmrand[m1].z;
            double rinv = 1./rmrandlen[m1];
            double r2inv = rinv*rinv;
            ylm[m1*lmsize + 0] = s74pi * 0.5 * z * rinv * (5.0 * z*z*r2inv - 3.0 );
            ylm[m1*lmsize + 1] = s2132pi * x * rinv * ( 5.0 * z*z*r2inv - 1.0 );
            ylm[m1*lmsize + 2] = s2132pi * y * rinv * ( 5.0 * z*z*r2inv - 1.0 );
            ylm[m1*lmsize + 3] = s1054pi * x * y * z * rinv*r2inv;
            ylm[m1*lmsize + 4] = s1054pi * 0.5 * z * rinv * (x*x*r2inv - y*y*r2inv);
            ylm[m1*lmsize + 5] = s3532pi * x * rinv * ( x*x*r2inv - 3.0*y*y*r2inv);
            ylm[m1*lmsize + 6] = s3532pi * y * rinv * ( 3.0*x*x*r2inv - y*y*r2inv);
            for (int isym=0; isym<nsym; isym++) {
              for (int m1=0; m1<lmsize; m1++) {
                x = rmsym[m1][isym][0];
                y = rmsym[m1][isym][1];
                z = rmsym[m1][isym][2];
                rinv = 1./rmsymlen[m1][isym];
                r2inv = rinv*rinv;
                ylmstar[isym][m1*lmsize + 0] = s74pi * 0.5 * z * rinv * (5.0 * z*z*r2inv - 3.0 );
                ylmstar[isym][m1*lmsize + 1] = s2132pi * x * rinv * ( 5.0 * z*z*r2inv - 1.0 );
                ylmstar[isym][m1*lmsize + 2] = s2132pi * y * rinv * ( 5.0 * z*z*r2inv - 1.0 );
                ylmstar[isym][m1*lmsize + 3] = s1054pi * x * y * z * rinv*r2inv;
                ylmstar[isym][m1*lmsize + 4] = s1054pi * 0.5 * z * rinv * (x*x*r2inv - y*y*r2inv);
                ylmstar[isym][m1*lmsize + 5] = s3532pi * x * rinv * ( x*x*r2inv - 3.0*y*y*r2inv);
                ylmstar[isym][m1*lmsize + 6] = s3532pi * y * rinv * ( 3.0*x*x*r2inv - y*y*r2inv);
              }
            }
          }
        }

        for (int i=0; i<lmsize*lmsize; i++)
          ylminv[i] = ylm[i];

        int info;
        vector<int> ipiv(lmsize);
        double* pylm = &ylminv[0];
        dgetrf(&lmsize,&lmsize,pylm,&lmsize,&ipiv[0],&info);

        vector<double> work(1);
        int lwork = -1;
        // First call to compute optimal size of work array, returned in work[0]
        dgetri(&lmsize, &ylminv[0], &lmsize, &ipiv[0], &work[0], &lwork, &info);
        lwork = (int) work[0] + 1;
        work.resize(lwork);
        dgetri(&lmsize, &ylminv[0], &lmsize, &ipiv[0], &work[0], &lwork, &info);

        for (int isym=0; isym<nsym; isym++) {
          vector<double> ylmmat(lmsize*lmsize);
          double one = 1.0, zero = 0.0;
          char cn='n';
          dgemm(&cn,&cn,&lmsize,&lmsize,&lmsize,&one,&ylmstar[isym][0],&lmsize,&ylminv[0],&lmsize,
                &zero,&ylmmat[0],&lmsize);

          for (int i=0; i<lmsize*lmsize; i++)
            ylmsym[hubl][isym][i] = ylmmat[i];
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double HubbardPotential::energy(bool compute_hpsi, Wavefunction& dwf)
{
  // define atom block size
  const int na_block_size = 32;

  vector<vector<double> > tau;
  atoms_.get_positions(tau,true);

  vector<vector<vector<vector<vector<double> > > > > hub_occ;
  hub_occ.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    hub_occ[ispin].resize(nsp);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int is = 0; is < nsp; is++ )
      hub_occ[ispin][is].resize(na[is]);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int is = 0; is < nsp; is++ )
      for (int ia=0; ia<na[is]; ia++)
        hub_occ[ispin][is][ia].resize(lmsize_[is]);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int is = 0; is < nsp; is++ )
      for (int ia=0; ia<na[is]; ia++)
        for (int m1=0; m1<lmsize_[is]; m1++) 
          hub_occ[ispin][is][ia][m1].resize(lmsize_[is]);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int is = 0; is < nsp; is++ )
      for (int ia=0; ia<na[is]; ia++)
        for (int m1=0; m1<lmsize_[is]; m1++)
          for (int m2=0; m2<lmsize_[is]; m2++) 
            hub_occ[ispin][is][ia][m1][m2] = 0.0;

  double ehub[2] = {0.0,0.0};
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {

      vector<vector<double> > hubproj_loc_gamma;
      vector<vector<complex<double> > > hubproj_loc;
      vector<double> hubproj_buf_gamma;
      vector<complex<double> > hubproj_buf;
      vector<vector<double> > wfhubloc_gamma;
      vector<vector<complex<double> > > wfhubloc;
      hubproj_loc_gamma.resize(wf_.nkptloc());
      hubproj_loc.resize(wf_.nkptloc());
      wfhubloc_gamma.resize(wf_.nkptloc());
      wfhubloc.resize(wf_.nkptloc());

      double wtsum = 0.0;
      for (int kloc=0; kloc<wf_.nkptloc(); kloc++) {
        double wt = wf_.weight(wf_.kptloc(kloc))/wf_.weightsum();
        wtsum += wt;
        
        assert(wf_.sdloc(ispin,kloc) != 0);
        const SlaterDet& sd_ = *(wf_.sdloc(ispin,kloc));
        const ComplexMatrix& c = sd_.c();
        const Basis& basis_ = wf_.sdloc(ispin,kloc)->basis();
        const int ngwl = basis_.localsize();
        const int nstloc = sd_.nstloc();
        const vector<double>& occ = sd_.occ();

        vector<double> kpgr(na_block_size*ngwl); // kpgr[ig+ia*ngwl]
        vector<double> ckpgr(na_block_size*ngwl); // ckpgr[ig+ia*ngwl]
        vector<double> skpgr(na_block_size*ngwl); // skpgr[ig+ia*ngwl]
  
        for ( int is = 0; is < nsp; is++ ) {  
          if ( lmsize_[is] > 0 ) { // species is has Hubbard l defined

            // define number of atom blocks                                        
            const int na_blocks = na[is] / na_block_size +                         
                            ( na[is] % na_block_size == 0 ? 0 : 1 );         

            // calculate localized wavefunction projection operators for each ion center I
            // wfhubloc(g,I,m) = (i)^l Y_lm(g) phi_l(g) exp(-i(k+G)*R_I)
            if (basis_.real())
              wfhubloc_gamma[kloc].resize(lmsize_[is]*na_block_size*2*ngwl);
            else 
              wfhubloc[kloc].resize(lmsize_[is]*na_block_size*ngwl);        

            // apply projection to Kohn-Sham orbitals for each (2l+1) m of Hubbard l
            // hubproj(istate,I,m) = wfhubloc(istate,I,m) * sd_.c(istate)
            if (basis_.real()) {
              hubproj_loc_gamma[kloc].resize(lmsize_[is]*na_block_size*nstloc);
              hubproj_buf_gamma.resize(lmsize_[is]*na_block_size*nstloc);
            }
            else {
              hubproj_loc[kloc].resize(lmsize_[is]*na_block_size*nstloc);
              hubproj_buf.resize(lmsize_[is]*na_block_size*nstloc);
            }

            // compute DFT+U occupation numbers for localized Hubbard wf orbitals
            // n(m1,m2,I) = SUM_states occ(istate) conj(proj(istate,I,m1)) * proj(istate,I,m2)
            for ( int ia_block = 0; ia_block < na_blocks; ia_block++ ) {
              // process projectors of atoms in block ia_block             

              const int iastart = ia_block * na_block_size;                
              const int iaend = (ia_block+1) * na_block_size < na[is] ?    
                  (ia_block+1) * na_block_size : na[is];                                      
              const int ia_block_size = iaend - iastart;                   
                                                                     
              // compute ckpgr[is][ia][ig], skpgr[is][ia][ig]                  
              int k = 3;                                                   
              double mone = -1.0, zero = 0.0;
              char cn='n';

              // next line: const cast is ok since dgemm_ does not modify argument 
              double* kpgx = const_cast<double*>(basis_.kpgx_ptr(0));
              dgemm(&cn,&cn,(int*)&ngwl,(int*)&ia_block_size,&k,&mone,             
                    kpgx,(int*)&ngwl, &tau[is][3*iastart],&k,                        
                    &zero,&kpgr[0],(int*)&ngwl);                                     
              
              int len = ia_block_size * ngwl;                                      
#if AIX || BGL  
              vsincos(&skpgr[0],&ckpgr[0],&kpgr[0],&len);                                
#else
              for ( int i = 0; i < len; i++ ) {
                const double arg = kpgr[i];
                skpgr[i] = sin(arg);                                                 
                ckpgr[i] = cos(arg);                                                 
              }                                                                    
#endif
              
              for ( int lm = 0; lm < lmsize_[is]; lm++ ) {

                // phiylm[ispin][kloc][is][ig+ngwl*lm]
                const double * p = &phiylm[ispin][kloc][is][ngwl*lm];                            

                for ( int ia = 0; ia < ia_block_size; ia++ ) {
                  double* wfh;
                  if (basis_.real())
                    wfh = &wfhubloc_gamma[kloc][2*(ia+lm*ia_block_size)*ngwl];
                  else
                    wfh = (double*)&wfhubloc[kloc][(ia+lm*ia_block_size)*ngwl];
                  const double* c = &ckpgr[ia*ngwl];
                  const double* s = &skpgr[ia*ngwl];
                  if ( hub_l_[is] == 0 ) {
                    for ( int ig = 0; ig < ngwl; ig++ ) {
                      wfh[2*ig]   = p[ig] * c[ig];                                   
                      wfh[2*ig+1] = p[ig] * s[ig];                                   
                    }
                  }
                  else if ( hub_l_[is] == 1 ) {
                    /* Next line: -i * eigr */                                   
                    /* -i * (a+i*b) = b - i*a */                                 
                    for ( int ig = 0; ig < ngwl; ig++ ) {
                      wfh[2*ig]   =  p[ig] * s[ig];
                      wfh[2*ig+1] = -p[ig] * c[ig];
                    }
                  }
                  else if ( hub_l_[is] == 2 ) {
                    // Next line: (-) sign for -eigr                             
                    /* -1 * (a+i*b) = -a - i*b */                                 
                    for ( int ig = 0; ig < ngwl; ig++ ) {
                      wfh[2*ig]   = -p[ig] * c[ig];
                      wfh[2*ig+1] = -p[ig] * s[ig];
                    }
                  }
                  else if ( hub_l_[is] == 3 ) {
                    /* i * (a+i*b) = -b + i*a */ 
                    for ( int ig = 0; ig < ngwl; ig++ ) {
                      wfh[2*ig]   = -p[ig] * s[ig];
                      wfh[2*ig+1] = p[ig] * c[ig];
                    }
                  }
                }
              }
                
              // compute proj[nlma][nstloc] = wfhub^T * c                             
              double one=1.0;
              char ct='t';
              int twongwl = 2 * ngwl;
              int nlmnaloc = ia_block_size * lmsize_[is]; 
              int c_lda;
              const complex<double>* cp = c.cvalptr();
        
              if (basis_.real()) {
                c_lda = 2*c.mloc();                                        
                dgemm(&ct,&cn,&nlmnaloc,(int*)&nstloc,&twongwl,&one,
                      &wfhubloc_gamma[kloc][0],&twongwl, (double*)cp, &c_lda,
                      &zero,&hubproj_loc_gamma[kloc][0],&nlmnaloc);
              }
              else {
                c_lda = c.mloc();                                        
                complex<double> zzero = complex<double>(0.0,0.0);
                complex<double> zone = complex<double>(1.0,0.0);
                char cc='c';
                zgemm(&cc,&cn,&nlmnaloc,(int*)&nstloc,(int*)&ngwl,&zone,
                      &wfhubloc[kloc][0],(int *)&ngwl, (complex<double> *)cp, &c_lda,
                      &zzero,&hubproj_loc[kloc][0],&nlmnaloc);
              }
              
              // for k=(0,0,0):  need to correct for double counting of G=0, i.e. if ctxt_.myrow() == 0
              if (basis_.real()) { 
                if ( ctxt_.myrow() == 0 ) {
                  // rank-one update                                                 
                  // dger(m,n,alpha,x,incx,y,incy,a,lda);                            
                  // a += alpha * x * transpose(y)                                   
                  // x = first row of wfhubloc                                        
                  // y^T = first row of c                                            
                  double alpha = -0.5;                                               
                  dger(&nlmnaloc,(int*)&nstloc,&alpha,&wfhubloc_gamma[kloc][0],&twongwl,
                       (double*)cp,&c_lda,&hubproj_loc_gamma[kloc][0],&nlmnaloc);
                }
              }
              
              if (nstloc > 0) {
                // Allreduce hubproj partial sum
                MPI_Comm basis_comm = basis_.context().comm();                       
                int hubproj_size = nlmnaloc*nstloc;                                   
                if (basis_.real()) {
                  MPI_Allreduce(&hubproj_loc_gamma[kloc][0],&hubproj_buf_gamma[0],hubproj_size,
                                MPI_DOUBLE,MPI_SUM,basis_comm);                        
                }
                else {
                  MPI_Allreduce(&hubproj_loc[kloc][0],&hubproj_buf[0],hubproj_size,
                                MPI_DOUBLE_COMPLEX,MPI_SUM,basis_comm);                        
                }
              }
                            
              // factor 2.0 in next line is: counting G, -G                        
              if (basis_.real())
                for (int i=0; i<hubproj_loc_gamma[kloc].size(); i++) 
                  hubproj_loc_gamma[kloc][i] = 2.0 * hubproj_buf_gamma[i];
              else
                for (int i=0; i<hubproj_loc[kloc].size(); i++) 
                  hubproj_loc[kloc][i] = hubproj_buf[i];

              // n(m1,m2,I) = SUM_states occ(istate) conj(proj(istate,I,m1)) * proj(istate,I,m2)
              if (nstloc > 0) {
                for (int ia=0; ia<ia_block_size; ia++) {
                  for (int m1=0; m1<lmsize_[is]; m1++) {
                    for (int m2=0; m2<lmsize_[is]; m2++) {
                      const int nbase = sd_.context().mycol() * c.nb();                      
                      if (basis_.real()) {
                        for ( int n = 0; n < nstloc; n++ ) {
                          const double sdocc = occ[n + nbase];
                          const int i1 = ia + m1*ia_block_size + n * nlmnaloc;           
                          const int i2 = ia + m2*ia_block_size + n * nlmnaloc;           
                          hub_occ[ispin][is][ia][m1][m2] += wt * sdocc * hubproj_loc_gamma[kloc][i1] * hubproj_loc_gamma[kloc][i2];
                        }
                      }
                      else {
                        for ( int n = 0; n < nstloc; n++ ) {
                          const double sdocc = occ[n + nbase];                        
                          const int i1 = ia + m1*ia_block_size + n * nlmnaloc;
                          const int i2 = ia + m2*ia_block_size + n * nlmnaloc;
                          hub_occ[ispin][is][ia][m1][m2] += wt * sdocc * real(conj(hubproj_loc[kloc][i1])*hubproj_loc[kloc][i2]);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }


      // sum hub_occ over k-points and states
      wf_.spincontext(ispin)->dsum('r',1,1,&wtsum,1);    
      for ( int is = 0; is < nsp; is++ ) {  
        if ( lmsize_[is] > 0 ) { // species is has Hubbard l defined
          int occsize = lmsize_[is]*lmsize_[is]*na[is];
          vector<double> occtmp(occsize);
          for (int ia=0; ia<na[is]; ia++) 
            for (int m1=0; m1<lmsize_[is]; m1++) 
              for (int m2=0; m2<lmsize_[is]; m2++) 
                occtmp[ia*lmsize_[is]*lmsize_[is] + lmsize_[is]*m1 + m2] = hub_occ[ispin][is][ia][m1][m2];                    
          double *occp = (double*) &occtmp[0];              
          wf_.spincontext(ispin)->dsum('r',occsize,1,&occp[0],occsize);    

          for (int ia=0; ia<na[is]; ia++) 
            for (int m1=0; m1<lmsize_[is]; m1++) 
              for (int m2=0; m2<lmsize_[is]; m2++) 
                hub_occ[ispin][is][ia][m1][m2] = occtmp[ia*lmsize_[is]*lmsize_[is] + lmsize_[is]*m1 + m2];
        }
      }

      // symmetrize hub_occ
      if (symset_.nsym() > 0) {
        atoms_.findSymmetricAtoms(symset_);
          
        vector<vector<vector<double> > > hub_occ_sym;
        hub_occ_sym.resize(nsp);
        for ( int is = 0; is < nsp; is++ )
          hub_occ_sym[is].resize(na[is]);
        for ( int is = 0; is < nsp; is++ )
          for (int ia=0; ia<na[is]; ia++)
            hub_occ_sym[is][ia].resize(lmsize_[is]*lmsize_[is]);
        for ( int is = 0; is < nsp; is++ )
          for (int ia=0; ia<na[is]; ia++)
            for (int m=0; m<lmsize_[is]*lmsize_[is]; m++) 
              hub_occ_sym[is][ia][m] = 0.0;

        for ( int is = 0; is < nsp; is++ ) {  
          if ( lmsize_[is] > 0 ) { // species is has Hubbard l defined
            const int hubl = hub_l_[is];
            const int lmsize = lmsize_[is];
            for (int ia=0; ia<na[is]; ia++) {
              for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) { 
                for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) { 
                  double symsum = 0.0;
                  for (int isym=0; isym<symset_.nsym(); isym++) {
                    int ja = atoms_.symatomid(is,ia,isym);
                    for ( int m3 = 0; m3 < lmsize_[is]; m3++ ) 
                      for ( int m4 = 0; m4 < lmsize_[is]; m4++ )
                        symsum += ylmsym[hubl][isym][lmsize*m1+m3]*hub_occ[ispin][is][ja][m3][m4]*ylmsym[hubl][isym][lmsize*m2+m4];
                  }
                  // add identity term
                  int nsym = symset_.nsym() + 1;
                  for ( int m3 = 0; m3 < lmsize_[is]; m3++ ) 
                    for ( int m4 = 0; m4 < lmsize_[is]; m4++ )
                      symsum += ylmsym[hubl][nsym-1][lmsize*m1+m3]*hub_occ[ispin][is][ia][m3][m4]*ylmsym[hubl][nsym-1][lmsize*m2+m4];
                  hub_occ_sym[is][ia][m1*lmsize+m2] = symsum/(double)nsym;
                }
              }
            }
            for (int ia=0; ia<na[is]; ia++)
              for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) 
                for ( int m2 = 0; m2 < lmsize_[is]; m2++ )
                  hub_occ[ispin][is][ia][m1][m2] = hub_occ_sym[is][ia][m1*lmsize_[is]+m2];
          }
        }
      }
            
      // print out hub_occ
      if (ctxt_.mype() == 0) { 
        //if (wf_.spincontext(ispin)->myproc()==0) { // output gets interleaved, need to send ispin 1 to pe 0
        for ( int is = 0; is < nsp; is++ ) {  
          if ( lmsize_[is] > 0 ) {     // species is has Hubbard l defined
            for (int ia=0; ia<na[is]; ia++) {
              cout << " Hubbard occ for species " << is << ", atom " << ia << ", spin " << ispin << endl;
              for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) { 
                cout << "    ";
                for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) { 
                  cout << hub_occ[ispin][is][ia][m1][m2] << "  ";
                }
                cout << endl;
              }
            }
          }
        }
      }
        


      
      // accumulate Ehub contribution
      for (int kloc=0; kloc<wf_.nkptloc(); kloc++) {
        double wt = wf_.weight(wf_.kptloc(kloc))/(wf_.weightsum()*wtsum);
        
        assert(wf_.sdloc(ispin,kloc) != 0);
        const SlaterDet& sd_ = *(wf_.sdloc(ispin,kloc));
        const ComplexMatrix& c = sd_.c();
        const Basis& basis_ = wf_.sdloc(ispin,kloc)->basis();
        const int ngwl = basis_.localsize();
        int twongwl = 2 * ngwl;
        const int nstloc = sd_.nstloc();

        for ( int is = 0; is < nsp; is++ ) {  
          if ( lmsize_[is] > 0 ) { // species is has Hubbard l defined
            // define number of atom blocks                                        
            const int na_blocks = na[is] / na_block_size +                         
                            ( na[is] % na_block_size == 0 ? 0 : 1 );         
            for ( int ia_block = 0; ia_block < na_blocks; ia_block++ ) {
              // process projectors of atoms in block ia_block             
              const int iastart = ia_block * na_block_size;                
              const int iaend = (ia_block+1) * na_block_size < na[is] ?    
                  (ia_block+1) * na_block_size : na[is];                                      
              const int ia_block_size = iaend - iastart;                   
        
              int nlmnaloc = ia_block_size * lmsize_[is]; 

              for (int ia=0; ia<ia_block_size; ia++) {
                for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) { 
                  ehub[ispin] += wt * ( hub_alpha_[is] + 0.5*hub_u_[is] ) * hub_occ[ispin][is][ia][m1][m1];
                  for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) { 
                    ehub[ispin] -= wt * 0.5*hub_u_[is] * hub_occ[ispin][is][ia][m1][m2] * hub_occ[ispin][is][ia][m2][m1];
                  }
                }
              }

              if ( compute_hpsi ) {
                SlaterDet& dsd_ = *(dwf.sdloc(ispin,kloc));
                ComplexMatrix& dc = dsd_.c();
                complex<double>* dcp = dc.valptr();
                vector<double> hpsiocc(lmsize_[is]);
                if (basis_.real()) {
                  for ( int n = 0; n < nstloc; n++ ) {
                    for (int ia=0; ia<ia_block_size; ia++) {
                      for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) { 

                        for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) 
                          hpsiocc[m2] = -1. * hub_u_[is] * hub_occ[ispin][is][ia][m1][m2];
                        hpsiocc[m1] += ( hub_alpha_[is] + 0.5*hub_u_[is] ) * hub_occ[ispin][is][ia][m1][m1];

                        double mult = 0.0;
                        for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) {
                          const int i2 = ia + m2*ia_block_size + n * nlmnaloc;
                          mult += hpsiocc[m2]*hubproj_loc_gamma[kloc][i2];
                        }
                        const int locind = 2*(ia+m1*ia_block_size)*ngwl;
                        const int cpind = dc.mloc()*n;
                        int inc1 = 1;
                        daxpy(&twongwl,&mult,(double*)&wfhubloc_gamma[kloc][locind],&inc1,(double*)&dcp[cpind],&inc1);
                      }
                    }
                  }
                }
                else {
                  for ( int n = 0; n < nstloc; n++ ) {
                    for (int ia=0; ia<ia_block_size; ia++) {
                      for ( int m1 = 0; m1 < lmsize_[is]; m1++ ) { 

                        for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) 
                          hpsiocc[m2] = -1. * hub_u_[is] * hub_occ[ispin][is][ia][m1][m2];
                        hpsiocc[m1] += ( hub_alpha_[is] + 0.5*hub_u_[is] ) * hub_occ[ispin][is][ia][m1][m1];

                        complex<double> mult = complex<double>(0.0,0.0);
                        for ( int m2 = 0; m2 < lmsize_[is]; m2++ ) {
                          const int i2 = ia + m2*ia_block_size + n * nlmnaloc;
                          mult += hpsiocc[m2]*hubproj_loc[kloc][i2];
                        }
                        const int locind = (ia+m1*ia_block_size)*ngwl;
                        const int cpind = dc.mloc()*n;
                        int inc1 = 1;
                        zaxpy((int*)&ngwl,&mult,&wfhubloc[kloc][locind],&inc1,&dcp[cpind],&inc1);
                      }
                    }
                  }
                }
              }
            } // for ia_block
          } // lmsize_[is]>0
        } // is
      }
    }
  }
  // reduction of ehub across rows (sum over distributed states, k-points and spin)
  wf_.wfcontext()->dsum('r',2,1,&ehub[0],2);
  
  double ehubtot = ehub[0] + ehub[1];
  
  return ehubtot;
}
