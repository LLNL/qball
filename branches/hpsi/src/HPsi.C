////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// HPsi.C
//
////////////////////////////////////////////////////////////////////////////////

#include "HPsi.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "Species.h"
#include "Matrix.h"
#include "PrintMem.h"
#include "blas.h"
#include <iomanip>
using namespace std;

#if USE_MASSV
extern "C" void vsincos(double *x, double *y, double *z, int *n);
#endif

////////////////////////////////////////////////////////////////////////////////
HPsi::HPsi(const Wavefunction& wfi, const AtomSet& as, const ChargeDensity& cd, vector<vector<double> >& v_r) : atoms_(as),vofr_(v_r)
{
   // wfi just used to initialize data structures with the correct basis
   // plane wave coefficients don't matter
   init(wfi);
   cell_moved(wfi);

  // FT for interpolation of wavefunctions on the fine grid
  ft.resize(wfi.nspin());
  for ( int ispin = 0; ispin < wfi.nspin(); ispin++ )
    ft[ispin].resize(wfi.nkp());

  for ( int ispin = 0; ispin < wfi.nspin(); ispin++ )
    for ( int ikp = 0; ikp < wfi.nkp(); ikp++ ) 
      ft[ispin][ikp] = cd.ft(ispin,ikp);

}

////////////////////////////////////////////////////////////////////////////////
HPsi::~HPsi(void) {
}

////////////////////////////////////////////////////////////////////////////////
void HPsi::init(const Wavefunction& wfi) {

   // max deviation from locality at upper quadrature point rcutloc[is]
   const double epsilon = 1.e-4;

   nsp = atoms_.nsp();
  
   lmax.resize(nsp);
   lloc.resize(nsp);
   lproj.resize(nsp);
   na.resize(nsp);
   npr.resize(nsp);
   nprna.resize(nsp);
   wt.resize(nsp);  
   nquad.resize(nsp);
   rquad.resize(nsp);
   wquad.resize(nsp);

   nspnl = 0;
   for ( int is = 0; is < nsp; is++ ) {
      Species *s = atoms_.species_list[is];
      npr[is] = 0;
      nprna[is] = 0;
    
      na[is] = atoms_.na(is);
      if ( s->non_local() && !s->ultrasoft()) {
         nspnl++;
         lmax[is] = s->lmax();
         lloc[is] = s->llocal();
         nquad[is] = s->nquad();
      
         // compute total number of projectors npr[is]
         // KB potentials have nlm projectors
         // semilocal potentials have nlm*nquad projectors
         if ( s->nquad() == 0 ) {
            npr[is] = s->nlm();
         }
         else {
            npr[is] = s->nlm() * nquad[is];
         }
         nprna[is] = npr[is] * na[is];
      
         // l value for projector ipr
         lproj[is].resize(npr[is]);
      
         // quadrature abcissae and weights
         rquad[is].resize(nquad[is]);
         wquad[is].resize(nquad[is]);
      
         enum quadrature_rule_type { TRAPEZOID, MODIF_TRAPEZOID, 
                                     TRAPEZOID_WITH_RIGHT_ENDPOINT, 
                                     SIMPSON };
      
         const quadrature_rule_type quad_rule = TRAPEZOID;
         //const quadrature_rule_type quad_rule = MODIF_TRAPEZOID;
         //const quadrature_rule_type quad_rule = TRAPEZOID_WITH_RIGHT_ENDPOINT;
         //const quadrature_rule_type quad_rule = SIMPSON;

         if ( quad_rule == TRAPEZOID ) {
            // trapezoidal rule with interior points only
            // (end points are zero)
            const double h = s->rquad() / (nquad[is]+1);
            for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
               rquad[is][iquad] = (iquad+1) * h;
               wquad[is][iquad] = h;
            }
            //cout << " HPsi::init: trapezoidal rule (interior)"
            //     << endl;
         }
         else if ( quad_rule == MODIF_TRAPEZOID ) {
            // use modified trapezoidal rule with interior points, and include
            // correction for first derivative at r=0 as
            // h^2/12 f'(0) where f'(0) is estimated with f(h)/h
            // i.e. add correction h/12) * f(h)
            // See Davis & Rabinowitz, p. 132
            const double h = s->rquad() / (nquad[is]+1);
            for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
               rquad[is][iquad] = (iquad+1) * h;
               wquad[is][iquad] = h;
            }
            wquad[is][0] += h / 12.0;
            //cout << " HPsi::init: modified trapezoidal rule"
            //     << endl;
         }
         else if ( quad_rule == TRAPEZOID_WITH_RIGHT_ENDPOINT ) {
            const double h = s->rquad() / nquad[is];
            for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
               rquad[is][iquad] = (iquad+1) * h;
               wquad[is][iquad] = h;
            }
            wquad[is][nquad[is]-1] = 0.5 * h;
            //cout << " HPsi::init: trapezoidal rule with right endpoint"
            //     << endl;
         }
         else if ( quad_rule == SIMPSON ) {
            // must have 2n+1 points
            assert(nquad[is]%2==1);
            const double h = s->rquad() / (nquad[is]-1);
            for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
               rquad[is][iquad] = iquad * h;
               if ( ( iquad == 0 ) || ( iquad == nquad[is]-1 ) )
                  wquad[is][iquad] = h / 3.0;
               else if ( iquad % 2 == 0 )
                  wquad[is][iquad] = h * 2.0 / 3.0;
               else
                  wquad[is][iquad] = h * 4.0 / 3.0;
            }
            //cout << " HPsi::init: Simpson rule" << endl;
         }
         else {
            assert(false);
         }
      
         // compute weights wt[is][ipr]
         wt[is].resize(npr[is]);
      
         // compute lproj[is][ipr] 
         int ipr_base = 0;
         for ( int l = 0; l <= lmax[is]; l++ ) {
            if ( l != lloc[is] ) {
               if ( nquad[is] == 0 ) {
                  // Kleinman-Bylander form
                  // wt[is][ipr]
                  // index = ipr_base+m
                  for ( int m = 0; m < 2*l+1; m++ ) {
                     const int ipr = ipr_base + m;
                     wt[is][ipr] = s->wsg(l);
                     lproj[is][ipr] = l;
                  }
                  ipr_base += 2*l+1;
               }
               else {
                  for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                     const double r = rquad[is][iquad];
                     double v,dv,vl,dvl;
                     s->dvpsr(l,r,v,dv); 
                     s->dvpsr(lloc[is],r,vl,dvl); 
                     // wt[is][iquad+ipr*nquad]
                     for ( int m = 0; m < 2*l+1; m++ ) {
                        const int ipr = ipr_base + iquad + nquad[is] * m;
                        wt[is][ipr] = ( v - vl ) * wquad[is][iquad];
                        lproj[is][ipr] = l;
                     }
                  }
                  ipr_base += (2*l+1) * nquad[is];
               }
            }
         }
         assert(ipr_base==npr[is]);
      } // if s->non_local()
   }

   twnl.resize(wfi.nkptloc());
   int ispin = 0; // basis sizes should agree between spin channels
   if (wfi.spinactive(ispin)) {
      for ( int kloc=0; kloc<wfi.nkptloc(); kloc++) {
         int ikp = wfi.kptloc(kloc);         // global index of local kpoint
         if (wfi.kptactive(ikp)) {
            assert(wfi.sd(ispin,ikp) != 0);
            const Basis& basis_ = wfi.sd(ispin,ikp)->basis();
            const int ngwl = basis_.localsize();

            twnl[kloc].resize(nsp);
            for ( int is = 0; is < nsp; is++ ) {
               Species *s = atoms_.species_list[is];
               if ( s->non_local() && !s->ultrasoft()) 
                  twnl[kloc][is].resize(npr[is]*ngwl);      
            }
         }
      }
   }
   
}

////////////////////////////////////////////////////////////////////////////////
void HPsi::cell_moved(const Wavefunction& wfi) {
   // update array twnl[ikp][is][ipr][ig]
   // following a change of cell dimensions
   // It is assumed that basis_ has been updated
   // It is assumed that nsp, npr[is], nquad[is] did not change since init

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

   int ispin = 0; // basis sizes should agree between spin channels
   if (wfi.spinactive(ispin)) {
      for ( int kloc=0; kloc<wfi.nkptloc(); kloc++) {
         int ikp = wfi.kptloc(kloc);         // global index of local kpoint
         if (wfi.kptactive(ikp)) {
            assert(wfi.sd(ispin,ikp) != 0);

            const Basis& basis_ = wfi.sd(ispin,ikp)->basis();
            const int ngwl = basis_.localsize();
   
            const double *kpg   = basis_.kpg_ptr();
            const double *kpg2  = basis_.kpg2_ptr();
            const double *kpgi  = basis_.kpgi_ptr();
            const double *kpg_x = basis_.kpgx_ptr(0);
            const double *kpg_y = basis_.kpgx_ptr(1);
            const double *kpg_z = basis_.kpgx_ptr(2);

            // compute twnl
            for ( int is = 0; is < atoms_.nsp(); is++ ) {
               Species *s = atoms_.species_list[is];
               int ilm = 0;
               for ( int l = 0; l <= lmax[is]; l++ ) {      
                  if ( l != lloc[is] ) {
                     if ( l == 0 ) {
                        if ( nquad[is] == 0 ) {
                           // Kleinman-Bylander
              
                           // twnl[kloc][is][ipr][ig]
                           // ipr = ilm = 0
                           // index = ig + ngwl*ipr, i.e. index = ig
                           double *t0  = &twnl[kloc][is][0];
                           for ( int ig = 0; ig < ngwl; ig++ ) {
                              double v;
                              s->vnlg(0,kpg[ig],v);
                              t0[ig] = s14pi * v;
                           }
                        }
                        else {
                           // semi-local
                           for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                              // twnl[kloc][is][ipr][ig]
                              // ipr = iquad + nquad[is]*ilm, where ilm=0
                              //     = iquad
                              // index = ig + ngwl*iquad
                              double *t0 = &twnl[kloc][is][ngwl*iquad];
                              const double r = rquad[is][iquad];

                              for ( int ig = 0; ig < ngwl; ig++ ) {
                                 // I(l=0) = 4 pi j_l(G r) r
                                 // twnl[kloc][is][ipr][l][ig] = 4 pi j_0(Gr_i) r_i Ylm
                                 // j_0(Gr) * r = sin(Gr) / G
                                 // Ylm = s14pi
                                 const double arg = kpg[ig] * r;
                                 // Note: for G=0, gi[0] = 0
                  
                                 const double tgi = kpgi[ig];
                                 const double ts = sin(arg);
                                 const double tc = cos(arg);
              
                                 t0[ig] = fpi * s14pi * ts * tgi;

                                 //ewd need to correct case when k+G = 0,
                                 //ewd i.e. sin((k+G)r)/(k+G)r -> 1, not 0
                                 if (kpgi[ig] == 0.0 && kpg[ig] == 0.0) 
                                    t0[ig] = fpi * s14pi * r;
                  
                              }
                           }
                        }
                        ilm += 2*l+1;
                     }
                     else if ( l == 1 ) {
                        if ( nquad[is] == 0 ) {
                           // Kleinman-Bylander
            
                           // twnl[kloc][is][ipr][ig]
                           // ipr = ilm
                           const int ipr1 = ilm;
                           const int ipr2 = ilm+1;
                           const int ipr3 = ilm+2;
                           // index = ig + ngwl*ilm
                           double *t1 = &twnl[kloc][is][ngwl*ipr1];
                           double *t2 = &twnl[kloc][is][ngwl*ipr2];
                           double *t3 = &twnl[kloc][is][ngwl*ipr3];
              
                           for ( int ig = 0; ig < ngwl; ig++ ) {
                              double v;
                              const double tg = kpg[ig];
                              s->vnlg(l,tg,v);
              
                              const double tgx = kpg_x[ig];
                              const double tgy = kpg_y[ig];
                              const double tgz = kpg_z[ig];
                              const double tgi = kpgi[ig];
              
                              const double y1 = s34pi * tgx * tgi;
                              const double y2 = s34pi * tgy * tgi;
                              const double y3 = s34pi * tgz * tgi;
              
                              t1[ig]  = y1 * v;
                              t2[ig]  = y2 * v;
                              t3[ig]  = y3 * v;
                           }
                        }
                        else {
                           // semi-local
                           for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                              // twnl[kloc][is][ipr][ig]
                              // index = ig + ngwl*(iquad+nquad[is]*ilm)
                              const int ipr1 = iquad+nquad[is]*ilm;
                              const int ipr2 = iquad+nquad[is]*(ilm+1);
                              const int ipr3 = iquad+nquad[is]*(ilm+2);
                              double *t1 = &twnl[kloc][is][ngwl*ipr1];
                              double *t2 = &twnl[kloc][is][ngwl*ipr2];
                              double *t3 = &twnl[kloc][is][ngwl*ipr3];
            
                              const double r = rquad[is][iquad];
                              for ( int ig = 0; ig < ngwl; ig++ ) {
                                 double v = 0.0, dv = 0.0;
                                 // j_1(Gr) = (1/(Gr))*(sin(Gr)/(Gr)-cos(Gr))
                                 const double tg = kpg[ig];
                                 const double z = tg * r;
                                 if ( z != 0.0 ) {
                                    const double zi = 1.0 / z;
                                    const double c = cos(z);
                                    const double s = sin(z);
                                    const double j1 = ( s * zi - c ) * zi;
                                    // v = 4 pi j1(Gr) r
                                    v = fpi * j1 * r;
                                 }
              
                                 const double tgx = kpg_x[ig];
                                 const double tgy = kpg_y[ig];
                                 const double tgz = kpg_z[ig];
                                 const double tgi = kpgi[ig];
 
                                 const double y1 = s34pi * tgx * tgi;
                                 const double y2 = s34pi * tgy * tgi;
                                 const double y3 = s34pi * tgz * tgi;
 
                                 t1[ig]  = y1 * v;
                                 t2[ig]  = y2 * v;
                                 t3[ig]  = y3 * v;
                              } // ig
                           }
                        }
                        ilm += 2*l+1;
                     }
                     else if ( l == 2 ) {
                        if ( nquad[is] == 0 ) {
                           // Kleinman-Bylander
                           const int ipr4 = ilm;
                           const int ipr5 = ilm+1;
                           const int ipr6 = ilm+2;
                           const int ipr7 = ilm+3;
                           const int ipr8 = ilm+4;
            
                           double *t4 = &twnl[kloc][is][ngwl*ipr4];
                           double *t5 = &twnl[kloc][is][ngwl*ipr5];
                           double *t6 = &twnl[kloc][is][ngwl*ipr6];
                           double *t7 = &twnl[kloc][is][ngwl*ipr7];
                           double *t8 = &twnl[kloc][is][ngwl*ipr8];
            
                           for ( int ig = 0; ig < ngwl; ig++ ) {
                              double v;
                              const double tg = kpg[ig];
              
                              s->vnlg(l,tg,v);
              
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
              
                              t4[ig]  = y4 * v;
                              t5[ig]  = y5 * v;
                              t6[ig]  = y6 * v;
                              t7[ig]  = y7 * v;
                              t8[ig]  = y8 * v;

                           }
                        }
                        else {
                           // semi-local
                           for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                              // twnl[kloc][is][ipr][ig]
                              // ipr = iquad+nquad[is]*ilm
                              // index = ig + ngwl*ipr
                              const int ipr4 = iquad+nquad[is]*ilm;
                              const int ipr5 = iquad+nquad[is]*(ilm+1);
                              const int ipr6 = iquad+nquad[is]*(ilm+2);
                              const int ipr7 = iquad+nquad[is]*(ilm+3);
                              const int ipr8 = iquad+nquad[is]*(ilm+4);
            
                              double *t4 = &twnl[kloc][is][ngwl*ipr4];
                              double *t5 = &twnl[kloc][is][ngwl*ipr5];
                              double *t6 = &twnl[kloc][is][ngwl*ipr6];
                              double *t7 = &twnl[kloc][is][ngwl*ipr7];
                              double *t8 = &twnl[kloc][is][ngwl*ipr8];

                              const double r = rquad[is][iquad];
                              for ( int ig = 0; ig < ngwl; ig++ ) {
                                 double v = 0.0;
                                 // j_2(z) = (3/z^3-1/z) sin(z) - 3/z^2 cos(z)
                                 // j_2(z) = (1/z)*((3/z^2-1)*sin(z) - (3/z) cos(z) )
                
                                 const double tg = kpg[ig];
                                 const double z = tg * r;
                                 if ( z != 0.0 ) {
                                    const double zi = 1.0 / z;
                                    const double c = cos(z);
                                    const double s = sin(z);
                                    const double j2 = ((3.0*zi*zi - 1.0) * s - 3.0*zi * c ) * zi;
                                    const double z2 = z * z;
                                    // v = 4 pi j2(Gr) r
                                    v = fpi * j2 * r;
                                 }
              
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
 
                                 t4[ig]  = y4 * v;
                                 t5[ig]  = y5 * v;
                                 t6[ig]  = y6 * v;
                                 t7[ig]  = y7 * v;
                                 t8[ig]  = y8 * v;

                              } // ig
                           } // iquad
                        }
                        ilm += 2*l+1;
                     }
                     else if ( l == 3 ) {
                        if ( nquad[is] == 0 ) {
                           // Kleinman-Bylander
                           const int ipr9 = ilm;
                           const int ipr10 = ilm+1;
                           const int ipr11 = ilm+2;
                           const int ipr12 = ilm+3;
                           const int ipr13 = ilm+4;
                           const int ipr14 = ilm+5;
                           const int ipr15 = ilm+6;
            
                           double *t9 = &twnl[kloc][is][ngwl*ipr9];
                           double *t10 = &twnl[kloc][is][ngwl*ipr10];
                           double *t11 = &twnl[kloc][is][ngwl*ipr11];
                           double *t12 = &twnl[kloc][is][ngwl*ipr12];
                           double *t13 = &twnl[kloc][is][ngwl*ipr13];
                           double *t14 = &twnl[kloc][is][ngwl*ipr14];
                           double *t15 = &twnl[kloc][is][ngwl*ipr15];
            
                           for ( int ig = 0; ig < ngwl; ig++ ) {
                              double v,dv;
                              const double tg = kpg[ig];
              
                              s->vnlg(l,tg,v);
              
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
                              const double tgxy = tgx * tgy * tgi2;
                              const double tgyz = tgy * tgz * tgi2;
                              const double tgxz = tgx * tgz * tgi2;
              
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

                              t9[ig]  = y9 * v;
                              t10[ig]  = y10 * v;
                              t11[ig]  = y11 * v;
                              t12[ig]  = y12 * v;
                              t13[ig]  = y13 * v;
                              t14[ig]  = y14 * v;
                              t15[ig]  = y15 * v;
                           }
                        }
                        else {
                           // semi-local
                           for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                              // twnl[kloc][is][ipr][ig]
                              // ipr = iquad+nquad[is]*ilm
                              // index = ig + ngwl*ipr
                              const int ipr9 = iquad+nquad[is]*ilm;
                              const int ipr10 = iquad+nquad[is]*(ilm+1);
                              const int ipr11 = iquad+nquad[is]*(ilm+2);
                              const int ipr12 = iquad+nquad[is]*(ilm+3);
                              const int ipr13 = iquad+nquad[is]*(ilm+4);
                              const int ipr14 = iquad+nquad[is]*(ilm+5);
                              const int ipr15 = iquad+nquad[is]*(ilm+6);
            
                              double *t9 = &twnl[kloc][is][ngwl*ipr9];
                              double *t10 = &twnl[kloc][is][ngwl*ipr10];
                              double *t11 = &twnl[kloc][is][ngwl*ipr11];
                              double *t12 = &twnl[kloc][is][ngwl*ipr12];
                              double *t13 = &twnl[kloc][is][ngwl*ipr13];
                              double *t14 = &twnl[kloc][is][ngwl*ipr14];
                              double *t15 = &twnl[kloc][is][ngwl*ipr15];

                              const double r = rquad[is][iquad];
                              for ( int ig = 0; ig < ngwl; ig++ ) {
                                 double v = 0.0, dv = 0.0;
                                 // j_3(z) = (15/z^4-6/z^2) sin(z) - (15/z^3 - 1/z) cos(z)
                
                                 const double tg = kpg[ig];
                                 const double z = tg * r;
                                 if ( z != 0.0 ) {
                                    const double zi = 1.0 / z;
                                    const double c = cos(z);
                                    const double s = sin(z);
                                    const double j3 = (15.0*zi*zi - 6.0) * zi*zi * s - (15.0*zi*zi-1.0) * zi * c;
                                    // v = 4 pi j3(Gr) r
                                    v = fpi * j3 * r;
                                 }
              
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
                                 const double tgxy = tgx * tgy * tgi2;
                                 const double tgyz = tgy * tgz * tgi2;
                                 const double tgxz = tgx * tgz * tgi2;
                                 const double tgx3 = tgx2*tgx * tgi3;
                                 const double tgy3 = tgy2*tgy * tgi3;
                                 const double tgz3 = tgz2*tgz * tgi3;
              
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
                
                                 t9[ig]  = y9 * v;
                                 t10[ig]  = y10 * v;
                                 t11[ig]  = y11 * v;
                                 t12[ig]  = y12 * v;
                                 t13[ig]  = y13 * v;
                                 t14[ig]  = y14 * v;
                                 t15[ig]  = y15 * v;

                              } // ig
                           } // iquad
                        }
                        ilm += 2*l+1;
                     }
                     else {
                        assert(false);
                     }
                  } // l != lloc[is]
               } // for l
               assert(ilm == s->nlm());
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void HPsi::compute(const Wavefunction& wf, Wavefunction& dwf)
{
   // apply Hamiltonian operator to input Wavefunction wf, store in dwf
   //
   // note:  this is current a hacky implementation, just to let us test some new
   // eigensolvers and preconditioners.  This won't work with ultrasoft or DFT+U.

   
   // define atom block size
   //const int na_block_size = 32;
   const int na_block_size = 256;
   vector<vector<double> > tau;
   atoms_.get_positions(tau,true);

   for ( int ispin = 0; ispin < wf.nspin(); ispin++ ) {
      if (wf.spinactive(ispin)) {
         for ( int kloc=0; kloc<wf.nkptloc(); kloc++) {
            assert(wf.nkptloc() == dwf.nkptloc());
            int ikp = wf.kptloc(kloc);         // global index of local kpoint
            if (wf.kptactive(ikp)) {
               assert(wf.sd(ispin,ikp) != 0);
               assert(dwf.sd(ispin,ikp) != 0);

               const SlaterDet* sdp = wf.sd(ispin,ikp);
               SlaterDet* dsdp = dwf.sd(ispin,ikp);
               const Basis& basis_ = sdp->basis();
               const vector<double>& occ = sdp->occ();
               const vector<double>& eig = sdp->eig();
               const int ngwl = basis_.localsize();
               const Context& ctxt_ = sdp->context();
               
               const ComplexMatrix& c = sdp->c();
               ComplexMatrix& cp = dsdp->c();
               const int mloc = c.mloc();
               assert(mloc == cp.mloc());
               
               cp.clear();
               double tsum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  
               const double omega = basis_.cell().volume();
               assert(omega != 0.0);
               const double omega_inv = 1.0 / omega;

               // regular norm-conserving nonlocal potential projection
               for ( int is = 0; is < nsp; is++ ) {  
                  Species *s = atoms_.species_list[is];
                  if (npr[is] > 0 ) { // species is is non-local, norm-conserving

                     valarray<double> kpgr(na_block_size*ngwl); // kpgr[ig+ia*ngwl]
                     valarray<double> ckpgr(na_block_size*ngwl); // ckpgr[ig+ia*ngwl]
                     valarray<double> skpgr(na_block_size*ngwl); // skpgr[ig+ia*ngwl]

                     valarray<double> tmpfion(3*na[is]);      
                     tmpfion = 0.0;                           
                     // define number of atom blocks                                        
                     const int na_blocks = na[is] / na_block_size +                         
                         ( na[is] % na_block_size == 0 ? 0 : 1 );         

                     valarray<double> anl_loc_gamma;
                     valarray<complex<double> > anl_loc;
                     //if (basis_.real())
                     //  anl_loc_gamma.resize(npr[is]*na_block_size*2*ngwl);
                     //else 
                     //  anl_loc.resize(npr[is]*na_block_size*ngwl);
                     //vector<double> anl_loc_gamma;
                     //vector<complex<double> > anl_loc;
                     if (basis_.real())
                        anl_loc_gamma.resize(npr[is]*na_block_size*2*mloc,0.0);
                     else 
                        anl_loc.resize(npr[is]*na_block_size*mloc,complex<double>(0.0,0.0));
        
                     const int nstloc = sdp->nstloc();
                     // fnl_loc[ipra][n]
                     valarray<double> fnl_loc_gamma;
                     valarray<double> fnl_buf_gamma;
                     valarray<complex<double> > fnl_loc;
                     valarray<complex<double> > fnl_buf;
                     if (basis_.real()) {
                        fnl_loc_gamma.resize(npr[is]*na_block_size*nstloc);
                        fnl_buf_gamma.resize(npr[is]*na_block_size*nstloc);
                     }
                     else {
                        fnl_loc.resize(npr[is]*na_block_size*nstloc);
                        fnl_buf.resize(npr[is]*na_block_size*nstloc);
                     }
                     for ( int ia_block = 0; ia_block < na_blocks; ia_block++ ) {
                        // process projectors of atoms in block ia_block             

                        const int iastart = ia_block * na_block_size;                
                        const int iaend = (ia_block+1) * na_block_size < na[is] ?    
                            (ia_block+1) * na_block_size :               
                            na[is];                                      
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
#if USE_MASSV  
                        vsincos(&skpgr[0],&ckpgr[0],&kpgr[0],&len);
#else
#pragma omp parallel for
                        for ( int i = 0; i < len; i++ ) {
                           const double arg = kpgr[i];
                           skpgr[i] = sin(arg);                                                 
                           ckpgr[i] = cos(arg);                                                 
                        }                                                                    
#endif

                        // compute anl_loc                                                   
                        for ( int ipr = 0; ipr < npr[is]; ipr++ ) {
                           // twnl[kloc][is][ig+ngwl*ipr]
                           const double * t = &twnl[kloc][is][ngwl*ipr];                            
                           const int l = lproj[is][ipr];                                      

                           // anl_loc[ig+ipra*ngwl]                                           
                           for ( int ia = 0; ia < ia_block_size; ia++ ) {
                              double* a;
                              if (basis_.real())
                                 a = &anl_loc_gamma[2*(ia+ipr*ia_block_size)*ngwl];             
                              else
                                 a = (double*)&anl_loc[(ia+ipr*ia_block_size)*ngwl];             
                              const double* c = &ckpgr[ia*ngwl];
                              const double* s = &skpgr[ia*ngwl];
                              if ( l == 0 ) {
#pragma omp parallel for
                                 for ( int ig = 0; ig < ngwl; ig++ ) {
                                    a[2*ig]   = t[ig] * c[ig];                                   
                                    a[2*ig+1] = t[ig] * s[ig];                                   
                                 }
                              }
                              else if ( l == 1 ) {
                                 /* Next line: -i * eigr */                                   
                                 /* -i * (a+i*b) = b - i*a */                                 
#pragma omp parallel for
                                 for ( int ig = 0; ig < ngwl; ig++ ) {
                                    a[2*ig]   =  t[ig] * s[ig];
                                    a[2*ig+1] = -t[ig] * c[ig];
                                 }
                              }
                              else if ( l == 2 ) {
                                 // Next line: (-) sign for -eigr                             
                                 /* -1 * (a+i*b) = -a - i*b */                                 
#pragma omp parallel for
                                 for ( int ig = 0; ig < ngwl; ig++ ) {
                                    a[2*ig]   = -t[ig] * c[ig];
                                    a[2*ig+1] = -t[ig] * s[ig];
                                 }
                              }
                              else if ( l == 3 ) {
                                 /* i * (a+i*b) = -b + i*a */ 
#pragma omp parallel for
                                 for ( int ig = 0; ig < ngwl; ig++ ) {
                                    a[2*ig]   = -t[ig] * s[ig];
                                    a[2*ig+1] = t[ig] * c[ig];
                                 }
                              }
                           }
                        } // ipr                                                             
                                                                             
                        // array anl_loc is complete                                         
                                                                             
                        // compute fnl[npra][nstloc] = anl^T * c                             
                        double one=1.0;                                                      
                        char ct='t';                                                         
                        int twongwl = 2 * ngwl;                                              
                        int twomloc = 2 * mloc;                                              
                        int nprnaloc = ia_block_size * npr[is];                              
                        int c_lda;
                        const complex<double>* pc = c.cvalptr();

                        if (basis_.real()) {
                           c_lda = 2*sdp->c().mloc();                                        
                           dgemm(&ct,&cn,&nprnaloc,(int*)&nstloc,&twongwl,&one,
                                 &anl_loc_gamma[0],&twongwl, (double*)pc, &c_lda,
                                 &zero,&fnl_loc_gamma[0],&nprnaloc);
                        }
                        else {
                           c_lda = mloc;
                           complex<double> zzero = complex<double>(0.0,0.0);
                           complex<double> zone = complex<double>(1.0,0.0);
                           char cc='c';
          
                           zgemm(&cc,&cn,&nprnaloc,(int*)&nstloc,(int*)&ngwl,&zone,
                                 &anl_loc[0],(int *)&ngwl, (complex<double> *)pc, &c_lda,
                                 &zzero,&fnl_loc[0],&nprnaloc);
                        }
          
                        // for k=(0,0,0):  need to correct for double counting of G=0, i.e. if ctxt_.myrow() == 0
                        if (basis_.real()) { 
                           if ( ctxt_.myrow() == 0 ) {
                              // rank-one update                                                 
                              // dger(m,n,alpha,x,incx,y,incy,a,lda);                            
                              // a += alpha * x * transpose(y)                                   
                              // x = first row of anl_loc                                        
                              // y^T = first row of c                                            
                              double alpha = -0.5;                                               
                              dger(&nprnaloc,(int*)&nstloc,&alpha,&anl_loc_gamma[0],&twongwl,
                                   (double*)pc,&c_lda,&fnl_loc_gamma[0],&nprnaloc);
                           }
                        }

                        if (nstloc > 0) {
                           // Allreduce fnl partial sum
                           MPI_Comm basis_comm = basis_.context().comm();                       
                           int fnl_size = nprnaloc*nstloc;                                   
                           if (basis_.real()) {
                              MPI_Allreduce(&fnl_loc_gamma[0],&fnl_buf_gamma[0],fnl_size,                      
                                            MPI_DOUBLE,MPI_SUM,basis_comm);                        
                           }
                           else {
                              MPI_Allreduce(&fnl_loc[0],&fnl_buf[0],fnl_size,                      
                                            MPI_DOUBLE_COMPLEX,MPI_SUM,basis_comm);                        
                           }
                        }
                      
                        // factor 2.0 in next line is: counting G, -G                        
                        if (basis_.real())
                           fnl_loc_gamma = 2.0 * fnl_buf_gamma;
                        else
                           fnl_loc = fnl_buf;
                                
                        // compute fnl
                        if (basis_.real()) {
                           for ( int ipr = 0; ipr < npr[is]; ipr++ ) { 
                              const double fac = wt[is][ipr] * omega_inv;                        
                              for ( int lj=0; lj < sdp->c().nblocks(); lj++ )
                              {
                                 for ( int jj=0; jj < sdp->c().nbs(lj); jj++ )
                                 {
                                    // global state index
                                    const int nglobal = sdp->c().j(lj,jj);
                                    const int norig = lj*sdp->c().nb()+jj;
                                    const double facn = fac * occ[nglobal];
                                    for ( int ia = 0; ia < ia_block_size; ia++ ) {
                                       const int i = ia + ipr*ia_block_size + norig * nprnaloc;           
                                       const double tmp = fnl_loc_gamma[i];
                                       fnl_loc_gamma[i] = fac * tmp;
                                    }
                                 }
                              }
                           }
                        }
                        else {
                           for ( int ipr = 0; ipr < npr[is]; ipr++ ) { 
                              const double fac = wt[is][ipr] * omega_inv;                        
                              for ( int lj=0; lj < sdp->c().nblocks(); lj++ )
                              {
                                 for ( int jj=0; jj < sdp->c().nbs(lj); jj++ )
                                 {
                                    // global state index
                                    const int nglobal = sdp->c().j(lj,jj);
                                    const int norig = lj*sdp->c().nb()+jj;
                                    const double facn = fac * occ[nglobal];
                                    for ( int ia = 0; ia < ia_block_size; ia++ ) {
                                       const int i = ia + ipr*ia_block_size + norig * nprnaloc;           
                                       const complex<double> tmp = fnl_loc[i];                                 
                                       fnl_loc[i] = fac * tmp;
                                    }
                                 }
                              }
                           }
                        }

                        // compute cp += anl * fnl                                         
                        complex<double>* pcp = cp.valptr();                            

                        int cp_lda;
                        if (basis_.real()) {
                           cp_lda = 2*dsdp->c().mloc();
                           dgemm(&cn,&cn,&twongwl,(int*)&nstloc,&nprnaloc,&one,
                                 &anl_loc_gamma[0],&twongwl, &fnl_loc_gamma[0],&nprnaloc,
                                 &one,(double*)pcp, &cp_lda);
                        }
                        else {
                           int cp_lda = dsdp->c().mloc();
                           complex<double> zone = complex<double>(1.0,0.0);
                           zgemm(&cn,&cn,(int*)&ngwl,(int*)&nstloc,&nprnaloc,&zone,
                                 &anl_loc[0],(int*)&ngwl, &fnl_loc[0],&nprnaloc,
                                 &zone,(complex<double>*)pcp, &cp_lda);
                        }
                     } // ia_block              
                  } // npr[is]>0
               } // is

               // nonlocal contribution finished, add other terms

               
               // Laplacian
               const double* kpg2 = basis_.kpg2_ptr();
               const int ngwloc = basis_.localsize();
               for ( int n = 0; n < sdp->nstloc(); n++ ) {
                  for ( int ig = 0; ig < ngwloc; ig++ ) {
                     cp[ig+mloc*n] += 0.5 * kpg2[ig] * c[ig+mloc*n];
                  }
               }
               sdp->rs_mul_add(*ft[ispin][ikp], &vofr_[ispin][0], *dsdp);

            } // if kptactive
         } // kloc
      } // if spinactive
   } // ispin

   return;

}
