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
// NonLocalPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "NonLocalPotential.h"
#include "Species.h"
#include "Matrix.h"
#include "PrintMem.h"
#include "blas.h"
#include <iomanip>
using namespace std;

#if HAVE_MASSV
extern "C" void vsincos(double *x, double *y, double *z, int *n);
#endif
extern "C" void myzdotc(const int size, complex<double>* v1, complex<double>* v2, complex<double>* vout);


////////////////////////////////////////////////////////////////////////////////
NonLocalPotential::~NonLocalPotential(void) {
  //for ( int is = 0; is < anl.size(); is++ )
  //{
  //  delete anl[is];
  //}
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::print_timing(void) {
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    uint64_t count = (*i).second.counts();
    //if ( ctxt_.myproc()==0 ) {
    if ( ctxt_.mype()==0 ) {
       cout << left << setw(34) << "<timing where=\"nonlocal\""
            << setw(8) << " name=\""
            << setw(15) << (*i).first << "\""
            << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
            << " max=\"" << setprecision(3) << setw(9) << tmax << "\""
            << " count=\"" << setw(9) << count << "\"/>"
            << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::init(const bool compute_stress) {

  // max deviation from locality at upper quadrature point rcutloc[is]
  const double epsilon = 1.e-4;

  const int ngwl = basis_.localsize();
  
  nsp = atoms_.nsp();
  
  lmax.resize(nsp);
  lloc.resize(nsp);
  lproj.resize(nsp);
  icproj.resize(nsp);
  iprojlm.resize(nsp);
  na.resize(nsp);
  npr.resize(nsp);
  nprna.resize(nsp);
  wt.resize(nsp);
  twnl.resize(nsp);
  dtwnl.resize(nsp);
  
  nquad.resize(nsp);
  rquad.resize(nsp);
  wquad.resize(nsp);
   
  nspnl = 0;
  ultrasoft_ = false;
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_.species_list[is];
    if (s->ultrasoft()) ultrasoft_ = true;
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
      icproj[is].resize(npr[is]);
      iprojlm[is].resize(s->lmax() + 1);
      
      twnl[is].resize(npr[is]*ngwl);      

      if (compute_stress)
         dtwnl[is].resize(npr[is]*6*ngwl);
      
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
        //cout << " NonLocalPotential::init: trapezoidal rule (interior)"
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
        //cout << " NonLocalPotential::init: modified trapezoidal rule"
        //     << endl;
      }
      else if ( quad_rule == TRAPEZOID_WITH_RIGHT_ENDPOINT ) {
        const double h = s->rquad() / nquad[is];
        for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
          rquad[is][iquad] = (iquad+1) * h;
          wquad[is][iquad] = h;
        }
        wquad[is][nquad[is]-1] = 0.5 * h;
        //cout << " NonLocalPotential::init: trapezoidal rule with right endpoint"
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
        //cout << " NonLocalPotential::init: Simpson rule" << endl;
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

	  iprojlm[is][l].resize(2*l + 1);
	  for ( int m = 0; m < 2*l+1; m++ ) {
	    iprojlm[is][l][m].resize(s->nchannels());
	  }
	  
          if ( nquad[is] == 0 ) {
	    for(int ic = 0; ic < s->nchannels(); ic++){
	      // Kleinman-Bylander form
	      // wt[is][ipr]
	      // index = ipr_base+m
	      for ( int m = 0; m < 2*l+1; m++ ) {
		const int ipr = ipr_base + m;
		wt[is][ipr] = s->wsg(l, ic);
		lproj[is][ipr] = l;
		icproj[is][ipr] = ic;
		iprojlm[is][l][m][ic] = ipr;
	      }
	      ipr_base += 2*l+1;
	    }
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
		iprojlm[is][l][m][0] = ipr;
              }
            }
            ipr_base += (2*l+1) * nquad[is];
          }
        }
      }
      assert(ipr_base==npr[is]);
    } // if s->non_local()
  }

  highmem_ = false;

}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::update_twnl(const bool compute_stress) {
  // update arrays twnl[is][ipr][ig], dtwnl[is][ipr][j][ig],
  // following a change of cell dimensions
  // It is assumed that basis_ has been updated
  // It is assumed that nsp, npr[is], nquad[is] did not change since init
  
  tmap["update_twnl"].start();
    
  const int ngwl = basis_.localsize();
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

  const double *kpg   = basis_.kpg_ptr();
  const double *kpg2  = basis_.kpg2_ptr();
  const double *kpgi  = basis_.kpgi_ptr();
  const double *kpg_x = basis_.kpgx_ptr(0);
  const double *kpg_y = basis_.kpgx_ptr(1);
  const double *kpg_z = basis_.kpgx_ptr(2);

  // compute twnl and dtwnl
  for ( int is = 0; is < nsp; is++ ) {
    Species *s = atoms_.species_list[is];
    if (!s->ultrasoft()) {
      int ilm = 0;
      for ( int l = 0; l <= lmax[is]; l++ ) {      
        if ( l != lloc[is] ) {
          if ( l == 0 ) {
            if ( nquad[is] == 0 ) {
	      for(int ic = 0; ic < s->nchannels(); ic++){
		// Kleinman-Bylander
		
		// twnl[is][ipr][ig]
		// ipr = ilm = 0
		// index = ig + ngwl*ipr, i.e. index = ig
		double *t0  = &twnl[is][iprojlm[is][l][0][ic]*ngwl];
            
		// dtwnl[is][ipr][ij][ngwl]
		// index = ig + ngwl * ( ij + 6 * ipr ), ipr = 0
		// i.e. index = ig + ij * ngwl
		double *dt0_xx,*dt0_yy,*dt0_zz,*dt0_xy,*dt0_yz,*dt0_xz;
		if (compute_stress)
		  {
		    dt0_xx = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 0*ngwl];
		    dt0_yy = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 1*ngwl];
		    dt0_zz = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 2*ngwl];
		    dt0_xy = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 3*ngwl];
		    dt0_yz = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 4*ngwl];
		    dt0_xz = &dtwnl[is][ngwl*6*iprojlm[is][l][0][ic] + 5*ngwl];
		  }
              
		for ( int ig = 0; ig < ngwl; ig++ ) {
		  double v,dv;
		  s->dvnlg(l, ic, kpg[ig], v, dv);

		  t0[ig] = s14pi * v;

		  if (compute_stress)
		    {
		      const double tgx = kpg_x[ig];
		      const double tgy = kpg_y[ig];
		      const double tgz = kpg_z[ig];
		      const double tgx2 = tgx * tgx;
		      const double tgy2 = tgy * tgy;
		      const double tgz2 = tgz * tgz;
              
		      const double tmp = kpgi[ig] * s14pi * dv;
		      dt0_xx[ig] = tmp * tgx * tgx;
		      dt0_yy[ig] = tmp * tgy * tgy;
		      dt0_zz[ig] = tmp * tgz * tgz;
		      dt0_xy[ig] = tmp * tgx * tgy;
		      dt0_yz[ig] = tmp * tgy * tgz;
		      dt0_xz[ig] = tmp * tgx * tgz;
		    }
		}
	      }
            }
            else {
              // semi-local
              for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                // twnl[is][ipr][ig]
                // ipr = iquad + nquad[is]*ilm, where ilm=0
                //     = iquad
                // index = ig + ngwl*iquad
                double *t0 = &twnl[is][ngwl*iquad];
                // dtwnl[is][ipr][j][ngwl]
                // index = ig + ngwl * ( ij + 6 * iquad)
                double *dt0_xx,*dt0_yy,*dt0_zz,*dt0_xy,*dt0_yz,*dt0_xz;
                if (compute_stress)
                {
                   dt0_xx = &dtwnl[is][ngwl*(0+6*iquad)];
                   dt0_yy = &dtwnl[is][ngwl*(1+6*iquad)];
                   dt0_zz = &dtwnl[is][ngwl*(2+6*iquad)];
                   dt0_xy = &dtwnl[is][ngwl*(3+6*iquad)];
                   dt0_yz = &dtwnl[is][ngwl*(4+6*iquad)];
                   dt0_xz = &dtwnl[is][ngwl*(5+6*iquad)];
                }
                const double r = rquad[is][iquad];

                for ( int ig = 0; ig < ngwl; ig++ ) {
                  // I(l=0) = 4 pi j_l(G r) r
                  // twnl[is][ipr][l][ig] = 4 pi j_0(Gr_i) r_i Ylm
                  // j_0(Gr) * r = sin(Gr) / G
                  // Ylm = s14pi
                  const double arg = kpg[ig] * r;
                  // Note: for G=0, gi[0] = 0
                  
                  const double tgx = kpg_x[ig];
                  const double tgy = kpg_y[ig];
                  const double tgz = kpg_z[ig];
                  const double tgi = kpgi[ig];
                  const double tgi2 = tgi * tgi;
                
                  const double ts = sin(arg);
                  const double tc = cos(arg);
              
                  t0[ig] = fpi * s14pi * ts * tgi;

                  //ewd need to correct case when k+G = 0, i.e. sin((k+G)r)/(k+G)r -> 1, not 0
                  if (kpgi[ig] == 0.0 && kpg[ig] == 0.0) 
                    t0[ig] = fpi * s14pi * r;
                  
                  // dtwnl = fpi s14pi G_i G_j / G (r cos(Gr)/G -sin(Gr)/G^2)
                  if (compute_stress)
                  {
                     const double tmp = fpi * s14pi * tgi2 * (r*tc - ts*tgi);
                     dt0_xx[ig] = tmp * tgx * tgx;
                     dt0_yy[ig] = tmp * tgy * tgy;
                     dt0_zz[ig] = tmp * tgz * tgz;
                     dt0_xy[ig] = tmp * tgx * tgy;
                     dt0_yz[ig] = tmp * tgy * tgz;
                     dt0_xz[ig] = tmp * tgx * tgz;
                  }
                }
              }
            }
            ilm += 2*l+1;
          }
          else if ( l == 1 ) {
            if ( nquad[is] == 0 ) {
	      for(int ic = 0; ic < s->nchannels(); ic++){
		// Kleinman-Bylander
		
		// twnl[is][ipr][ig]
		const int ipr1 = iprojlm[is][l][0][ic];
		const int ipr2 = iprojlm[is][l][1][ic];
		const int ipr3 = iprojlm[is][l][2][ic];
		// index = ig + ngwl*ilm
		double *t1 = &twnl[is][ngwl*ipr1];
		double *t2 = &twnl[is][ngwl*ipr2];
		double *t3 = &twnl[is][ngwl*ipr3];
		
		// dtwnl[is][ipr][ij][ngwl]
		// index = ig + ngwl * ( ij + 6 * ipr )
		double *dt1_xx,*dt1_yy,*dt1_zz,*dt1_xy,*dt1_yz,*dt1_xz;
		double *dt2_xx,*dt2_yy,*dt2_zz,*dt2_xy,*dt2_yz,*dt2_xz;
		double *dt3_xx,*dt3_yy,*dt3_zz,*dt3_xy,*dt3_yz,*dt3_xz;
		if (compute_stress)
		  {
		    dt1_xx = &dtwnl[is][ngwl*(0+6*ipr1)];
		    dt1_yy = &dtwnl[is][ngwl*(1+6*ipr1)];
		    dt1_zz = &dtwnl[is][ngwl*(2+6*ipr1)];
		    dt1_xy = &dtwnl[is][ngwl*(3+6*ipr1)];
		    dt1_yz = &dtwnl[is][ngwl*(4+6*ipr1)];
		    dt1_xz = &dtwnl[is][ngwl*(5+6*ipr1)];

		    dt2_xx = &dtwnl[is][ngwl*(0+6*ipr2)];
		    dt2_yy = &dtwnl[is][ngwl*(1+6*ipr2)];
		    dt2_zz = &dtwnl[is][ngwl*(2+6*ipr2)];
		    dt2_xy = &dtwnl[is][ngwl*(3+6*ipr2)];
		    dt2_yz = &dtwnl[is][ngwl*(4+6*ipr2)];
		    dt2_xz = &dtwnl[is][ngwl*(5+6*ipr2)];
              
		    dt3_xx = &dtwnl[is][ngwl*(0+6*ipr3)];
		    dt3_yy = &dtwnl[is][ngwl*(1+6*ipr3)];
		    dt3_zz = &dtwnl[is][ngwl*(2+6*ipr3)];
		    dt3_xy = &dtwnl[is][ngwl*(3+6*ipr3)];
		    dt3_yz = &dtwnl[is][ngwl*(4+6*ipr3)];
		    dt3_xz = &dtwnl[is][ngwl*(5+6*ipr3)];
		  }
              
		for ( int ig = 0; ig < ngwl; ig++ ) {
		  double v,dv;
		  const double tg = kpg[ig];
		  s->dvnlg(l, ic, tg, v, dv);
              
		  const double tgx = kpg_x[ig];
		  const double tgy = kpg_y[ig];
		  const double tgz = kpg_z[ig];
		  const double tgx2 = tgx * tgx;
		  const double tgy2 = tgy * tgy;
		  const double tgz2 = tgz * tgz;
              
		  const double tgi = kpgi[ig];
		  const double tgi2 = tgi * tgi;
              
		  const double y1 = s34pi * tgx * tgi;
		  const double y2 = s34pi * tgy * tgi;
		  const double y3 = s34pi * tgz * tgi;
              
		  t1[ig]  = y1 * v;
		  t2[ig]  = y2 * v;
		  t3[ig]  = y3 * v;

		  if (compute_stress)
		    {
		      const double fac1 = - y1 * ( v - tg * dv ) * tgi2;
		      // m=x
		      dt1_xx[ig] = fac1 * tgx2 + v * y1;
		      dt1_yy[ig] = fac1 * tgy2;
		      dt1_zz[ig] = fac1 * tgz2;
		      dt1_xy[ig] = fac1 * tgx * tgy;
		      dt1_yz[ig] = fac1 * tgy * tgz;
		      dt1_xz[ig] = fac1 * tgx * tgz;
              
		      const double fac2 = - y2 * ( v - tg * dv ) * tgi2;
		      // m=y
		      dt2_xx[ig] = fac2 * tgx2;
		      dt2_yy[ig] = fac2 * tgy2 + v * y2;
		      dt2_zz[ig] = fac2 * tgz2;
		      dt2_xy[ig] = fac2 * tgx * tgy + v * y1;
		      dt2_yz[ig] = fac2 * tgy * tgz;
		      dt2_xz[ig] = fac2 * tgx * tgz;
              
		      const double fac3 = - y3 * ( v - tg * dv ) * tgi2;
		      // m=z
		      dt3_xx[ig] = fac3 * tgx2;
		      dt3_yy[ig] = fac3 * tgy2;
		      dt3_zz[ig] = fac3 * tgz2 + v * y3;
		      dt3_xy[ig] = fac3 * tgx * tgy;
		      dt3_yz[ig] = fac3 * tgy * tgz + v * y2;
		      dt3_xz[ig] = fac3 * tgx * tgz + v * y1;
		    }
		}
              }
            }
            else {
              // semi-local
              for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                // twnl[is][ipr][ig]
                // index = ig + ngwl*(iquad+nquad[is]*ilm)
                const int ipr1 = iquad+nquad[is]*ilm;
                const int ipr2 = iquad+nquad[is]*(ilm+1);
                const int ipr3 = iquad+nquad[is]*(ilm+2);
                double *t1 = &twnl[is][ngwl*ipr1];
                double *t2 = &twnl[is][ngwl*ipr2];
                double *t3 = &twnl[is][ngwl*ipr3];
            
                // dtwnl[is][ipr][j][ngwl]
                // index = ig + ngwl * ( ij + 6 * ipr )
                double *dt1_xx,*dt1_yy,*dt1_zz,*dt1_xy,*dt1_yz,*dt1_xz;
                double *dt2_xx,*dt2_yy,*dt2_zz,*dt2_xy,*dt2_yz,*dt2_xz;
                double *dt3_xx,*dt3_yy,*dt3_zz,*dt3_xy,*dt3_yz,*dt3_xz;
                if (compute_stress)
                {
                   dt1_xx = &dtwnl[is][ngwl*(0+6*ipr1)];
                   dt1_yy = &dtwnl[is][ngwl*(1+6*ipr1)];
                   dt1_zz = &dtwnl[is][ngwl*(2+6*ipr1)];
                   dt1_xy = &dtwnl[is][ngwl*(3+6*ipr1)];
                   dt1_yz = &dtwnl[is][ngwl*(4+6*ipr1)];
                   dt1_xz = &dtwnl[is][ngwl*(5+6*ipr1)];

                   dt2_xx = &dtwnl[is][ngwl*(0+6*ipr2)];
                   dt2_yy = &dtwnl[is][ngwl*(1+6*ipr2)];
                   dt2_zz = &dtwnl[is][ngwl*(2+6*ipr2)];
                   dt2_xy = &dtwnl[is][ngwl*(3+6*ipr2)];
                   dt2_yz = &dtwnl[is][ngwl*(4+6*ipr2)];
                   dt2_xz = &dtwnl[is][ngwl*(5+6*ipr2)];

                   dt3_xx = &dtwnl[is][ngwl*(0+6*ipr3)];
                   dt3_yy = &dtwnl[is][ngwl*(1+6*ipr3)];
                   dt3_zz = &dtwnl[is][ngwl*(2+6*ipr3)];
                   dt3_xy = &dtwnl[is][ngwl*(3+6*ipr3)];
                   dt3_yz = &dtwnl[is][ngwl*(4+6*ipr3)];
                   dt3_xz = &dtwnl[is][ngwl*(5+6*ipr3)];
                }
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
                    const double dj1 = 
                        ( 2.0 * z * c + ( z*z - 2.0 ) * s ) * zi*zi*zi;
                    // v = 4 pi j1(Gr) r
                    v = fpi * j1 * r;
                    // dv = d/dG v = 4 pi dj1(Gr)/d(Gr) d(Gr)/dG r
                    //    = 4 pi dj1 r^2
                    dv = fpi * dj1 * r * r;
                  }
              
                  const double tgx = kpg_x[ig];
                  const double tgy = kpg_y[ig];
                  const double tgz = kpg_z[ig];
                  const double tgx2 = tgx * tgx;
                  const double tgy2 = tgy * tgy;
                  const double tgz2 = tgz * tgz;
 
                  const double tgi = kpgi[ig];
                  const double tgi2 = tgi * tgi;
 
                  const double y1 = s34pi * tgx * tgi;
                  const double y2 = s34pi * tgy * tgi;
                  const double y3 = s34pi * tgz * tgi;
 
                  t1[ig]  = y1 * v;
                  t2[ig]  = y2 * v;
                  t3[ig]  = y3 * v;

                  if (compute_stress)
                  {
                     const double fac1 = - y1 * ( v - tg * dv ) * tgi2;
                     // m=x
                     dt1_xx[ig] = fac1 * tgx2 + v * y1;
                     dt1_yy[ig] = fac1 * tgy2;
                     dt1_zz[ig] = fac1 * tgz2;
                     dt1_xy[ig] = fac1 * tgx * tgy;
                     dt1_yz[ig] = fac1 * tgy * tgz;
                     dt1_xz[ig] = fac1 * tgx * tgz;
 
                     const double fac2 = - y2 * ( v - tg * dv ) * tgi2;
                     // m=y
                     dt2_xx[ig] = fac2 * tgx2;
                     dt2_yy[ig] = fac2 * tgy2 + v * y2;
                     dt2_zz[ig] = fac2 * tgz2;
                     dt2_xy[ig] = fac2 * tgx * tgy + v * y1;
                     dt2_yz[ig] = fac2 * tgy * tgz;
                     dt2_xz[ig] = fac2 * tgx * tgz;
                     
                     const double fac3 = - y3 * ( v - tg * dv ) * tgi2;
                     // m=z
                     dt3_xx[ig] = fac3 * tgx2;
                     dt3_yy[ig] = fac3 * tgy2;
                     dt3_zz[ig] = fac3 * tgz2 + v * y3;
                     dt3_xy[ig] = fac3 * tgx * tgy;
                     dt3_yz[ig] = fac3 * tgy * tgz + v * y2;
                     dt3_xz[ig] = fac3 * tgx * tgz + v * y1;
                  }
                } // ig
              }
            }
            ilm += 2*l+1;
          }
          else if ( l == 2 ) {
            if ( nquad[is] == 0 ) {
	      for(int ic = 0; ic < s->nchannels(); ic++){
		// Kleinman-Bylander
		const int ipr4 = iprojlm[is][l][0][ic];
		const int ipr5 = iprojlm[is][l][1][ic];
		const int ipr6 = iprojlm[is][l][2][ic];
		const int ipr7 = iprojlm[is][l][3][ic];
		const int ipr8 = iprojlm[is][l][4][ic];
            
		double *t4 = &twnl[is][ngwl*ipr4];
		double *t5 = &twnl[is][ngwl*ipr5];
		double *t6 = &twnl[is][ngwl*ipr6];
		double *t7 = &twnl[is][ngwl*ipr7];
		double *t8 = &twnl[is][ngwl*ipr8];
            
		// dtwnl[is][ipr][ij][ngwl]
		// index = ig + ngwl * ( ij + 6 * ipr )
		double *dt4_xx,*dt4_yy,*dt4_zz,*dt4_xy,*dt4_yz,*dt4_xz;
		double *dt5_xx,*dt5_yy,*dt5_zz,*dt5_xy,*dt5_yz,*dt5_xz;
		double *dt6_xx,*dt6_yy,*dt6_zz,*dt6_xy,*dt6_yz,*dt6_xz;
		double *dt7_xx,*dt7_yy,*dt7_zz,*dt7_xy,*dt7_yz,*dt7_xz;
		double *dt8_xx,*dt8_yy,*dt8_zz,*dt8_xy,*dt8_yz,*dt8_xz;
		if (compute_stress)
		  {
		    dt4_xx = &dtwnl[is][ngwl*(0+6*ipr4)];
		    dt4_yy = &dtwnl[is][ngwl*(1+6*ipr4)];
		    dt4_zz = &dtwnl[is][ngwl*(2+6*ipr4)];
		    dt4_xy = &dtwnl[is][ngwl*(3+6*ipr4)];
		    dt4_yz = &dtwnl[is][ngwl*(4+6*ipr4)];
		    dt4_xz = &dtwnl[is][ngwl*(5+6*ipr4)];

		    dt5_xx = &dtwnl[is][ngwl*(0+6*ipr5)];
		    dt5_yy = &dtwnl[is][ngwl*(1+6*ipr5)];
		    dt5_zz = &dtwnl[is][ngwl*(2+6*ipr5)];
		    dt5_xy = &dtwnl[is][ngwl*(3+6*ipr5)];
		    dt5_yz = &dtwnl[is][ngwl*(4+6*ipr5)];
		    dt5_xz = &dtwnl[is][ngwl*(5+6*ipr5)];
            
		    dt6_xx = &dtwnl[is][ngwl*(0+6*ipr6)];
		    dt6_yy = &dtwnl[is][ngwl*(1+6*ipr6)];
		    dt6_zz = &dtwnl[is][ngwl*(2+6*ipr6)];
		    dt6_xy = &dtwnl[is][ngwl*(3+6*ipr6)];
		    dt6_yz = &dtwnl[is][ngwl*(4+6*ipr6)];
		    dt6_xz = &dtwnl[is][ngwl*(5+6*ipr6)];
            
		    dt7_xx = &dtwnl[is][ngwl*(0+6*ipr7)];
		    dt7_yy = &dtwnl[is][ngwl*(1+6*ipr7)];
		    dt7_zz = &dtwnl[is][ngwl*(2+6*ipr7)];
		    dt7_xy = &dtwnl[is][ngwl*(3+6*ipr7)];
		    dt7_yz = &dtwnl[is][ngwl*(4+6*ipr7)];
		    dt7_xz = &dtwnl[is][ngwl*(5+6*ipr7)];
            
		    dt8_xx = &dtwnl[is][ngwl*(0+6*ipr8)];
		    dt8_yy = &dtwnl[is][ngwl*(1+6*ipr8)];
		    dt8_zz = &dtwnl[is][ngwl*(2+6*ipr8)];
		    dt8_xy = &dtwnl[is][ngwl*(3+6*ipr8)];
		    dt8_yz = &dtwnl[is][ngwl*(4+6*ipr8)];
		    dt8_xz = &dtwnl[is][ngwl*(5+6*ipr8)];
		  }
              
		for ( int ig = 0; ig < ngwl; ig++ ) {
		  double v,dv;
		  const double tg = kpg[ig];
              
		  s->dvnlg(l, ic, tg, v, dv);
              
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
              
		  const double y1x = s34pi * tgx * tgi;
		  const double y1y = s34pi * tgy * tgi;
		  const double y1z = s34pi * tgz * tgi;
              
		  const double dx_xx = y1x * tgxx - y1x;
		  const double dx_yy = y1x * tgyy;
		  const double dx_zz = y1x * tgzz;
		  const double dx_xy = y1x * tgxy;
		  const double dx_yz = y1x * tgyz;
		  const double dx_xz = y1x * tgxz;
              
		  const double dy_xx = y1y * tgxx;
		  const double dy_yy = y1y * tgyy - y1y;
		  const double dy_zz = y1y * tgzz;
		  const double dy_xy = y1y * tgxy - y1x;
		  const double dy_yz = y1y * tgyz;
		  const double dy_xz = y1y * tgxz;
              
		  const double dz_xx = y1z * tgxx;
		  const double dz_yy = y1z * tgyy;
		  const double dz_zz = y1z * tgzz - y1z;
		  const double dz_xy = y1z * tgxy;
		  const double dz_yz = y1z * tgyz - y1y;
		  const double dz_xz = y1z * tgxz - y1x;
              
		  t4[ig]  = y4 * v;
		  t5[ig]  = y5 * v;
		  t6[ig]  = y6 * v;
		  t7[ig]  = y7 * v;
		  t8[ig]  = y8 * v;

		  if (compute_stress)
		    {
		      // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
		      dt4_xx[ig] = -(v * s20pi * dz_xx * y1z - y4 * dv * tg * tgxx);
		      dt4_yy[ig] = -(v * s20pi * dz_yy * y1z - y4 * dv * tg * tgyy);
		      dt4_zz[ig] = -(v * s20pi * dz_zz * y1z - y4 * dv * tg * tgzz);
		      dt4_xy[ig] = -(v * s20pi * dz_xy * y1z - y4 * dv * tg * tgxy);
		      dt4_yz[ig] = -(v * s20pi * dz_yz * y1z - y4 * dv * tg * tgyz);
		      dt4_xz[ig] = -(v * s20pi * dz_xz * y1z - y4 * dv * tg * tgxz);
              
		      // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
		      dt5_xx[ig] = -(v * s20pi3 * (y1x * dx_xx - y1y * dy_xx) - y5 * dv * tg * tgxx);
		      dt5_yy[ig] = -(v * s20pi3 * (y1x * dx_yy - y1y * dy_yy) - y5 * dv * tg * tgyy);
		      dt5_zz[ig] = -(v * s20pi3 * (y1x * dx_zz - y1y * dy_zz) - y5 * dv * tg * tgzz);
		      dt5_xy[ig] = -(v * s20pi3 * (y1x * dx_xy - y1y * dy_xy) - y5 * dv * tg * tgxy);
		      dt5_yz[ig] = -(v * s20pi3 * (y1x * dx_yz - y1y * dy_yz) - y5 * dv * tg * tgyz);
		      dt5_xz[ig] = -(v * s20pi3 * (y1x * dx_xz - y1y * dy_xz) - y5 * dv * tg * tgxz);

		      // y6 = s54pi sqrt(3) x y / r^2
		      dt6_xx[ig] = -(v * s20pi3 * (dx_xx * y1y + y1x * dy_xx) - y6 * dv * tg * tgxx);
		      dt6_yy[ig] = -(v * s20pi3 * (dx_yy * y1y + y1x * dy_yy) - y6 * dv * tg * tgyy);
		      dt6_zz[ig] = -(v * s20pi3 * (dx_zz * y1y + y1x * dy_zz) - y6 * dv * tg * tgzz);
		      dt6_xy[ig] = -(v * s20pi3 * (dx_xy * y1y + y1x * dy_xy) - y6 * dv * tg * tgxy);
		      dt6_yz[ig] = -(v * s20pi3 * (dx_yz * y1y + y1x * dy_yz) - y6 * dv * tg * tgyz);
		      dt6_xz[ig] = -(v * s20pi3 * (dx_xz * y1y + y1x * dy_xz) - y6 * dv * tg * tgxz);

		      // y7 = s54pi sqrt(3) y z / r^2
		      dt7_xx[ig] = -(v * s20pi3 * (dy_xx * y1z + y1y * dz_xx) - y7 * dv * tg * tgxx);
		      dt7_yy[ig] = -(v * s20pi3 * (dy_yy * y1z + y1y * dz_yy) - y7 * dv * tg * tgyy);
		      dt7_zz[ig] = -(v * s20pi3 * (dy_zz * y1z + y1y * dz_zz) - y7 * dv * tg * tgzz);
		      dt7_xy[ig] = -(v * s20pi3 * (dy_xy * y1z + y1y * dz_xy) - y7 * dv * tg * tgxy);
		      dt7_yz[ig] = -(v * s20pi3 * (dy_yz * y1z + y1y * dz_yz) - y7 * dv * tg * tgyz);
		      dt7_xz[ig] = -(v * s20pi3 * (dy_xz * y1z + y1y * dz_xz) - y7 * dv * tg * tgxz);
              
		      // y8 = s54pi sqrt(3) z x / r^2
		      dt8_xx[ig] = -(v * s20pi3 * (dx_xx * y1z + y1x * dz_xx) - y8 * dv * tg * tgxx);
		      dt8_yy[ig] = -(v * s20pi3 * (dx_yy * y1z + y1x * dz_yy) - y8 * dv * tg * tgyy);
		      dt8_zz[ig] = -(v * s20pi3 * (dx_zz * y1z + y1x * dz_zz) - y8 * dv * tg * tgzz);
		      dt8_xy[ig] = -(v * s20pi3 * (dx_xy * y1z + y1x * dz_xy) - y8 * dv * tg * tgxy);
		      dt8_yz[ig] = -(v * s20pi3 * (dx_yz * y1z + y1x * dz_yz) - y8 * dv * tg * tgyz);
		      dt8_xz[ig] = -(v * s20pi3 * (dx_xz * y1z + y1x * dz_xz) - y8 * dv * tg * tgxz);
		    }
		}
              }
            }
            else {
              // semi-local
              for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                // twnl[is][ipr][ig]
                // ipr = iquad+nquad[is]*ilm
                // index = ig + ngwl*ipr
                const int ipr4 = iquad+nquad[is]*ilm;
                const int ipr5 = iquad+nquad[is]*(ilm+1);
                const int ipr6 = iquad+nquad[is]*(ilm+2);
                const int ipr7 = iquad+nquad[is]*(ilm+3);
                const int ipr8 = iquad+nquad[is]*(ilm+4);
            
                double *t4 = &twnl[is][ngwl*ipr4];
                double *t5 = &twnl[is][ngwl*ipr5];
                double *t6 = &twnl[is][ngwl*ipr6];
                double *t7 = &twnl[is][ngwl*ipr7];
                double *t8 = &twnl[is][ngwl*ipr8];

                // dtwnl[is][ipr][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * ipr )
                double *dt4_xx,*dt4_yy,*dt4_zz,*dt4_xy,*dt4_yz,*dt4_xz;
                double *dt5_xx,*dt5_yy,*dt5_zz,*dt5_xy,*dt5_yz,*dt5_xz;
                double *dt6_xx,*dt6_yy,*dt6_zz,*dt6_xy,*dt6_yz,*dt6_xz;
                double *dt7_xx,*dt7_yy,*dt7_zz,*dt7_xy,*dt7_yz,*dt7_xz;
                double *dt8_xx,*dt8_yy,*dt8_zz,*dt8_xy,*dt8_yz,*dt8_xz;
                if (compute_stress)
                {
                   dt4_xx = &dtwnl[is][ngwl*(0+6*ipr4)];
                   dt4_yy = &dtwnl[is][ngwl*(1+6*ipr4)];
                   dt4_zz = &dtwnl[is][ngwl*(2+6*ipr4)];
                   dt4_xy = &dtwnl[is][ngwl*(3+6*ipr4)];
                   dt4_yz = &dtwnl[is][ngwl*(4+6*ipr4)];
                   dt4_xz = &dtwnl[is][ngwl*(5+6*ipr4)];

                   dt5_xx = &dtwnl[is][ngwl*(0+6*ipr5)];
                   dt5_yy = &dtwnl[is][ngwl*(1+6*ipr5)];
                   dt5_zz = &dtwnl[is][ngwl*(2+6*ipr5)];
                   dt5_xy = &dtwnl[is][ngwl*(3+6*ipr5)];
                   dt5_yz = &dtwnl[is][ngwl*(4+6*ipr5)];
                   dt5_xz = &dtwnl[is][ngwl*(5+6*ipr5)];
 
                   dt6_xx = &dtwnl[is][ngwl*(0+6*ipr6)];
                   dt6_yy = &dtwnl[is][ngwl*(1+6*ipr6)];
                   dt6_zz = &dtwnl[is][ngwl*(2+6*ipr6)];
                   dt6_xy = &dtwnl[is][ngwl*(3+6*ipr6)];
                   dt6_yz = &dtwnl[is][ngwl*(4+6*ipr6)];
                   dt6_xz = &dtwnl[is][ngwl*(5+6*ipr6)];
 
                   dt7_xx = &dtwnl[is][ngwl*(0+6*ipr7)];
                   dt7_yy = &dtwnl[is][ngwl*(1+6*ipr7)];
                   dt7_zz = &dtwnl[is][ngwl*(2+6*ipr7)];
                   dt7_xy = &dtwnl[is][ngwl*(3+6*ipr7)];
                   dt7_yz = &dtwnl[is][ngwl*(4+6*ipr7)];
                   dt7_xz = &dtwnl[is][ngwl*(5+6*ipr7)];
 
                   dt8_xx = &dtwnl[is][ngwl*(0+6*ipr8)];
                   dt8_yy = &dtwnl[is][ngwl*(1+6*ipr8)];
                   dt8_zz = &dtwnl[is][ngwl*(2+6*ipr8)];
                   dt8_xy = &dtwnl[is][ngwl*(3+6*ipr8)];
                   dt8_yz = &dtwnl[is][ngwl*(4+6*ipr8)];
                   dt8_xz = &dtwnl[is][ngwl*(5+6*ipr8)];
                }
                
                const double r = rquad[is][iquad];
                for ( int ig = 0; ig < ngwl; ig++ ) {
                  double v = 0.0, dv = 0.0;
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
                    const double dj2 = 
                        ( (4.0 * z2 - 9.0) * s + z*(9.0-z2) * c ) / (z2*z2) ;
                    // v = 4 pi j2(Gr) r
                    v = fpi * j2 * r;
                    // dv = d/dG v = 4 pi dj2(Gr)/d(Gr) d(Gr)/dG r
                    //    = 4 pi dj2 r^2
                    dv = fpi * dj2 * r * r;                  
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
 
                  const double y1x = s34pi * tgx * tgi;
                  const double y1y = s34pi * tgy * tgi;
                  const double y1z = s34pi * tgz * tgi;
 
                  const double dx_xx = y1x * tgxx - y1x;
                  const double dx_yy = y1x * tgyy;
                  const double dx_zz = y1x * tgzz;
                  const double dx_xy = y1x * tgxy;
                  const double dx_yz = y1x * tgyz;
                  const double dx_xz = y1x * tgxz;
 
                  const double dy_xx = y1y * tgxx;
                  const double dy_yy = y1y * tgyy - y1y;
                  const double dy_zz = y1y * tgzz;
                  const double dy_xy = y1y * tgxy - y1x;
                  const double dy_yz = y1y * tgyz;
                  const double dy_xz = y1y * tgxz;
 
                  const double dz_xx = y1z * tgxx;
                  const double dz_yy = y1z * tgyy;
                  const double dz_zz = y1z * tgzz - y1z;
                  const double dz_xy = y1z * tgxy;
                  const double dz_yz = y1z * tgyz - y1y;
                  const double dz_xz = y1z * tgxz - y1x;
 
                  t4[ig]  = y4 * v;
                  t5[ig]  = y5 * v;
                  t6[ig]  = y6 * v;
                  t7[ig]  = y7 * v;
                  t8[ig]  = y8 * v;

                  if (compute_stress)
                  {
                     // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
                     dt4_xx[ig] = -(v * s20pi * dz_xx * y1z - y4 * dv * tg * tgxx);
                     dt4_yy[ig] = -(v * s20pi * dz_yy * y1z - y4 * dv * tg * tgyy);
                     dt4_zz[ig] = -(v * s20pi * dz_zz * y1z - y4 * dv * tg * tgzz);
                     dt4_xy[ig] = -(v * s20pi * dz_xy * y1z - y4 * dv * tg * tgxy);
                     dt4_yz[ig] = -(v * s20pi * dz_yz * y1z - y4 * dv * tg * tgyz);
                     dt4_xz[ig] = -(v * s20pi * dz_xz * y1z - y4 * dv * tg * tgxz);
 
                     // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
                     dt5_xx[ig] = -(v * s20pi3 * (y1x * dx_xx - y1y * dy_xx) - y5 * dv * tg * tgxx);
                     dt5_yy[ig] = -(v * s20pi3 * (y1x * dx_yy - y1y * dy_yy) - y5 * dv * tg * tgyy);
                     dt5_zz[ig] = -(v * s20pi3 * (y1x * dx_zz - y1y * dy_zz) - y5 * dv * tg * tgzz);
                     dt5_xy[ig] = -(v * s20pi3 * (y1x * dx_xy - y1y * dy_xy) - y5 * dv * tg * tgxy);
                     dt5_yz[ig] = -(v * s20pi3 * (y1x * dx_yz - y1y * dy_yz) - y5 * dv * tg * tgyz);
                     dt5_xz[ig] = -(v * s20pi3 * (y1x * dx_xz - y1y * dy_xz) - y5 * dv * tg * tgxz);

                     // y6 = s54pi sqrt(3) x y / r^2
                     dt6_xx[ig] = -(v * s20pi3 * (dx_xx * y1y + y1x * dy_xx) - y6 * dv * tg * tgxx);
                     dt6_yy[ig] = -(v * s20pi3 * (dx_yy * y1y + y1x * dy_yy) - y6 * dv * tg * tgyy);
                     dt6_zz[ig] = -(v * s20pi3 * (dx_zz * y1y + y1x * dy_zz) - y6 * dv * tg * tgzz);
                     dt6_xy[ig] = -(v * s20pi3 * (dx_xy * y1y + y1x * dy_xy) - y6 * dv * tg * tgxy);
                     dt6_yz[ig] = -(v * s20pi3 * (dx_yz * y1y + y1x * dy_yz) - y6 * dv * tg * tgyz);
                     dt6_xz[ig] = -(v * s20pi3 * (dx_xz * y1y + y1x * dy_xz) - y6 * dv * tg * tgxz);

                     // y7 = s54pi sqrt(3) y z / r^2
                     dt7_xx[ig] = -(v * s20pi3 * (dy_xx * y1z + y1y * dz_xx) - y7 * dv * tg * tgxx);
                     dt7_yy[ig] = -(v * s20pi3 * (dy_yy * y1z + y1y * dz_yy) - y7 * dv * tg * tgyy);
                     dt7_zz[ig] = -(v * s20pi3 * (dy_zz * y1z + y1y * dz_zz) - y7 * dv * tg * tgzz);
                     dt7_xy[ig] = -(v * s20pi3 * (dy_xy * y1z + y1y * dz_xy) - y7 * dv * tg * tgxy);
                     dt7_yz[ig] = -(v * s20pi3 * (dy_yz * y1z + y1y * dz_yz) - y7 * dv * tg * tgyz);
                     dt7_xz[ig] = -(v * s20pi3 * (dy_xz * y1z + y1y * dz_xz) - y7 * dv * tg * tgxz);
 
                     // y8 = s54pi sqrt(3) z x / r^2
                     dt8_xx[ig] = -(v * s20pi3 * (dx_xx * y1z + y1x * dz_xx) - y8 * dv * tg * tgxx);
                     dt8_yy[ig] = -(v * s20pi3 * (dx_yy * y1z + y1x * dz_yy) - y8 * dv * tg * tgyy);
                     dt8_zz[ig] = -(v * s20pi3 * (dx_zz * y1z + y1x * dz_zz) - y8 * dv * tg * tgzz);
                     dt8_xy[ig] = -(v * s20pi3 * (dx_xy * y1z + y1x * dz_xy) - y8 * dv * tg * tgxy);
                     dt8_yz[ig] = -(v * s20pi3 * (dx_yz * y1z + y1x * dz_yz) - y8 * dv * tg * tgyz);
                     dt8_xz[ig] = -(v * s20pi3 * (dx_xz * y1z + y1x * dz_xz) - y8 * dv * tg * tgxz);
                  }
                } // ig
              } // iquad
            }
            ilm += 2*l+1;
          }
          else if ( l == 3 ) {
            if ( nquad[is] == 0 ) {
	      for(int ic = 0; ic < s->nchannels(); ic++){
              // Kleinman-Bylander
		const int ipr9  = iprojlm[is][l][0][ic];
		const int ipr10 = iprojlm[is][l][1][ic];
		const int ipr11 = iprojlm[is][l][2][ic];
		const int ipr12 = iprojlm[is][l][3][ic];
		const int ipr13 = iprojlm[is][l][4][ic];
		const int ipr14 = iprojlm[is][l][5][ic];
		const int ipr15 = iprojlm[is][l][6][ic];
            
		double *t9 = &twnl[is][ngwl*ipr9];
		double *t10 = &twnl[is][ngwl*ipr10];
		double *t11 = &twnl[is][ngwl*ipr11];
		double *t12 = &twnl[is][ngwl*ipr12];
		double *t13 = &twnl[is][ngwl*ipr13];
		double *t14 = &twnl[is][ngwl*ipr14];
		double *t15 = &twnl[is][ngwl*ipr15];
            
		// dtwnl[is][ipr][ij][ngwl]
		// index = ig + ngwl * ( ij + 6 * ipr )
		double *dt9_xx,*dt9_yy,*dt9_zz,*dt9_xy,*dt9_yz,*dt9_xz;
		double *dt10_xx,*dt10_yy,*dt10_zz,*dt10_xy,*dt10_yz,*dt10_xz;
		double *dt11_xx,*dt11_yy,*dt11_zz,*dt11_xy,*dt11_yz,*dt11_xz;
		double *dt12_xx,*dt12_yy,*dt12_zz,*dt12_xy,*dt12_yz,*dt12_xz;
		double *dt13_xx,*dt13_yy,*dt13_zz,*dt13_xy,*dt13_yz,*dt13_xz;
		double *dt14_xx,*dt14_yy,*dt14_zz,*dt14_xy,*dt14_yz,*dt14_xz;
		double *dt15_xx,*dt15_yy,*dt15_zz,*dt15_xy,*dt15_yz,*dt15_xz;
		if (compute_stress)
		  {
		    dt9_xx = &dtwnl[is][ngwl*(0+6*ipr9)];
		    dt9_yy = &dtwnl[is][ngwl*(1+6*ipr9)];
		    dt9_zz = &dtwnl[is][ngwl*(2+6*ipr9)];
		    dt9_xy = &dtwnl[is][ngwl*(3+6*ipr9)];
		    dt9_yz = &dtwnl[is][ngwl*(4+6*ipr9)];
		    dt9_xz = &dtwnl[is][ngwl*(5+6*ipr9)];

		    dt10_xx = &dtwnl[is][ngwl*(0+6*ipr10)];
		    dt10_yy = &dtwnl[is][ngwl*(1+6*ipr10)];
		    dt10_zz = &dtwnl[is][ngwl*(2+6*ipr10)];
		    dt10_xy = &dtwnl[is][ngwl*(3+6*ipr10)];
		    dt10_yz = &dtwnl[is][ngwl*(4+6*ipr10)];
		    dt10_xz = &dtwnl[is][ngwl*(5+6*ipr10)];
            
		    dt11_xx = &dtwnl[is][ngwl*(0+6*ipr11)];
		    dt11_yy = &dtwnl[is][ngwl*(1+6*ipr11)];
		    dt11_zz = &dtwnl[is][ngwl*(2+6*ipr11)];
		    dt11_xy = &dtwnl[is][ngwl*(3+6*ipr11)];
		    dt11_yz = &dtwnl[is][ngwl*(4+6*ipr11)];
		    dt11_xz = &dtwnl[is][ngwl*(5+6*ipr11)];
            
		    dt12_xx = &dtwnl[is][ngwl*(0+6*ipr12)];
		    dt12_yy = &dtwnl[is][ngwl*(1+6*ipr12)];
		    dt12_zz = &dtwnl[is][ngwl*(2+6*ipr12)];
		    dt12_xy = &dtwnl[is][ngwl*(3+6*ipr12)];
		    dt12_yz = &dtwnl[is][ngwl*(4+6*ipr12)];
		    dt12_xz = &dtwnl[is][ngwl*(5+6*ipr12)];
            
		    dt13_xx = &dtwnl[is][ngwl*(0+6*ipr13)];
		    dt13_yy = &dtwnl[is][ngwl*(1+6*ipr13)];
		    dt13_zz = &dtwnl[is][ngwl*(2+6*ipr13)];
		    dt13_xy = &dtwnl[is][ngwl*(3+6*ipr13)];
		    dt13_yz = &dtwnl[is][ngwl*(4+6*ipr13)];
		    dt13_xz = &dtwnl[is][ngwl*(5+6*ipr13)];

		    dt14_xx = &dtwnl[is][ngwl*(0+6*ipr14)];
		    dt14_yy = &dtwnl[is][ngwl*(1+6*ipr14)];
		    dt14_zz = &dtwnl[is][ngwl*(2+6*ipr14)];
		    dt14_xy = &dtwnl[is][ngwl*(3+6*ipr14)];
		    dt14_yz = &dtwnl[is][ngwl*(4+6*ipr14)];
		    dt14_xz = &dtwnl[is][ngwl*(5+6*ipr14)];
            
		    dt15_xx = &dtwnl[is][ngwl*(0+6*ipr15)];
		    dt15_yy = &dtwnl[is][ngwl*(1+6*ipr15)];
		    dt15_zz = &dtwnl[is][ngwl*(2+6*ipr15)];
		    dt15_xy = &dtwnl[is][ngwl*(3+6*ipr15)];
		    dt15_yz = &dtwnl[is][ngwl*(4+6*ipr15)];
		    dt15_xz = &dtwnl[is][ngwl*(5+6*ipr15)];
		  }
              
		for ( int ig = 0; ig < ngwl; ig++ ) {
		  double v,dv;
		  const double tg = kpg[ig];
              
		  s->dvnlg(l, ic, tg, v, dv);
              
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
              
		  const double y1x = tgx * tgi;
		  const double y1y = tgy * tgi;
		  const double y1z = tgz * tgi;
              
		  const double dx_xx = y1x * tgxx - y1x;
		  const double dx_yy = y1x * tgyy;
		  const double dx_zz = y1x * tgzz;
		  const double dx_xy = y1x * tgxy;
		  const double dx_yz = y1x * tgyz;
		  const double dx_xz = y1x * tgxz;
              
		  const double dy_xx = y1y * tgxx;
		  const double dy_yy = y1y * tgyy - y1y;
		  const double dy_zz = y1y * tgzz;
		  const double dy_xy = y1y * tgxy - y1x;
		  const double dy_yz = y1y * tgyz;
		  const double dy_xz = y1y * tgxz;
              
		  const double dz_xx = y1z * tgxx;
		  const double dz_yy = y1z * tgyy;
		  const double dz_zz = y1z * tgzz - y1z;
		  const double dz_xy = y1z * tgxy;
		  const double dz_yz = y1z * tgyz - y1y;
		  const double dz_xz = y1z * tgxz - y1x;
              
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

		  //ewd: notes on l=3 stress derivation
		  //
		  // stress derivatives use change of variables:  x = tgx*tgi, y = tgy*tgi, z = tgz*tgi
		  //
		  // the nonlocal stress contribution is given by (e.g.):
		  //   dy9/deps_ab = dy9/dx dx/deps_ab + dy9/dy dy/deps_ab + dy9/dz dz/deps_ab
		  //
		  // where (e.g.) 
		  //   dx/deps_ab = dx/dtgx dtgx/deps_ab + dx/dtgy dtgy/deps_ab + dx/dtgz dtgz/deps_ab
		  // and
		  //   dtgx/deps_ab = - delta_xa * tgb
		  //
		  // note:  we follow Francois' convention of pulling the minus sign out of 
		  //   the dx_xx,dz_yy, etc. terms

		  if (compute_stress)
		    {
		      // y9 = s74pi * 0.5 * (5*z^3 - 3z)
		      const double dy9z = s74pi * 0.5 * (15.*tgzz - 3.);
		      dt9_xx[ig] = -dy9z * dz_xx * v + y9 * tg * tgxx * dv;
		      dt9_yy[ig] = -dy9z * dz_yy * v + y9 * tg * tgyy * dv;
		      dt9_zz[ig] = -dy9z * dz_zz * v + y9 * tg * tgzz * dv;
		      dt9_xy[ig] = -dy9z * dz_xy * v + y9 * tg * tgxy * dv;
		      dt9_yz[ig] = -dy9z * dz_yz * v + y9 * tg * tgyz * dv;
		      dt9_xz[ig] = -dy9z * dz_xz * v + y9 * tg * tgxz * dv;

		      // y10 = s2132pi * x * (5z^2 - 1)
		      const double dy10x = s2132pi * (5.*tgzz-1.0);
		      const double dy10z = s2132pi * 10.*tgxz;
		      dt10_xx[ig] = -(dy10x*dx_xx + dy10z*dz_xx) * v + y10 * tg * tgxx * dv;
		      dt10_yy[ig] = -(dy10x*dx_yy + dy10z*dz_yy) * v + y10 * tg * tgyy * dv;
		      dt10_zz[ig] = -(dy10x*dx_zz + dy10z*dz_zz) * v + y10 * tg * tgzz * dv;
		      dt10_xy[ig] = -(dy10x*dx_xy + dy10z*dz_xy) * v + y10 * tg * tgxy * dv;
		      dt10_yz[ig] = -(dy10x*dx_yz + dy10z*dz_yz) * v + y10 * tg * tgyz * dv;
		      dt10_xz[ig] = -(dy10x*dx_xz + dy10z*dz_xz) * v + y10 * tg * tgxz * dv;


		      // y11 = s2132pi * y * (5z^2 - 1)
		      const double dy11y = s2132pi * (5.*tgzz-1.0);
		      const double dy11z = s2132pi * 10.*tgyz;
		      dt11_xx[ig] = -(dy11y*dy_xx + dy11z*dz_xx) * v + y11 * tg * tgxx * dv;
		      dt11_yy[ig] = -(dy11y*dy_yy + dy11z*dz_yy) * v + y11 * tg * tgyy * dv;
		      dt11_zz[ig] = -(dy11y*dy_zz + dy11z*dz_zz) * v + y11 * tg * tgzz * dv;
		      dt11_xy[ig] = -(dy11y*dy_xy + dy11z*dz_xy) * v + y11 * tg * tgxy * dv;
		      dt11_yz[ig] = -(dy11y*dy_yz + dy11z*dz_yz) * v + y11 * tg * tgyz * dv;
		      dt11_xz[ig] = -(dy11y*dy_xz + dy11z*dz_xz) * v + y11 * tg * tgxz * dv;

		      // y12 = s1054pi * x * y * z;
		      const double dy12x = s1054pi * tgyz;
		      const double dy12y = s1054pi * tgxz;
		      const double dy12z = s1054pi * tgxy;
		      dt12_xx[ig] = -(dy12x*dx_xx + dy12y*dy_xx + dy12z*dz_xx) * v + y12 * tg * tgxx * dv;
		      dt12_yy[ig] = -(dy12x*dx_yy + dy12y*dy_yy + dy12z*dz_yy) * v + y12 * tg * tgyy * dv;
		      dt12_zz[ig] = -(dy12x*dx_zz + dy12y*dy_zz + dy12z*dz_zz) * v + y12 * tg * tgzz * dv;
		      dt12_xy[ig] = -(dy12x*dx_xy + dy12y*dy_xy + dy12z*dz_xy) * v + y12 * tg * tgxy * dv;
		      dt12_yz[ig] = -(dy12x*dx_yz + dy12y*dy_yz + dy12z*dz_yz) * v + y12 * tg * tgyz * dv;
		      dt12_xz[ig] = -(dy12x*dx_xz + dy12y*dy_xz + dy12z*dz_xz) * v + y12 * tg * tgxz * dv;

		      // y13 = s1054pi * 0.5 * z * (x^2 - y^2)
		      const double dy13x = s1054pi * tgxz;
		      const double dy13y = -s1054pi * tgyz;
		      const double dy13z = s1054pi * 0.5 * (tgxx - tgyy);
		      dt13_xx[ig] = -(dy13x*dx_xx + dy13y*dy_xx + dy13z*dz_xx) * v + y13 * tg * tgxx * dv;
		      dt13_yy[ig] = -(dy13x*dx_yy + dy13y*dy_yy + dy13z*dz_yy) * v + y13 * tg * tgyy * dv;
		      dt13_zz[ig] = -(dy13x*dx_zz + dy13y*dy_zz + dy13z*dz_zz) * v + y13 * tg * tgzz * dv;
		      dt13_xy[ig] = -(dy13x*dx_xy + dy13y*dy_xy + dy13z*dz_xy) * v + y13 * tg * tgxy * dv;
		      dt13_yz[ig] = -(dy13x*dx_yz + dy13y*dy_yz + dy13z*dz_yz) * v + y13 * tg * tgyz * dv;
		      dt13_xz[ig] = -(dy13x*dx_xz + dy13y*dy_xz + dy13z*dz_xz) * v + y13 * tg * tgxz * dv;

		      // y14 = s3532pi * (x^3 - 3*x*y^2)
		      const double dy14x = s3532pi * 3. * (tgxx - tgyy);
		      const double dy14y = -s3532pi * 6. * tgxy;
		      dt14_xx[ig] = -(dy14x*dx_xx + dy14y*dy_xx) * v + y14 * tg * tgxx * dv;
		      dt14_yy[ig] = -(dy14x*dx_yy + dy14y*dy_yy) * v + y14 * tg * tgyy * dv;
		      dt14_zz[ig] = -(dy14x*dx_zz + dy14y*dy_zz) * v + y14 * tg * tgzz * dv;
		      dt14_xy[ig] = -(dy14x*dx_xy + dy14y*dy_xy) * v + y14 * tg * tgxy * dv;
		      dt14_yz[ig] = -(dy14x*dx_yz + dy14y*dy_yz) * v + y14 * tg * tgyz * dv;
		      dt14_xz[ig] = -(dy14x*dx_xz + dy14y*dy_xz) * v + y14 * tg * tgxz * dv;

		      // y15 = s3532pi * (3*y*x^2 - y^3)
		      const double dy15x = s3532pi * 6.* tgxy;
		      const double dy15y = s3532pi * 3. * (tgxx - tgyy);
		      dt15_xx[ig] = -(dy15x*dx_xx + dy15y*dy_xx) * v + y15 * tg * tgxx * dv;
		      dt15_yy[ig] = -(dy15x*dx_yy + dy15y*dy_yy) * v + y15 * tg * tgyy * dv;
		      dt15_zz[ig] = -(dy15x*dx_zz + dy15y*dy_zz) * v + y15 * tg * tgzz * dv;
		      dt15_xy[ig] = -(dy15x*dx_xy + dy15y*dy_xy) * v + y15 * tg * tgxy * dv;
		      dt15_yz[ig] = -(dy15x*dx_yz + dy15y*dy_yz) * v + y15 * tg * tgyz * dv;
		      dt15_xz[ig] = -(dy15x*dx_xz + dy15y*dy_xz) * v + y15 * tg * tgxz * dv;
		    }
		}
              }
            }
            else {
              // semi-local
              for ( int iquad = 0; iquad < nquad[is]; iquad++ ) {
                // twnl[is][ipr][ig]
                // ipr = iquad+nquad[is]*ilm
                // index = ig + ngwl*ipr
                const int ipr9 = iquad+nquad[is]*ilm;
                const int ipr10 = iquad+nquad[is]*(ilm+1);
                const int ipr11 = iquad+nquad[is]*(ilm+2);
                const int ipr12 = iquad+nquad[is]*(ilm+3);
                const int ipr13 = iquad+nquad[is]*(ilm+4);
                const int ipr14 = iquad+nquad[is]*(ilm+5);
                const int ipr15 = iquad+nquad[is]*(ilm+6);
            
                double *t9 = &twnl[is][ngwl*ipr9];
                double *t10 = &twnl[is][ngwl*ipr10];
                double *t11 = &twnl[is][ngwl*ipr11];
                double *t12 = &twnl[is][ngwl*ipr12];
                double *t13 = &twnl[is][ngwl*ipr13];
                double *t14 = &twnl[is][ngwl*ipr14];
                double *t15 = &twnl[is][ngwl*ipr15];

                // dtwnl[is][ipr][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * ipr )
                double *dt9_xx,*dt9_yy,*dt9_zz,*dt9_xy,*dt9_yz,*dt9_xz;
                double *dt10_xx,*dt10_yy,*dt10_zz,*dt10_xy,*dt10_yz,*dt10_xz;
                double *dt11_xx,*dt11_yy,*dt11_zz,*dt11_xy,*dt11_yz,*dt11_xz;
                double *dt12_xx,*dt12_yy,*dt12_zz,*dt12_xy,*dt12_yz,*dt12_xz;
                double *dt13_xx,*dt13_yy,*dt13_zz,*dt13_xy,*dt13_yz,*dt13_xz;
                double *dt14_xx,*dt14_yy,*dt14_zz,*dt14_xy,*dt14_yz,*dt14_xz;
                double *dt15_xx,*dt15_yy,*dt15_zz,*dt15_xy,*dt15_yz,*dt15_xz;

                if (compute_stress)
                {
                   dt9_xx = &dtwnl[is][ngwl*(0+6*ipr9)];
                   dt9_yy = &dtwnl[is][ngwl*(1+6*ipr9)];
                   dt9_zz = &dtwnl[is][ngwl*(2+6*ipr9)];
                   dt9_xy = &dtwnl[is][ngwl*(3+6*ipr9)];
                   dt9_yz = &dtwnl[is][ngwl*(4+6*ipr9)];
                   dt9_xz = &dtwnl[is][ngwl*(5+6*ipr9)];

                   dt10_xx = &dtwnl[is][ngwl*(0+6*ipr10)];
                   dt10_yy = &dtwnl[is][ngwl*(1+6*ipr10)];
                   dt10_zz = &dtwnl[is][ngwl*(2+6*ipr10)];
                   dt10_xy = &dtwnl[is][ngwl*(3+6*ipr10)];
                   dt10_yz = &dtwnl[is][ngwl*(4+6*ipr10)];
                   dt10_xz = &dtwnl[is][ngwl*(5+6*ipr10)];
 
                   dt11_xx = &dtwnl[is][ngwl*(0+6*ipr11)];
                   dt11_yy = &dtwnl[is][ngwl*(1+6*ipr11)];
                   dt11_zz = &dtwnl[is][ngwl*(2+6*ipr11)];
                   dt11_xy = &dtwnl[is][ngwl*(3+6*ipr11)];
                   dt11_yz = &dtwnl[is][ngwl*(4+6*ipr11)];
                   dt11_xz = &dtwnl[is][ngwl*(5+6*ipr11)];
 
                   dt12_xx = &dtwnl[is][ngwl*(0+6*ipr12)];
                   dt12_yy = &dtwnl[is][ngwl*(1+6*ipr12)];
                   dt12_zz = &dtwnl[is][ngwl*(2+6*ipr12)];
                   dt12_xy = &dtwnl[is][ngwl*(3+6*ipr12)];
                   dt12_yz = &dtwnl[is][ngwl*(4+6*ipr12)];
                   dt12_xz = &dtwnl[is][ngwl*(5+6*ipr12)];
 
                   dt13_xx = &dtwnl[is][ngwl*(0+6*ipr13)];
                   dt13_yy = &dtwnl[is][ngwl*(1+6*ipr13)];
                   dt13_zz = &dtwnl[is][ngwl*(2+6*ipr13)];
                   dt13_xy = &dtwnl[is][ngwl*(3+6*ipr13)];
                   dt13_yz = &dtwnl[is][ngwl*(4+6*ipr13)];
                   dt13_xz = &dtwnl[is][ngwl*(5+6*ipr13)];

                   dt14_xx = &dtwnl[is][ngwl*(0+6*ipr14)];
                   dt14_yy = &dtwnl[is][ngwl*(1+6*ipr14)];
                   dt14_zz = &dtwnl[is][ngwl*(2+6*ipr14)];
                   dt14_xy = &dtwnl[is][ngwl*(3+6*ipr14)];
                   dt14_yz = &dtwnl[is][ngwl*(4+6*ipr14)];
                   dt14_xz = &dtwnl[is][ngwl*(5+6*ipr14)];

                   dt15_xx = &dtwnl[is][ngwl*(0+6*ipr15)];
                   dt15_yy = &dtwnl[is][ngwl*(1+6*ipr15)];
                   dt15_zz = &dtwnl[is][ngwl*(2+6*ipr15)];
                   dt15_xy = &dtwnl[is][ngwl*(3+6*ipr15)];
                   dt15_yz = &dtwnl[is][ngwl*(4+6*ipr15)];
                   dt15_xz = &dtwnl[is][ngwl*(5+6*ipr15)];
                }
                
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
                    const double z2 = z * z;
                    const double dj3 = 
                        ( (-60.*zi + 27.*z - z*z2) * s + (60.0 - 7.*z2) * c ) / (z2*z2) ;
                    // v = 4 pi j3(Gr) r
                    v = fpi * j3 * r;
                    // dv = d/dG v = 4 pi dj3(Gr)/d(Gr) d(Gr)/dG r
                    //    = 4 pi dj3 r^2
                    dv = fpi * dj3 * r * r;                  
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
              
                  const double y1x = tgx * tgi;
                  const double y1y = tgy * tgi;
                  const double y1z = tgz * tgi;
                
                  const double dx_xx = y1x * tgxx - y1x;
                  const double dx_yy = y1x * tgyy;
                  const double dx_zz = y1x * tgzz;
                  const double dx_xy = y1x * tgxy;
                  const double dx_yz = y1x * tgyz;
                  const double dx_xz = y1x * tgxz;
              
                  const double dy_xx = y1y * tgxx;
                  const double dy_yy = y1y * tgyy - y1y;
                  const double dy_zz = y1y * tgzz;
                  const double dy_xy = y1y * tgxy - y1x;
                  const double dy_yz = y1y * tgyz;
                  const double dy_xz = y1y * tgxz;
              
                  const double dz_xx = y1z * tgxx;
                  const double dz_yy = y1z * tgyy;
                  const double dz_zz = y1z * tgzz - y1z;
                  const double dz_xy = y1z * tgxy;
                  const double dz_yz = y1z * tgyz - y1y;
                  const double dz_xz = y1z * tgxz - y1x;
              
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

                  //ewd: notes on l=3 stress derivation
                  //
                  // stress derivatives use change of variables:  x = tgx*tgi, y = tgy*tgi, z = tgz*tgi
                  //
                  // the nonlocal stress contribution is given by (e.g.):
                  //   dy9/deps_ab = dy9/dx dx/deps_ab + dy9/dy dy/deps_ab + dy9/dz dz/deps_ab
                  //
                  // where (e.g.) 
                  //   dx/deps_ab = dx/dtgx dtgx/deps_ab + dx/dtgy dtgy/deps_ab + dx/dtgz dtgz/deps_ab
                  // and
                  //   dtgx/deps_ab = - delta_xa * tgb
                  //
                  // note:  we follow Francois' convention of pulling the minus sign out of 
                  //   the dx_xx,dz_yy, etc. terms
                
                  if (compute_stress)
                  {
                     // y9 = s74pi * 0.5 * (5*z^3 - 3z)
                     const double dy9z = s74pi * 0.5 * (15.*tgzz - 3.);
                     dt9_xx[ig] = -dy9z * dz_xx * v + y9 * tg * tgxx * dv;
                     dt9_yy[ig] = -dy9z * dz_yy * v + y9 * tg * tgyy * dv;
                     dt9_zz[ig] = -dy9z * dz_zz * v + y9 * tg * tgzz * dv;
                     dt9_xy[ig] = -dy9z * dz_xy * v + y9 * tg * tgxy * dv;
                     dt9_yz[ig] = -dy9z * dz_yz * v + y9 * tg * tgyz * dv;
                     dt9_xz[ig] = -dy9z * dz_xz * v + y9 * tg * tgxz * dv;
                
                     // y10 = s2132pi * x * (5z^2 - 1)
                     const double dy10x = s2132pi * (5.*tgzz-1.0);
                     const double dy10z = s2132pi * 10.*tgxz;
                     dt10_xx[ig] = -(dy10x*dx_xx + dy10z*dz_xx) * v + y10 * tg * tgxx * dv;
                     dt10_yy[ig] = -(dy10x*dx_yy + dy10z*dz_yy) * v + y10 * tg * tgyy * dv;
                     dt10_zz[ig] = -(dy10x*dx_zz + dy10z*dz_zz) * v + y10 * tg * tgzz * dv;
                     dt10_xy[ig] = -(dy10x*dx_xy + dy10z*dz_xy) * v + y10 * tg * tgxy * dv;
                     dt10_yz[ig] = -(dy10x*dx_yz + dy10z*dz_yz) * v + y10 * tg * tgyz * dv;
                     dt10_xz[ig] = -(dy10x*dx_xz + dy10z*dz_xz) * v + y10 * tg * tgxz * dv;
                
                
                     // y11 = s2132pi * y * (5z^2 - 1)
                     const double dy11y = s2132pi * (5.*tgzz-1.0);
                     const double dy11z = s2132pi * 10.*tgyz;
                     dt11_xx[ig] = -(dy11y*dy_xx + dy11z*dz_xx) * v + y11 * tg * tgxx * dv;
                     dt11_yy[ig] = -(dy11y*dy_yy + dy11z*dz_yy) * v + y11 * tg * tgyy * dv;
                     dt11_zz[ig] = -(dy11y*dy_zz + dy11z*dz_zz) * v + y11 * tg * tgzz * dv;
                     dt11_xy[ig] = -(dy11y*dy_xy + dy11z*dz_xy) * v + y11 * tg * tgxy * dv;
                     dt11_yz[ig] = -(dy11y*dy_yz + dy11z*dz_yz) * v + y11 * tg * tgyz * dv;
                     dt11_xz[ig] = -(dy11y*dy_xz + dy11z*dz_xz) * v + y11 * tg * tgxz * dv;
                
                     // y12 = s1054pi * x * y * z;
                     const double dy12x = s1054pi * tgyz;
                     const double dy12y = s1054pi * tgxz;
                     const double dy12z = s1054pi * tgxy;
                     dt12_xx[ig] = -(dy12x*dx_xx + dy12y*dy_xx + dy12z*dz_xx) * v + y12 * tg * tgxx * dv;
                     dt12_yy[ig] = -(dy12x*dx_yy + dy12y*dy_yy + dy12z*dz_yy) * v + y12 * tg * tgyy * dv;
                     dt12_zz[ig] = -(dy12x*dx_zz + dy12y*dy_zz + dy12z*dz_zz) * v + y12 * tg * tgzz * dv;
                     dt12_xy[ig] = -(dy12x*dx_xy + dy12y*dy_xy + dy12z*dz_xy) * v + y12 * tg * tgxy * dv;
                     dt12_yz[ig] = -(dy12x*dx_yz + dy12y*dy_yz + dy12z*dz_yz) * v + y12 * tg * tgyz * dv;
                     dt12_xz[ig] = -(dy12x*dx_xz + dy12y*dy_xz + dy12z*dz_xz) * v + y12 * tg * tgxz * dv;
                
                     // y13 = s1054pi * 0.5 * z * (x^2 - y^2)
                     const double dy13x = s1054pi * tgxz;
                     const double dy13y = -s1054pi * tgyz;
                     const double dy13z = s1054pi * 0.5 * (tgxx - tgyy);
                     dt13_xx[ig] = -(dy13x*dx_xx + dy13y*dy_xx + dy13z*dz_xx) * v + y13 * tg * tgxx * dv;
                     dt13_yy[ig] = -(dy13x*dx_yy + dy13y*dy_yy + dy13z*dz_yy) * v + y13 * tg * tgyy * dv;
                     dt13_zz[ig] = -(dy13x*dx_zz + dy13y*dy_zz + dy13z*dz_zz) * v + y13 * tg * tgzz * dv;
                     dt13_xy[ig] = -(dy13x*dx_xy + dy13y*dy_xy + dy13z*dz_xy) * v + y13 * tg * tgxy * dv;
                     dt13_yz[ig] = -(dy13x*dx_yz + dy13y*dy_yz + dy13z*dz_yz) * v + y13 * tg * tgyz * dv;
                     dt13_xz[ig] = -(dy13x*dx_xz + dy13y*dy_xz + dy13z*dz_xz) * v + y13 * tg * tgxz * dv;
                
                     // y14 = s3532pi * (x^3 - 3*x*y^2)
                     const double dy14x = s3532pi * 3. * (tgxx - tgyy);
                     const double dy14y = -s3532pi * 6. * tgxy;
                     dt14_xx[ig] = -(dy14x*dx_xx + dy14y*dy_xx) * v + y14 * tg * tgxx * dv;
                     dt14_yy[ig] = -(dy14x*dx_yy + dy14y*dy_yy) * v + y14 * tg * tgyy * dv;
                     dt14_zz[ig] = -(dy14x*dx_zz + dy14y*dy_zz) * v + y14 * tg * tgzz * dv;
                     dt14_xy[ig] = -(dy14x*dx_xy + dy14y*dy_xy) * v + y14 * tg * tgxy * dv;
                     dt14_yz[ig] = -(dy14x*dx_yz + dy14y*dy_yz) * v + y14 * tg * tgyz * dv;
                     dt14_xz[ig] = -(dy14x*dx_xz + dy14y*dy_xz) * v + y14 * tg * tgxz * dv;
                
                     // y15 = s3532pi * (3*y*x^2 - y^3)
                     const double dy15x = s3532pi * 6.* tgxy;
                     const double dy15y = s3532pi * 3. * (tgxx - tgyy);
                     dt15_xx[ig] = -(dy15x*dx_xx + dy15y*dy_xx) * v + y15 * tg * tgxx * dv;
                     dt15_yy[ig] = -(dy15x*dx_yy + dy15y*dy_yy) * v + y15 * tg * tgyy * dv;
                     dt15_zz[ig] = -(dy15x*dx_zz + dy15y*dy_zz) * v + y15 * tg * tgzz * dv;
                     dt15_xy[ig] = -(dy15x*dx_xy + dy15y*dy_xy) * v + y15 * tg * tgxy * dv;
                     dt15_yz[ig] = -(dy15x*dx_yz + dy15y*dy_yz) * v + y15 * tg * tgyz * dv;
                     dt15_xz[ig] = -(dy15x*dx_xz + dy15y*dy_xz) * v + y15 * tg * tgxz * dv;
                  }
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
      
      if(s->nchannels() == 1){
	assert(ilm == s->nlm());
      }
    }
  }
  tmap["update_twnl"].stop();
}

////////////////////////////////////////////////////////////////////////////////
double NonLocalPotential::energy(SlaterDet& sd, bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl,
    vector<complex<double> >& veff) {
  const vector<double>& occ = sd.occ();
  const vector<double>& eig = sd.eig();
  const int ngwl = basis_.localsize();
  //const int mloc = basis_.maxlocalsize();
  const int mloc = sd.c().mloc();
  // define atom block size

  //const int na_block_size = 32;
  int namax = 0;
  for (int is=0; is<nsp; is++)
     if (atoms_.na(is) > namax) namax = atoms_.na(is);
  //const int na_block_size = (namax > 256) ? 256 : namax;
  const int na_block_size = (namax > 128) ? 128 : namax;

  vector<vector<double> > tau;
  atoms_.get_positions(tau,true);

  double enl = 0.0;
  double tsum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  
  if ( nspnl == 0 && !ultrasoft_) return 0.0;
  const double omega = basis_.cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  if (ultrasoft_) {
    tmap["usnl_calcbpsi"].start();
    sd.calc_betapsi();  // update betapsi w. latest wf
    tmap["usnl_calcbpsi"].stop();
    const int ngloc = cdbasis_->localsize();
    assert(veff.size() >= ngloc);
    int nsp = atoms_.nsp();
    for (int is=0; is<nsp; is++) {
      Species *s = atoms_.species_list[is];
      if (s->ultrasoft()) { 
        const int na = atoms_.na(is);
        const int naloc = atoms_.usloc_nat[is];
        const int naloc_t = atoms_.usloc_nat_t[is];
        const int naloc_max = atoms_.naloc_max[is];
        const int nqtot = s->nqtot();
        const int nbeta = s->nbeta();
        const int nbetalm = s->nbetalm();
        const int nstloc = sd.nstloc();

        const ComplexMatrix* betag = sd.betag(is);
        const ComplexMatrix* betapsi = sd.betapsi(is);
        const complex<double>* bg = betag->cvalptr();
        const complex<double>* bp = betapsi->cvalptr();
        int qsize = nqtot*na;
        vector<double> dmat(qsize);
        for (int i=0; i<qsize; i++)
          dmat[i] = 0.0;

        // calculate D_nm^I = D_nm^0 + cell_vol * SUM V_eff(G) * Q_nm^I(G)
        tmap["usnl_dmat"].start();
        if (highmem_) {
           #pragma omp parallel for
           for (int ibl = 0; ibl < naloc; ibl++) {
              int ia = atoms_.usloc_atind[is][ibl];     // ia = absolute atom index of qnmg
              for (int qind=0; qind < nqtot; qind++) {
                 const int qi0 = ibl*nqtot*ngloc + qind*ngloc;
                 const int dind = ia*nqtot + qind;
                 const int one = 1;

                 //double qv = 0.0;
                 //for (int ig=0; ig<ngloc; ig++)
                 //  qv += real(conj(qnmg_[is][qi0+ig])*veff[ig]);

                 //zdotc gives seg fault from C++
                 //complex<double> qv = zdotc((int*)&ngloc,&qnmg_[is][qi0],&one,&veff[0],&one);

                 complex<double> qv = complex<double>(0.0,0.0);
                 myzdotc(ngloc,&qnmg_[is][qi0],&veff[0],&qv);
                 
                 dmat[dind] = real(qv);  // missing omega_inv from qnmg FFT cancels out omega
              }
           }
        }
        else {
           #pragma omp parallel for
          for (int ibl = 0; ibl < naloc; ibl++) {
            int ia = atoms_.usloc_atind[is][ibl];     // ia = absolute atom index of qnmg
            for (int qind=0; qind < nqtot; qind++) {
              const int dind = ia*nqtot + qind;
              double qv = 0.0;
              for (int ig=0; ig<ngloc; ig++)
                qv += real(conj(sfactcd_[is][ibl*ngloc+ig]*qnmg_[is][qind*ngloc+ig])*veff[ig]);
              dmat[dind] = qv;  // missing omega_inv from qnmg FFT cancels out omega
            }
          }
        }
        tmap["usnl_dmat"].stop();
        
        // sum over local atoms
        ctxt_.dsum('r',qsize,1,&dmat[0],qsize);    
        // sum over plane waves
        ctxt_.dsum('c',qsize,1,&dmat[0],qsize);    

        // add dzero term
        //#pragma omp parallel for schedule(guided)
        for (int ia = 0; ia < na; ia++) {
          for (int qind=0; qind < nqtot; qind++) {
            const int dind = ia*nqtot + qind;
            dmat[dind] += s->dzero(qind);
          }
        }
        
        // accumulate Enl contribution
        tmap["usnl_enl"].start();
        if (betapsi->size() > 0) { 
           const int bp_mloc = betapsi->mloc(); // nbetalm*naloc_t
           for ( int lj=0; lj < sd.c().nblocks(); lj++ )
           {
              for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
              {
                 // global state index
                 const int nglobal = sd.c().j(lj,jj);
                 const int norig = lj*sd.c().nb()+jj;
                 const double occn = occ[nglobal];
                 for (int ibl = 0; ibl < naloc_t; ibl++) {
                    for (int qind=0; qind < nqtot; qind++) {
                       int lm1,lm2;
                       s->qind_to_betalm(qind,lm1,lm2);
                       int ind1 = bp_mloc*norig + ibl*nbetalm + lm1;
                       int ind2 = bp_mloc*norig + ibl*nbetalm + lm2;
                       double mult = 2.0;
                       if (lm1 == lm2) mult = 1.0;
                       double enltmp = omega_inv*mult*occn*s->dzero(qind)*real(conj(bp[ind1])*bp[ind2]);
                       enl += enltmp;
                    }
                 }
              }
           }
        }
        tmap["usnl_enl"].stop();
        
        if ( compute_hpsi ) {
          tmap["usnl_hpsi"].start();
          // calculate A_m^I = SUM_n D_nm^I <beta_n^I|psi>
          ComplexMatrix amat(betapsi->context(),betapsi->m(),betapsi->n(),betapsi->mb(),betapsi->nb());
          amat.clear();
          complex<double>* ap = amat.valptr();
          const int bp_mloc = betapsi->mloc();
          const int bp_nloc = betapsi->nloc();
          const int bg_mloc = betag->mloc();
          if (betapsi->size() > 0) { 
             //#pragma omp parallel for schedule(guided)
            for (int ibl = 0; ibl < naloc_t; ibl++) {
              int ia = atoms_.usloc_atind_t[is][ibl];
              for (int qind=0; qind < nqtot; qind++) {
                int lm1,lm2;
                s->qind_to_betalm(qind,lm1,lm2);
                for (int n=0; n<bp_nloc; n++) {
                  int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                  int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                  ap[ind1] += omega_inv*dmat[qind+ia*nqtot]*bp[ind2];
                  if (lm1 != lm2)
                    ap[ind2] += omega_inv*dmat[qind+ia*nqtot]*bp[ind1];
                }
              }
            }
          }
          
          if (sd.highmem()) {
            // add contribution to hpsi from amat*betapsi    
            ComplexMatrix& cp = dsd.c();
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "NLP.hpsi, species " << is << ", calling gemm of betag (" << betag->m() << " x " << betag->n() << ", " << betag->mloc() << " x " << betag->nloc() << ") and amat (" << amat.m() << " x " << amat.n() << ", " << amat.mloc() << " x " << amat.nloc() << ")" << endl;
#endif
            cp.gemm('n','n',1.0,*betag,amat,1.0);
          }
          else {
            // need to multiply betag by structure factor before gemm w. amat
            const int ngwl = basis_.localsize();
            const complex<double>* bgp = betag->cvalptr();

            // create temporary matrices to do parallel gemms for each local atom
            ComplexMatrix bgsf(ctxt_,betag->m(),betag->n(),betag->mb(),betag->nb());
            complex<double>* bgsfp = bgsf.valptr();
            int tbp_m = nbetalm*ctxt_.npcol();
            ComplexMatrix tbpsum(ctxt_,tbp_m,sd.c().n(),nbetalm,sd.c().nb());
            complex<double>* tbpp = tbpsum.valptr();
            for (int ibl = 0; ibl < naloc_max; ibl++) {
              vector<complex<double> > bpsum(nstloc*nbetalm);
              for (int i=0; i<bpsum.size(); i++)
                bpsum[i] = complex<double>(0.0,0.0);

              if (ibl < naloc) {
                for (int lm=0; lm<nbetalm; lm++)
                  for ( int ig = 0; ig < ngwl; ig++ )
                    bgsfp[lm*bg_mloc+ig] = bgp[lm*bg_mloc+ig]*sfactwf_[is][ibl*ngwl+ig];
              }
              if (ibl < naloc_t) { 
                int ia = atoms_.usloc_atind_t[is][ibl];
                // copy local atom's betapsi data to tbpsum (betapsi and betag
                // have atoms distributed differently (across cols vs. rows)
                // but matrix mult will use same local indices
                if (betapsi->size() > 0) { 
                  for (int n=0; n<nstloc; n++) {
                    for (int qind=0; qind<nqtot; qind++) {
                      int lm1 = s->qnm_lm1(qind);
                      int lm2 = s->qnm_lm2(qind);
                      int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                      int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                      int tind1 = nbetalm*n + lm1;
                      int tind2 = nbetalm*n + lm2;
                      bpsum[tind1] += dmat[qind+ia*nqtot]*bp[ind2];
                      if (lm1 != lm2) 
                        bpsum[tind2] += dmat[qind+ia*nqtot]*bp[ind1];
                    }
                  }
                }
              }
              if (tbpsum.size() >  0) {
                for (int i=0; i<nstloc*nbetalm; i++)
                  tbpp[i] = omega_inv*bpsum[i];
              }
              ComplexMatrix& cp = dsd.c();
#ifdef PRINTALL
      if (ctxt_.mype() == 0)
         cout << "NLP.hpsi, species " << is << ", calling gemm of bgsf (" << bgsf.m() << " x " << bgsf.n() << ", " << bgsf.mloc() << " x " << bgsf.nloc() << ") and tbpsum (" << tbpsum.m() << " x " << tbpsum.n() << ", " << tbpsum.mloc() << " x " << tbpsum.nloc() << ")" << endl;
#endif
              cp.gemm('n','n',1.0,bgsf,tbpsum,1.0);
            }
          }
          tmap["usnl_hpsi"].stop();
        }

        // ultrasoft force terms:  -dmat*dbpsisq/dR_I - veff*dQnm/dR_I*bpsisq 
        if ( compute_forces ) {
          tmap["usnl_fion"].start();
          valarray<double> tmpfion(3*na);
          tmpfion = 0.0;
          vector<double> ddmat(qsize);
          
          for ( int j = 0; j < 3; j++ ) {
            sd.calc_dbetapsi(j);
            const ComplexMatrix* dbetapsi = sd.dbetapsi(is);
            const complex<double>* dbpj = dbetapsi->cvalptr();
            if (dbetapsi->size() > 0) { 
               const int bp_mloc = dbetapsi->mloc(); // nbetalm*naloc_t
               for ( int lj=0; lj < sd.c().nblocks(); lj++ )
               {
                  for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                  {
                     // global state index
                     const int nglobal = sd.c().j(lj,jj);
                     const int norig = lj*sd.c().nb()+jj;
                     const double occn = occ[nglobal];
                     const double eign = eig[nglobal];
                     for (int ibl = 0; ibl < naloc_t; ibl++) {
                        int ia = atoms_.usloc_atind_t[is][ibl];
                        for (int qind=0; qind < nqtot; qind++) {
                           int lm1,lm2;
                           s->qind_to_betalm(qind,lm1,lm2);
                           int ind1 = bp_mloc*norig + ibl*nbetalm + lm1;
                           int ind2 = bp_mloc*norig + ibl*nbetalm + lm2;
                           double mult = 2.0;
                           if (lm1 == lm2) mult = 1.0;
                           double dmateq = dmat[qind+ia*nqtot] - eign*sd.qaug(is,qind);
                           tmpfion[3*ia+j] -= omega_inv * mult * occn * dmateq *
                               real(conj(dbpj[ind1])*bp[ind2] + conj(bp[ind1])*dbpj[ind2]);
                        }
                     }
                  }
               }
            }

            // calculate derivative of D_ij with respect to ion positions
            const double *const gxj = cdbasis_->gx_ptr(j);                      
            for (int i=0; i<qsize; i++)
              ddmat[i] = 0.0;
            if (highmem_) {
              for (int ibl = 0; ibl < naloc; ibl++) {
                int ia = atoms_.usloc_atind[is][ibl];     // ia = absolute atom index of qnmg
                for (int qind=0; qind < nqtot; qind++) {
                  const int qi0 = ibl*nqtot*ngloc + qind*ngloc;
                  double dqv = 0.0;
                  for (int ig=0; ig<ngloc; ig++) {
                    complex<double> dqnmg = conj(qnmg_[is][qi0+ig]*complex<double>(0.0,-gxj[ig]));
                    dqv += real(dqnmg*veff[ig]);
                  }
                  const int dind = ia*nqtot + qind;
                  ddmat[dind] = dqv;
                }
              }
            }
            else {
              for (int ibl = 0; ibl < naloc; ibl++) {
                int ia = atoms_.usloc_atind[is][ibl];     // ia = absolute atom index of qnmg
                for (int qind=0; qind < nqtot; qind++) {
                  double dqv = 0.0;
                  for (int ig=0; ig<ngloc; ig++) {
                    complex<double> dqnmg = conj(sfactcd_[is][ibl*ngloc+ig]*qnmg_[is][qind*ngloc+ig]*complex<double>(0.0,-gxj[ig]));
                    dqv += real(dqnmg*veff[ig]);
                  }
                  const int dind = ia*nqtot + qind;
                  ddmat[dind] = dqv;
                }                  
              }
            }

            // sum over local atoms
            ctxt_.dsum('r',qsize,1,&ddmat[0],qsize);    
            // sum over plane waves
            ctxt_.dsum('c',qsize,1,&ddmat[0],qsize);    
            
            if (dbetapsi->size() > 0) { 
              const int bp_mloc = dbetapsi->mloc(); // nbetalm*naloc_t
              if (highmem_) {
                for (int ibl = 0; ibl < naloc_t; ibl++) {
                  int ia = atoms_.usloc_atind_t[is][ibl];
                  for (int qind=0; qind < nqtot; qind++) {
                    const int dind = ia*nqtot+qind;
                    int lm1,lm2;
                    s->qind_to_betalm(qind,lm1,lm2);
                    double mult = 2.0;
                    if (lm1 == lm2) mult = 1.0;
                    for ( int lj=0; lj < sd.c().nblocks(); lj++ )
                    {
                       for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                       {
                          // global state index
                          const int nglobal = sd.c().j(lj,jj);
                          const int norig = lj*sd.c().nb()+jj;
                          const double occn = occ[nglobal];
                          int ind1 = bp_mloc*norig + ibl*nbetalm + lm1;
                          int ind2 = bp_mloc*norig + ibl*nbetalm + lm2;
                          tmpfion[3*ia+j] += omega_inv*mult*occn*ddmat[dind]*real(conj(bp[ind1])*bp[ind2]);
                       }
                    }
                  }
                }
              }
              else {
                for (int ibl = 0; ibl < naloc_t; ibl++) {
                  int ia = atoms_.usloc_atind_t[is][ibl];
                  for (int qind=0; qind < nqtot; qind++) {
                    const int dind = ia*nqtot+qind;
                    int lm1,lm2;
                    s->qind_to_betalm(qind,lm1,lm2);
                    double mult = 2.0;
                    if (lm1 == lm2) mult = 1.0;
                    for ( int lj=0; lj < sd.c().nblocks(); lj++ )
                    {
                       for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                       {
                          // global state index
                          const int nglobal = sd.c().j(lj,jj);
                          const int norig = lj*sd.c().nb()+jj;
                          const double occn = occ[nglobal];
                          int ind1 = bp_mloc*norig + ibl*nbetalm + lm1;
                          int ind2 = bp_mloc*norig + ibl*nbetalm + lm2;
                          tmpfion[3*ia+j] -= omega_inv*mult*occn*ddmat[dind]*real(conj(bp[ind1])*bp[ind2]);
                       }
                    }
                  }                  
                }
              }
            }
            
          }
          ctxt_.dsum(3*na,1,&tmpfion[0],3*na);

          //#pragma omp parallel for schedule(guided)
          for ( int ia = 0; ia < na; ia++ ) {
            fion[is][3*ia+0] = tmpfion[3*ia];
            fion[is][3*ia+1] = tmpfion[3*ia+1];
            fion[is][3*ia+2] = tmpfion[3*ia+2];
          }
          tmap["usnl_fion"].stop();          
        }
      }
    }
    // need to reduce enl contribution along columns (then rows, at end)
    ctxt_.dsum('c',1,1,&enl,1);    
  }

  // regular norm-conserving nonlocal potential projection
  for ( int is = 0; is < nsp; is++ ) {  
    Species *s = atoms_.species_list[is];
    if (!s->ultrasoft() && npr[is] > 0 ) { // species is is non-local, norm-conserving

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
        
      const int nstloc = sd.nstloc();
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
        tmap["comp_eigr"].start();                                   
        int k = 3;                                                   
        double mone = -1.0, zero = 0.0;
        char cn='n';

        // next line: const cast is ok since dgemm_ does not modify argument 
        double* kpgx = const_cast<double*>(basis_.kpgx_ptr(0));
        dgemm(&cn,&cn,(int*)&ngwl,(int*)&ia_block_size,&k,&mone,             
              kpgx,(int*)&ngwl, &tau[is][3*iastart],&k,                        
              &zero,&kpgr[0],(int*)&ngwl);                                     

        int len = ia_block_size * ngwl;                                      
#if HAVE_MASSV  
        vsincos(&skpgr[0],&ckpgr[0],&kpgr[0],&len);
#else
        #pragma omp parallel for
        for ( int i = 0; i < len; i++ ) {
          const double arg = kpgr[i];
          skpgr[i] = sin(arg);                                                 
          ckpgr[i] = cos(arg);                                                 
        }                                                                    
#endif
        tmap["comp_eigr"].stop();                                            

        // compute anl_loc                                                   
        tmap["comp_anl"].start();
        for ( int ipr = 0; ipr < npr[is]; ipr++ ) {
          // twnl[is][ig+ngwl*ipr]
          const double * t = &twnl[is][ngwl*ipr];                            
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
        tmap["comp_anl"].stop();                                             
                                                                             
        // array anl_loc is complete                                         
                                                                             
        // compute fnl[npra][nstloc] = anl^T * c                             
        double one=1.0;                                                      
        char ct='t';                                                         
        int twongwl = 2 * ngwl;                                              
        int twomloc = 2 * mloc;                                              
        int nprnaloc = ia_block_size * npr[is];                              
        int c_lda;
        const complex<double>* c = sd.c().cvalptr();                        

        if (basis_.real()) {
          c_lda = 2*sd.c().mloc();                                        
          tmap["fnl_gemm"].start();                                            
          dgemm(&ct,&cn,&nprnaloc,(int*)&nstloc,&twongwl,&one,
                &anl_loc_gamma[0],&twongwl, (double*)c, &c_lda,
                &zero,&fnl_loc_gamma[0],&nprnaloc);
          tmap["fnl_gemm"].stop();
        }
        else {
           c_lda = mloc;
          complex<double> zzero = complex<double>(0.0,0.0);
          complex<double> zone = complex<double>(1.0,0.0);
          char cc='c';
          tmap["fnl_gemm"].start();
          zgemm(&cc,&cn,&nprnaloc,(int*)&nstloc,(int*)&ngwl,&zone,
                &anl_loc[0],(int *)&ngwl, (complex<double> *)c, &c_lda,
                &zzero,&fnl_loc[0],&nprnaloc);
          tmap["fnl_gemm"].stop();
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
               (double*)c,&c_lda,&fnl_loc_gamma[0],&nprnaloc);
          }
        }

        tmap["fnl_allreduce"].start();                                       
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
        tmap["fnl_allreduce"].stop();                                        
                      
        // factor 2.0 in next line is: counting G, -G                        
        if (basis_.real())
          fnl_loc_gamma = 2.0 * fnl_buf_gamma;
        else
          fnl_loc = fnl_buf;
                                
        // accumulate Enl contribution                                       
        if (basis_.real()) {
          for ( int ipr = 0; ipr < npr[is]; ipr++ ) { 
            const double fac = wt[is][ipr] * omega_inv;                        
            for ( int lj=0; lj < sd.c().nblocks(); lj++ )
            {
               for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
               {
                  // global state index
                  const int nglobal = sd.c().j(lj,jj);
                  const int norig = lj*sd.c().nb()+jj;
                  const double facn = fac * occ[nglobal];
                  for ( int ia = 0; ia < ia_block_size; ia++ ) {
                     const int i = ia + ipr*ia_block_size + norig * nprnaloc;           
                     const double tmp = fnl_loc_gamma[i];
                     enl += facn * tmp * tmp;
                     fnl_loc_gamma[i] = fac * tmp;
                  }
               }
            }
          }
        }
        else {
          for ( int ipr = 0; ipr < npr[is]; ipr++ ) { 
            const double fac = wt[is][ipr] * omega_inv;                        
            for ( int lj=0; lj < sd.c().nblocks(); lj++ )
            {
               for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
               {
                  // global state index
                  const int nglobal = sd.c().j(lj,jj);
                  const int norig = lj*sd.c().nb()+jj;
                  const double facn = fac * occ[nglobal];
                  for ( int ia = 0; ia < ia_block_size; ia++ ) {
                     const int i = ia + ipr*ia_block_size + norig * nprnaloc;           
                     const complex<double> tmp = fnl_loc[i];                                 
                     enl += facn * norm(tmp);
                     fnl_loc[i] = fac * tmp;
                  }
               }
            }
          }
        }

        if ( compute_hpsi ) {
          tmap["enl_hpsi"].start();                                          
          // compute cp += anl * fnl                                         
          complex<double>* cp = dsd.c().valptr();                            

          int cp_lda;
          if (basis_.real()) {
            cp_lda = 2*dsd.c().mloc();
            dgemm(&cn,&cn,&twongwl,(int*)&nstloc,&nprnaloc,&one,
                  &anl_loc_gamma[0],&twongwl, &fnl_loc_gamma[0],&nprnaloc,
                  &one,(double*)cp, &cp_lda);
          }
          else {
            int cp_lda = dsd.c().mloc();
            complex<double> zone = complex<double>(1.0,0.0);
            zgemm(&cn,&cn,(int*)&ngwl,(int*)&nstloc,&nprnaloc,&zone,
                    &anl_loc[0],(int*)&ngwl, &fnl_loc[0],&nprnaloc,
                    &zone,(complex<double>*)cp, &cp_lda);
          }

          tmap["enl_hpsi"].stop();
        }                       
                                
        // ionic forces
         
        if ( compute_forces ) {
          tmap["enl_fion"].start();
          valarray<double> dfnl_loc_gamma;
          valarray<complex<double> > dfnl_loc;
          if (basis_.real()) 
            dfnl_loc_gamma.resize(npr[is]*na_block_size*nstloc);
          else
            dfnl_loc.resize(npr[is]*na_block_size*nstloc);

          for ( int j = 0; j < 3; j++ ) {
            const double *const kpgxj = basis_.kpgx_ptr(j);                      
                                
            // compute anl_loc  
            for ( int ipr = 0; ipr < npr[is]; ipr++ ) {
              // twnl[is][ig+ngwl*ipr]                                       
              const double * t = &twnl[is][ngwl*ipr];                        
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
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    const double tt = kpgxj[ig] * t[ig];                       
                    // Next lines: -i * ( a + ib ) = b - ia                  
                    a[2*ig]   =  tt * s[ig];                                 
                    a[2*ig+1] = -tt * c[ig];                                 
                  }             
                }               
                else if ( l == 1 ) {
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    // Next lines: (-i)**2 * ( a + ib ) = - a - ib           
                    const double tt = kpgxj[ig] * t[ig];                     
                    a[2*ig]   = -tt * c[ig];                                  
                    a[2*ig+1] = -tt * s[ig];                                  
                  }             
                }               
                else if ( l == 2 ) {
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    // Next lines: (-i) * - ( a + ib ) = i*(a+ib) = - b + ia 
                    const double tt = kpgxj[ig] * t[ig];                       
                    a[2*ig]   = -tt * s[ig];                                 
                    a[2*ig+1] =  tt * c[ig];                                 
                  }             
                }               
                else if ( l == 3 ) {
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    // Next lines: (-i)**4 * ( a + ib ) = a + ib
                    const double tt = kpgxj[ig] * t[ig];                       
                    a[2*ig]   =  tt * c[ig];                                 
                    a[2*ig+1] =  tt * s[ig];                                 
                  }             
                }               
              }                 
            } // ipr            
 
            // array anl_loc is complete                                     
                                
            // compute dfnl     
            // dfnl.gemm('t','n',2.0,anl,c_proxy,0.0);                       
            // compute dfnl[npra][nstloc] = anl^T * c                         
            const int nprnaloc = ia_block_size * npr[is];                    
            const complex<double>* c = sd.c().cvalptr();                    
            if (basis_.real()) {
              double one=1.0;
              char ct='t';        
              const int twongwl = 2 * ngwl;                                    
              dgemm(&ct,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&twongwl,&one, 
                    &anl_loc_gamma[0],(int*)&twongwl, (double*)c,(int*)&c_lda,       
                    &zero,&dfnl_loc_gamma[0],(int*)&nprnaloc);
            }
            else {
              complex<double> zzero = complex<double>(0.0,0.0);
              complex<double> zone = complex<double>(1.0,0.0);
              char cc='c';
              zgemm(&cc,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&ngwl,&zone, 
                    &anl_loc[0],(int*)&ngwl, (complex<double> *)c,(int*)&c_lda,       
                    &zzero,&dfnl_loc[0],(int*)&nprnaloc);                       
            }

            // Note: no need to correct for double counting of the           
            // G=0 component which is always zero                            

            // factor 2.0 in next line is: counting G, -G                    
            if (basis_.real()) 
              dfnl_loc_gamma *= 2.0;
                                
            // accumulate non-local contributions to forces                  
            if (basis_.real()) {
               for ( int ipr = 0; ipr < npr[is]; ipr++ )
               {
                  for ( int lj=0; lj < sd.c().nblocks(); lj++ )
                  {
                     for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                     {
                        // global state index
                        const int nglobal = sd.c().j(lj,jj);
                        const int norig = lj*sd.c().nb()+jj;
                        // Factor 2.0 in next line from derivative of |Fnl|^2        
                        const double facn = 2.0 * occ[nglobal];                    
                        for ( int ia = 0; ia < ia_block_size; ia++ ) {
                           const int ia_global = ia + iastart;                        
                           const int i = ia + ipr*ia_block_size + norig * nprnaloc;       
                           tmpfion[3*ia_global+j] -= facn *                           
                               fnl_loc_gamma[i] * dfnl_loc_gamma[i];        
                        }
                     }
                  }
               }
            }
            else {
               for ( int ipr = 0; ipr < npr[is]; ipr++ )
               {
                  for ( int lj=0; lj < sd.c().nblocks(); lj++ )
                  {
                     for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                     {
                        // global state index
                        const int nglobal = sd.c().j(lj,jj);
                        const int norig = lj*sd.c().nb()+jj;
                        // Factor 2.0 in next line from derivative of |Fnl|^2        
                        const double facn = 2.0 * occ[nglobal];                    
                        for ( int ia = 0; ia < ia_block_size; ia++ ) {
                           const int ia_global = ia + iastart;                        
                           const int i = ia + ipr*ia_block_size + norig * nprnaloc;       
                           //cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia            
                           //     << ",n=" << n << "]: " << fnl_loc[i] << endl;       
                           tmpfion[3*ia_global+j] -= facn *
                               real(fnl_loc[i] * conj(dfnl_loc[i]));
                        }
                     }
                  }
               }
            }
          } // j                
 
          tmap["enl_fion"].stop();                                           
        } // compute_forces     
                                
        if ( compute_stress ) {                       
          tmap["enl_sigma"].start();                                         
          valarray<double> dfnl_loc_gamma;
          valarray<complex<double> > dfnl_loc;
          if (basis_.real()) 
            dfnl_loc_gamma.resize(npr[is]*na_block_size*nstloc);
          else
            dfnl_loc.resize(npr[is]*na_block_size*nstloc);

          for ( int ij = 0; ij < 6; ij++ ) {
            // compute anl_loc  
            int ipr = 0;        
            while ( ipr < npr[is] ) {
              // twnl[is][ig+ngwl*ipr]                                       
              const int l = lproj[is][ipr];                                  
              if ( l == 0 ) {
                // dtwnl[is][ipr][ij][ngwl]                                  
                // index = ig + ngwl * ( ij + 6 * ipr))                      
                // ipr = iquad + nquad[is] * ilm, where ilm = 0              
                const double *const dt0 = &dtwnl[is][ngwl*(ij+6*ipr)];       
                for ( int ia = 0; ia < ia_block_size; ia++ ) {
                  double* a0;
                  if  (basis_.real()) 
                    a0 = &anl_loc_gamma[2*(ia+ipr*ia_block_size)*ngwl];      
                  else
                    a0 = (double*)&anl_loc[(ia+ipr*ia_block_size)*ngwl];      
                  const double* c = &ckpgr[ia*ngwl];                           
                  const double* s = &skpgr[ia*ngwl];                           
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    const double d0 = dt0[ig];                               
                    // Next lines: (-i)^0 * ( a + ib ) = a + ib                  
                    a0[2*ig]   =  d0 * c[ig];                                
                    a0[2*ig+1] =  d0 * s[ig];                                
                  }             
                }               
              }                 
              else if ( l == 1 ) {
                const int ipr1 = ipr;                                        
                const int ipr2 = ipr + 1;                                    
                const int ipr3 = ipr + 2;                                    
                // dtwnl[is][ipr][ij][ngwl]                                  
                // index = ig + ngwl * ( ij + 6 * iprx ))                    
                const double *dt1 = &dtwnl[is][ngwl*(ij+6*ipr1)];            
                const double *dt2 = &dtwnl[is][ngwl*(ij+6*ipr2)];            
                const double *dt3 = &dtwnl[is][ngwl*(ij+6*ipr3)];            
                for ( int ia = 0; ia < ia_block_size; ia++ ) {
                  double* a1;
                  double* a2;
                  double* a3;
                  if (basis_.real()) {
                    a1 = &anl_loc_gamma[2*(ia+ipr1*ia_block_size)*ngwl];     
                    a2 = &anl_loc_gamma[2*(ia+ipr2*ia_block_size)*ngwl];     
                    a3 = &anl_loc_gamma[2*(ia+ipr3*ia_block_size)*ngwl];     
                  }
                  else {
                    a1 = (double*)&anl_loc[(ia+ipr1*ia_block_size)*ngwl];     
                    a2 = (double*)&anl_loc[(ia+ipr2*ia_block_size)*ngwl];     
                    a3 = (double*)&anl_loc[(ia+ipr3*ia_block_size)*ngwl];     
                  }
                  const double* c = &ckpgr[ia*ngwl];                           
                  const double* s = &skpgr[ia*ngwl];                           
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    const double d1 = dt1[ig];                               
                    const double d2 = dt2[ig];                               
                    const double d3 = dt3[ig];                               
                    // Next line: (-i)^l factor is -i                        
                    // Next line: -i * eigr
                    // -i * (a+i*b) = b - i*a                                
                    const double tc = -c[ig]; //  -cosgr[ia][ig]             
                    const double ts =  s[ig]; //   singr[ia][ig]             
                    a1[2*ig]   = d1 * ts;                                    
                    a1[2*ig+1] = d1 * tc;                                    
                    a2[2*ig]   = d2 * ts;                                    
                    a2[2*ig+1] = d2 * tc;                                    
                    a3[2*ig]   = d3 * ts;                                    
                    a3[2*ig+1] = d3 * tc;                                    
                  }             
                }               
              }                  
              else if ( l == 2 ) {
                const int ipr4 = ipr;                                         
                const int ipr5 = ipr + 1;                                     
                const int ipr6 = ipr + 2;                                     
                const int ipr7 = ipr + 3;                                     
                const int ipr8 = ipr + 4;                                     
                // dtwnl[is][ipr][iquad][ij][ngwl]                            
                // index = ig + ngwl * ( ij + 6 * ( iquad + nquad[is] * ipr ))
                const double *dt4 = &dtwnl[is][ngwl*(ij+6*ipr4)];             
                const double *dt5 = &dtwnl[is][ngwl*(ij+6*ipr5)];             
                const double *dt6 = &dtwnl[is][ngwl*(ij+6*ipr6)];             
                const double *dt7 = &dtwnl[is][ngwl*(ij+6*ipr7)];             
                const double *dt8 = &dtwnl[is][ngwl*(ij+6*ipr8)];             
                for ( int ia = 0; ia < ia_block_size; ia++ ) {
                  double* a4;
                  double* a5;
                  double* a6;
                  double* a7;
                  double* a8;
                  if (basis_.real()) {
                    a4 = &anl_loc_gamma[2*(ia+ipr4*ia_block_size)*ngwl];      
                    a5 = &anl_loc_gamma[2*(ia+ipr5*ia_block_size)*ngwl];      
                    a6 = &anl_loc_gamma[2*(ia+ipr6*ia_block_size)*ngwl];      
                    a7 = &anl_loc_gamma[2*(ia+ipr7*ia_block_size)*ngwl];      
                    a8 = &anl_loc_gamma[2*(ia+ipr8*ia_block_size)*ngwl];      
                  }
                  else {
                    a4 = (double*)&anl_loc[(ia+ipr4*ia_block_size)*ngwl];      
                    a5 = (double*)&anl_loc[(ia+ipr5*ia_block_size)*ngwl];      
                    a6 = (double*)&anl_loc[(ia+ipr6*ia_block_size)*ngwl];      
                    a7 = (double*)&anl_loc[(ia+ipr7*ia_block_size)*ngwl];      
                    a8 = (double*)&anl_loc[(ia+ipr8*ia_block_size)*ngwl];      
                  }
                  const double* c = &ckpgr[ia*ngwl];                            
                  const double* s = &skpgr[ia*ngwl];                            
 
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    const double d4 = dt4[ig];                                
                    const double d5 = dt5[ig];                                
                    const double d6 = dt6[ig];                                
                    const double d7 = dt7[ig];                                
                    const double d8 = dt8[ig];                                
                    // Next lines: (-i)^2 * ( a + ib ) =  - ( a + ib )        
                    const double tc = -c[ig]; //  -cosgr[ia][ig]              
                    const double ts = -s[ig]; //  -singr[ia][ig]              
                    a4[2*ig]   = d4 * tc;                                     
                    a4[2*ig+1] = d4 * ts;                                     
                    a5[2*ig]   = d5 * tc;                                     
                    a5[2*ig+1] = d5 * ts;                                     
                    a6[2*ig]   = d6 * tc;                                     
                    a6[2*ig+1] = d6 * ts;                                     
                    a7[2*ig]   = d7 * tc;                                     
                    a7[2*ig+1] = d7 * ts;                                     
                    a8[2*ig]   = d8 * tc;                                     
                    a8[2*ig+1] = d8 * ts;                                     
                  }              
                }                
              }                  
              else if ( l == 3 ) {
                const int ipr9 = ipr;                                         
                const int ipr10 = ipr + 1;
                const int ipr11 = ipr + 2;
                const int ipr12 = ipr + 3;
                const int ipr13 = ipr + 4;
                const int ipr14 = ipr + 5;
                const int ipr15 = ipr + 6;

                // dtwnl[is][ipr][iquad][ij][ngwl]                            
                // index = ig + ngwl * ( ij + 6 * ( iquad + nquad[is] * ipr ))
                const double *dt9 = &dtwnl[is][ngwl*(ij+6*ipr9)];
                const double *dt10 = &dtwnl[is][ngwl*(ij+6*ipr10)];
                const double *dt11 = &dtwnl[is][ngwl*(ij+6*ipr11)];
                const double *dt12 = &dtwnl[is][ngwl*(ij+6*ipr12)];
                const double *dt13 = &dtwnl[is][ngwl*(ij+6*ipr13)];
                const double *dt14 = &dtwnl[is][ngwl*(ij+6*ipr14)];
                const double *dt15 = &dtwnl[is][ngwl*(ij+6*ipr15)];
                for ( int ia = 0; ia < ia_block_size; ia++ ) {
                  double* a9;
                  double* a10;
                  double* a11;
                  double* a12;
                  double* a13;
                  double* a14;
                  double* a15;
                  if (basis_.real()) {
                    a9 = &anl_loc_gamma[2*(ia+ipr9*ia_block_size)*ngwl];      
                    a10 = &anl_loc_gamma[2*(ia+ipr10*ia_block_size)*ngwl];      
                    a11 = &anl_loc_gamma[2*(ia+ipr11*ia_block_size)*ngwl];      
                    a12 = &anl_loc_gamma[2*(ia+ipr12*ia_block_size)*ngwl];      
                    a13 = &anl_loc_gamma[2*(ia+ipr13*ia_block_size)*ngwl];      
                    a14 = &anl_loc_gamma[2*(ia+ipr14*ia_block_size)*ngwl];      
                    a15 = &anl_loc_gamma[2*(ia+ipr15*ia_block_size)*ngwl];      
                  }
                  else {
                    a9 = (double*)&anl_loc[(ia+ipr9*ia_block_size)*ngwl];      
                    a10 = (double*)&anl_loc[(ia+ipr10*ia_block_size)*ngwl];      
                    a11 = (double*)&anl_loc[(ia+ipr11*ia_block_size)*ngwl];      
                    a12 = (double*)&anl_loc[(ia+ipr12*ia_block_size)*ngwl];      
                    a13 = (double*)&anl_loc[(ia+ipr13*ia_block_size)*ngwl];      
                    a14 = (double*)&anl_loc[(ia+ipr14*ia_block_size)*ngwl];      
                    a15 = (double*)&anl_loc[(ia+ipr15*ia_block_size)*ngwl];      
                  }
                  const double* c = &ckpgr[ia*ngwl];                            
                  const double* s = &skpgr[ia*ngwl];                            
 
                  for ( int ig = 0; ig < ngwl; ig++ ) {
                    const double d9 = dt9[ig];
                    const double d10 = dt10[ig];
                    const double d11 = dt11[ig];
                    const double d12 = dt12[ig];
                    const double d13 = dt13[ig];
                    const double d14 = dt14[ig];
                    const double d15 = dt15[ig];
                    // Next lines: (-i)^3 * ( a + ib ) =  -b + ia        
                    const double ts = -s[ig]; //  -singr[ia][ig]              
                    const double tc = c[ig]; //  cosgr[ia][ig]

                    a9[2*ig]   = d9 * ts;
                    a9[2*ig+1] = d9 * tc;
                    a10[2*ig]   = d10 * ts;
                    a10[2*ig+1] = d10 * tc;
                    a11[2*ig]   = d11 * ts;
                    a11[2*ig+1] = d11 * tc;
                    a12[2*ig]   = d12 * ts;
                    a12[2*ig+1] = d12 * tc;
                    a13[2*ig]   = d13 * ts;
                    a13[2*ig+1] = d13 * tc;
                    a14[2*ig]   = d14 * ts;
                    a14[2*ig+1] = d14 * tc;
                    a15[2*ig]   = d15 * ts;
                    a15[2*ig+1] = d15 * tc;
                  }              
                }                
              }
              else {                  
                assert(false);   
              } // l             
              ipr += 2*l+1;      
            } // while ipr       
      
            // array anl_loc is complete                                      
                                 
            // compute dfnl      
            // dfnl.gemm('t','n',2.0,anl,c_proxy,0.0);                        
            // compute dfnl[npra][nstloc] = anl^T * c                         
            const int nprnaloc = ia_block_size * npr[is];                     
            const complex<double>* c = sd.c().cvalptr();                     
            if (basis_.real()) {
              double one=1.0;      
              char ct='t';         
              const int twongwl = 2 * ngwl;                                     
              dgemm(&ct,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&twongwl,&one,  
                    &anl_loc_gamma[0],(int*)&twongwl, (double*)c,(int*)&c_lda,        
                    &zero,&dfnl_loc_gamma[0],(int*)&nprnaloc);                        
            }
            else {
              complex<double> zzero = complex<double>(0.0,0.0);
              complex<double> zone = complex<double>(1.0,0.0);
              char cc='c';
              zgemm(&cc,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&ngwl,&zone, 
                    &anl_loc[0],(int*)&ngwl, (complex<double> *)c,(int*)&c_lda,       
                    &zzero,&dfnl_loc[0],(int*)&nprnaloc);                       
            }

            // Note: no need to correct for double counting of the            
            // G=0 component which is always zero                             

            // factor 2.0 in next line is: counting G, -G                     
            if (basis_.real()) 
              dfnl_loc_gamma *= 2.0;
                                 
            // accumulate non-local contributions to sigma_ij                 
            if (basis_.real())
            {
               for ( int lj=0; lj < sd.c().nblocks(); lj++ )
               {
                  for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                  {
                     // global state index
                     const int nglobal = sd.c().j(lj,jj);
                     const int norig = lj*sd.c().nb()+jj;
                     // Factor 2.0 in next line from derivative of |Fnl|^2           
                     const double facn = 2.0 * occ[nglobal];                       
                     for ( int ipra = 0; ipra < npr[is]*ia_block_size; ipra++ ) {
                        const int i = ipra + norig * nprnaloc;                            
                        tsum[ij] += facn * fnl_loc_gamma[i] * dfnl_loc_gamma[i];
                     }                  
                  }
               }
            }
            else {
               for ( int lj=0; lj < sd.c().nblocks(); lj++ )
               {
                  for ( int jj=0; jj < sd.c().nbs(lj); jj++ )
                  {
                     // global state index
                     const int nglobal = sd.c().j(lj,jj);
                     const int norig = lj*sd.c().nb()+jj;
                     // Factor 2.0 in next line from derivative of |Fnl|^2           
                     const double facn = 2.0 * occ[nglobal];                       
                     for ( int ipra = 0; ipra < npr[is]*ia_block_size; ipra++ ) {
                        const int i = ipra + norig * nprnaloc;                            
                        tsum[ij] += facn * real(fnl_loc[i] * conj(dfnl_loc[i]));
                     }
                  }
               }
            }
          } // ij                
          tmap["enl_sigma"].stop();                                           
        } // compute_stress      
      } // ia_block              
     
      if ( compute_forces ) {
        ctxt_.dsum(3*na[is],1,&tmpfion[0],3*na[is]);
        //#pragma omp parallel for schedule(guided)
        for ( int ia = 0; ia < na[is]; ia++ ) {
          fion[is][3*ia+0] = tmpfion[3*ia];
          fion[is][3*ia+1] = tmpfion[3*ia+1];
          fion[is][3*ia+2] = tmpfion[3*ia+2];
        }
      }
    } // npr[is]>0
  } // is

  // reduction of enl across rows
  ctxt_.dsum('r',1,1,&enl,1);    

  sigma_enl = 0.0;
  if ( compute_stress ) {
    ctxt_.dsum(6,1,&tsum[0],6);
    sigma_enl[0] = ( enl + tsum[0] ) * omega_inv;
    sigma_enl[1] = ( enl + tsum[1] ) * omega_inv;
    sigma_enl[2] = ( enl + tsum[2] ) * omega_inv;
    sigma_enl[3] = + tsum[3] * omega_inv;
    sigma_enl[4] = + tsum[4] * omega_inv;
    sigma_enl[5] = + tsum[5] * omega_inv;
  }

  return enl;
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::update_usfns(SlaterDet& sd, Basis* cdbasis) {

  // update ultrasoft potential functions Q_nm^I(G), beta^I(G) whenever atoms
  // move or basis changes
  cdbasis_ = cdbasis;
  sd.init_usfns(&atoms_);
  //ewd:12-21-2011 sd.calc_betag();
  
  vector<vector<double> > tau;
  atoms_.get_positions(tau,true);
  
  // fill Q_nm(G), multiply structure factor term if highmem_
  const int ngloc = cdbasis_->localsize();
  const int ngwl = basis_.localsize();
  const int nsp = atoms_.nsp();
  qnmg_.resize(nsp);
  sfactcd_.resize(nsp);
  sfactwf_.resize(nsp);
  for (int is=0; is<nsp; is++) {
    Species *s = atoms_.species_list[is];
    if (s->ultrasoft()) { 
      int naloc = atoms_.usloc_nat[is];
      int nqtot = s->nqtot();
      int nbeta = s->nbeta();

      if (highmem_) {

         //ewd DEBUG
         if (ctxt_.mype() == 0)
            cout << "NonLocalPotential::update_usfns, qnmg_.size = " << nqtot*naloc*ngloc << ", naloc = " << naloc << ", nqtot = " << nqtot << ", ngloc = " << ngloc << endl;
      
        qnmg_[is].resize(nqtot*naloc*ngloc);
        vector<complex<double> > qnm;
        vector<double> qaug;
        s->calc_qnmg(cdbasis_,qnm,qaug);
        sd.set_qaug(is,qaug);          // store qaug in SlaterDet
      
        for (int ibl = 0; ibl < naloc; ibl++) {
          int ia = atoms_.usloc_atind[is][ibl];

          vector<double> ckpgr(ngloc);
          vector<double> skpgr(ngloc);
          const double *const kpgx = cdbasis_->kpgx_ptr(0);
          const double *const kpgy = cdbasis_->kpgx_ptr(1);
          const double *const kpgz = cdbasis_->kpgx_ptr(2);

          // calculate structure factor
          for ( int ig = 0; ig < ngloc; ig++ ) {
            const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
            skpgr[ig] = sin(arg);
            ckpgr[ig] = cos(arg);
          }

          for (int qind=0; qind < nqtot; qind++) {
            const int ind0 = ibl*nqtot*ngloc + qind*ngloc;
            for ( int ig = 0; ig < ngloc; ig++ )
              qnmg_[is][ind0+ig] = qnm[qind*ngloc+ig]*complex<double>(ckpgr[ig],-skpgr[ig]);
          }
        }        
      }
      else {
         //ewd DEBUG
         if (ctxt_.mype() == 0)
            cout << "NonLocalPotential::update_usfns, qnmg_.size = " << nqtot*ngloc << ", nqtot = " << nqtot << ", ngloc = " << ngloc << endl;
         
        qnmg_[is].resize(nqtot*ngloc);
        sfactcd_[is].resize(naloc*ngloc);
        sfactwf_[is].resize(naloc*ngwl);
        vector<double> qaug;
        s->calc_qnmg(cdbasis_,qnmg_[is],qaug);
        sd.set_qaug(is,qaug);          // store qaug in SlaterDet
    
        for (int ibl = 0; ibl < naloc; ibl++) {
          int ia = atoms_.usloc_atind[is][ibl];
          // calculate structure factor on cdbasis
          const double *const cdkpgx = cdbasis_->kpgx_ptr(0);
          const double *const cdkpgy = cdbasis_->kpgx_ptr(1);
          const double *const cdkpgz = cdbasis_->kpgx_ptr(2);
          for ( int ig = 0; ig < ngloc; ig++ ) {
            const double arg = tau[is][3*ia]*cdkpgx[ig] + tau[is][3*ia+1]*cdkpgy[ig] + tau[is][3*ia+2]*cdkpgz[ig];
            sfactcd_[is][ibl*ngloc+ig] = complex<double>(cos(arg),-sin(arg));
          }
        }

        for (int ibl = 0; ibl < naloc; ibl++) {
          int ia = atoms_.usloc_atind[is][ibl];

          // calculate structure factor on wf basis
          const double *const wfkpgx = basis_.kpgx_ptr(0);
          const double *const wfkpgy = basis_.kpgx_ptr(1);
          const double *const wfkpgz = basis_.kpgx_ptr(2);
          for ( int ig = 0; ig < ngwl; ig++ ) {
            const double arg = tau[is][3*ia]*wfkpgx[ig] + tau[is][3*ia+1]*wfkpgy[ig] + tau[is][3*ia+2]*wfkpgz[ig];
            sfactwf_[is][ibl*ngwl+ig] = complex<double>(cos(arg),-sin(arg));
          }
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::print_memory(ostream& os, int kmult, int kmultloc, double& totsum, double& locsum) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setprecision(3);
  PrintMem pm;
    
  double twnl_size = 0.;
  double twnl_locsize = 0.;
  double qnmg_size = 0.;
  double qnmg_locsize = 0.;

  const int nsp = atoms_.nsp();
  const int ngw = basis_.size();
  const int ngwl = basis_.localsize();
  for (int is=0; is<nsp; is++) {
     //twnl_size += (double)(7*npr[is]*ngw*sizeof(double));      // includes dtwnl
     //twnl_locsize += (double)(7*npr[is]*ngwl*sizeof(double));
     twnl_size += (double)(npr[is]*ngw*sizeof(double));
     twnl_locsize += (double)(npr[is]*ngwl*sizeof(double));
     Species *s = atoms_.species_list[is];
     if (s->ultrasoft()) { 
        int nqtot = s->nqtot();
        if (highmem_) {
           int na = atoms_.na(is);
           int naloc = atoms_.usloc_atind[is].size();
           qnmg_size += (double)(nqtot*na*ngw*sizeof(complex<double>));
           qnmg_locsize += (double)(nqtot*naloc*ngwl*sizeof(complex<double>));
        }
        else {
           qnmg_size += (double)(nqtot*ngw*sizeof(complex<double>));
           qnmg_locsize += (double)(nqtot*ngwl*sizeof(complex<double>));
        }
     }
  }

  qnmg_size *= kmult;
  qnmg_locsize *= kmultloc;
  twnl_size *= kmult;
  twnl_locsize *= kmultloc;
  totsum += qnmg_size + twnl_size;
  locsum += qnmg_locsize + twnl_locsize;
  
  if (ultrasoft_) {
    string qnmg_unit = pm.memunit(qnmg_size);
    string qnmg_locunit = pm.memunit(qnmg_locsize);
    os << "<!-- memory nlp.qnmg    :  " << setw(7) << qnmg_size << qnmg_unit << "  (" << qnmg_locsize << qnmg_locunit << " local) -->" << endl;
  }
  if (twnl_size > 0) { 
    string twnl_unit = pm.memunit(twnl_size);
    string twnl_locunit = pm.memunit(twnl_locsize);
    os << "<!-- memory nlp.twnl    :  " << setw(7) << twnl_size << twnl_unit << "  (" << twnl_locsize << twnl_locunit << " local) -->" << endl;
  }
}
  
