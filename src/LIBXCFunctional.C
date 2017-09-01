////////////////////////////////////////////////////////////////////////////////
//
// LIBXCFunctional.C
//
// Interface of Libxc Exchange-correlation energy and potential
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LIBXCFunctional.C,v  Yi Yao $


#include <cmath>
#include <cassert>
#include <vector>
#include "LIBXCFunctional.h"
#include <iostream>

/* include Libxc's header file */
#ifdef HAVE_LIBXC
#include "xc.h"   
#endif


void LIBXCFunctional::setxc(void)
{ 

#ifdef HAVE_LIBXC
  int numberofxc = func_ids_.size();
  double coeff_temp;
  switch (_xcfamily)
  {
    case XC_FAMILY_LDA:
    {
      if ( _np == 0 ) return;
      if ( _nspin == 1 )
      {
        double exc_tmp;
        double vxc_tmp;
        assert(rho != 0);
        assert(exc != 0);
        assert(vxc1 != 0);
        for ( int ir = 0; ir < _np; ir++ )
        {
          exc[ir] = 0.0;
          vxc1[ir] = 0.0;
        }
        for ( int ixc = 0; ixc < numberofxc; ixc ++ ) {
          xc_func_init(&func_, func_ids_[ixc], XC_UNPOLARIZED);
          coeff_temp = func_coeffs_[ixc];

          for ( int ir = 0; ir < _np; ir++ )
          {
            if ( rho[ir] < 0 ) continue;
            xc_lda_exc_vxc(&func_, 1, &rho[ir], &exc_tmp, &vxc_tmp);
            exc[ir] += exc_tmp * coeff_temp;
            vxc1[ir] += vxc_tmp * coeff_temp;
            //xc_unpolarized(rho[ir],exc[ir],vxc1[ir]);
          }
          xc_func_end(&func_);
        }
      }
      else
      {
        // spin polarized
        assert(rho_up != 0);
        assert(rho_dn != 0);
        assert(exc != 0);
        assert(vxc1_up != 0);
        assert(vxc1_dn != 0);
        double exc_spin_temp;
        double vxc_spin_temp[2];
        double rho_spin_temp[2];

        for ( int ir = 0; ir < _np; ir++ )
        {
          exc[ir] = 0.0;
          vxc1_up[ir] = 0.0;
          vxc1_dn[ir] = 0.0;
        }

        for ( int ixc = 0; ixc < numberofxc; ixc ++ ) {
          xc_func_init(&func_, func_ids_[ixc], XC_POLARIZED);
          coeff_temp = func_coeffs_[ixc];

          for ( int ir = 0; ir < _np; ir++ )
          {
            rho_spin_temp[0] = rho_up[ir];
            rho_spin_temp[1] = rho_dn[ir];
            if (rho_spin_temp[0] < 0 ) rho_spin_temp[0] = 0.0;
            if (rho_spin_temp[1] < 0 ) rho_spin_temp[1] = 0.0;
            if (rho_spin_temp[0] + rho_spin_temp[1] < 0) continue;

            xc_lda_exc_vxc(&func_, 1, &rho_spin_temp[0], &exc_spin_temp, &vxc_spin_temp[0]);

            exc[ir] += exc_spin_temp * coeff_temp;
            vxc1_up[ir] += vxc_spin_temp[0] * coeff_temp;
            vxc1_dn[ir] += vxc_spin_temp[1] * coeff_temp;
            //std::cout << exc_spin_temp << " " << vxc_spin_temp[0] << " " << vxc_spin_temp[1] << std::endl;
          }
          xc_func_end(&func_);
        }

      }
    }  // end case XC_FAMILY_LDA
    break;
    case XC_FAMILY_GGA:
    {
      if ( _np == 0 ) return;
      if ( _nspin == 1 )
      {
        assert( rho != 0 );
        assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
        assert( exc != 0 );
        assert( vxc1 != 0 );
        assert( vxc2 != 0 );
    
        for ( int ir = 0; ir < _np; ir++ )
        {
          exc[ir] = 0.0;
          vxc1[ir] = 0.0;
          vxc2[ir] = 0.0;
        }
        double exc_temp, vxc_rho_temp, vxc_sigma_temp, sigma;
        for ( int ixc = 0; ixc < numberofxc; ixc ++ ) {

          xc_func_init(&func_, func_ids_[ixc], XC_UNPOLARIZED);
          coeff_temp = func_coeffs_[ixc];
          switch ( func_families_[ixc] )
          {
            case XC_FAMILY_LDA:
              for ( int ir = 0; ir < _np; ir++ )
              {
                if ( rho[ir] < 0 ) continue;
                xc_lda_exc_vxc(&func_, 1, &rho[ir], &exc_temp, &vxc_rho_temp);
                exc[ir] += exc_temp * coeff_temp;
                vxc1[ir] += vxc_rho_temp * coeff_temp;
              }
            break;
            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
              for ( int ir = 0; ir < _np; ir++ )
              {
                if (rho[ir] < 1.e-18 ) continue;
                //if (rho[ir] >= 1.e-18) {
                sigma =  (grad_rho[0][ir]*grad_rho[0][ir] +
                          grad_rho[1][ir]*grad_rho[1][ir] +
                          grad_rho[2][ir]*grad_rho[2][ir] );
                xc_gga_exc_vxc(&func_, 1, &rho[ir], &sigma, &exc_temp, &vxc_rho_temp, &vxc_sigma_temp);
                exc[ir] += exc_temp * coeff_temp;
                vxc1[ir] += vxc_rho_temp * coeff_temp;
                vxc2[ir] += - vxc_sigma_temp * 2.0 * coeff_temp;
                //}
              }
            break;
          }

          xc_func_end(&func_);
        }
      }
      else
      {
        // spin polarized
        assert( rho_up != 0 );
        assert( rho_dn != 0 );
        assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
        assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
        assert( exc_up != 0 );
        assert( exc_dn != 0 );
        assert( vxc1_up != 0 );
        assert( vxc1_dn != 0 );
        assert( vxc2_upup != 0 );
        assert( vxc2_updn != 0 );
        assert( vxc2_dnup != 0 );
        assert( vxc2_dndn != 0 );

        double exc_spin_temp;
        double vxc_spin_temp[2];
        double vsigma_spin_temp[3];
        double rho_spin_temp[2];
        double sigma_spin_temp[3];

        for ( int ir = 0; ir < _np; ir++ )
        {
          exc_up[ir] = 0.0;
          exc_dn[ir] = 0.0;
          vxc1_up[ir] = 0.0;
          vxc1_dn[ir] = 0.0;
          vxc2_upup[ir] = 0.0;
          vxc2_dndn[ir] = 0.0;
          vxc2_updn[ir] = 0.0;
          vxc2_dnup[ir] = 0.0;
        }
        for ( int ixc = 0; ixc < numberofxc; ixc ++ ) {
          xc_func_init(&func_, func_ids_[ixc], XC_POLARIZED);
          coeff_temp = func_coeffs_[ixc];
          switch ( func_families_[ixc] )
          {
            case XC_FAMILY_LDA:
            {
              for ( int ir = 0; ir < _np; ir++ )
              {
                rho_spin_temp[0] = rho_up[ir];
                rho_spin_temp[1] = rho_dn[ir];
                if (rho_spin_temp[0] < 0 ) rho_spin_temp[0] = 0.0;
                if (rho_spin_temp[1] < 0 ) rho_spin_temp[1] = 0.0;
                if (rho_spin_temp[0] + rho_spin_temp[1] < 0) continue;

                xc_lda_exc_vxc(&func_, 1, &rho_spin_temp[0], &exc_spin_temp, &vxc_spin_temp[0]);

                exc_up[ir] += exc_spin_temp * coeff_temp;
                exc_dn[ir] += exc_spin_temp * coeff_temp;
                vxc1_up[ir] += vxc_spin_temp[0] * coeff_temp;
                vxc1_dn[ir] += vxc_spin_temp[1] * coeff_temp;
              }
            }
            break;

            case XC_FAMILY_GGA:
            case XC_FAMILY_HYB_GGA:
            {
              for ( int ir = 0; ir < _np; ir++ )
              {
                if (rho_up[ir] < 1.e-10 && rho_dn[ir] < 1.e-10 ) continue;
                if (rho_up[ir] + rho_dn[ir] < 1.e-10 ) continue;
                double grx_up = grad_rho_up[0][ir];
                double gry_up = grad_rho_up[1][ir];
                double grz_up = grad_rho_up[2][ir];
                double grx_dn = grad_rho_dn[0][ir];
                double gry_dn = grad_rho_dn[1][ir];
                double grz_dn = grad_rho_dn[2][ir];
                double grad_up2 = grx_up*grx_up + gry_up*gry_up + grz_up*grz_up;
                double grad_dn2 = grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn;
                double grad_up_grad_dn = grx_up*grx_dn + gry_up*gry_dn + grz_up*grz_dn;
                //double grx = grx_up + grx_dn;
                //double gry = gry_up + gry_dn;
                //double grz = grz_up + grz_dn;
                //double grad2 = grx*grx + gry*gry + grz*grz;

                rho_spin_temp[0] = rho_up[ir];
                rho_spin_temp[1] = rho_dn[ir];
                if ( rho_spin_temp[0] < 0.0 ) rho_spin_temp[0] = 0.0;
                if ( rho_spin_temp[1] < 0.0 ) rho_spin_temp[1] = 0.0;
                sigma_spin_temp[0] = grad_up2;
                sigma_spin_temp[1] = grad_up_grad_dn;
                sigma_spin_temp[2] = grad_dn2;

                xc_gga_exc_vxc(&func_, 1, &rho_spin_temp[0], &sigma_spin_temp[0],
                               &exc_spin_temp, &vxc_spin_temp[0], &vsigma_spin_temp[0]);


                exc_up[ir] += exc_spin_temp * coeff_temp;
                exc_dn[ir] += exc_spin_temp * coeff_temp;
                vxc1_up[ir] += vxc_spin_temp[0] * coeff_temp;
                vxc1_dn[ir] += vxc_spin_temp[1] * coeff_temp;

                // I am not sure about these lines. They are not the same as native functionals.
                // However for PBE and BLYP the actual results are identical in real calculations.
                // YY
                vxc2_upup[ir] += -vsigma_spin_temp[0] * 2.0 * coeff_temp;
                vxc2_updn[ir] += -vsigma_spin_temp[1] * 1.0 * coeff_temp;
                vxc2_dnup[ir] += -vsigma_spin_temp[1] * 1.0 * coeff_temp;
                vxc2_dndn[ir] += -vsigma_spin_temp[2] * 2.0 * coeff_temp;
              }
            }
            break;
          }
          xc_func_end(&func_);
        }
      }

    } // end case XC_FAMILY_GGA
    break;

  } // end switch
#endif
}

