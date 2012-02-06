/*******************************************************************************
 *
 * spline.h
 *
 ******************************************************************************/

#define SPLINE_FLAT_BC 0.0       /* Flat boundary condition (y'=0) */
#define SPLINE_NATURAL_BC 1.e31  /* Natural boundary condition (Y"=0) */

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint (double *xa, double *ya, double *y2a, int n, double x, double *y);
void splintd (double *xa, double *ya, double *y2a, 
 int n, double x, double *y, double *dy);
