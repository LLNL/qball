////////////////////////////////////////////////////////////////////////////////
//
// LineMinimizer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LineMinimizer.h,v 1.2 2009/03/27 00:53:24 draeger1 Exp $

#ifndef LINEMINIMIZER_H
#define LINEMINIMIZER_H

#include <iostream>
#include <cassert>
#include <cmath>

class LineMinimizer
{
  private:

  double sigma1_;
  double alpha_start, alpha_low, alpha_high;
  double f_low, fp_low, f_high, fp_high;
  double f0, fp0;
  bool bracket, first_use;
  bool debug_print;
  double interpolate(void)
  {
    // select value of alpha in [alpha_low,alpha_high] using the
    // values of f_low, f_high, fp_low, fp_high

    if ( fp_low * fp_high < 0.0 )
    {
      // quadratic interpolation
      if ( debug_print )
      {
        std::cout << " LineMinimizer: quadratic interpolation: "
                  << std::endl;
        std::cout << " LineMinimizer: alpha_low/high: "
                  << alpha_low << " " << alpha_high << std::endl;
        std::cout << " LineMinimizer: f_low/high: "
                  << f_low << " " << f_high << std::endl;
        std::cout << " LineMinimizer: fp_low/high: "
                  << fp_low << " " << fp_high << std::endl;
      }
      const double dalpha = alpha_high - alpha_low;
      return alpha_low -
             0.5 * ( fp_low * dalpha * dalpha ) /
                   ( f_high - f_low - fp_low * dalpha );
    }
    else
    {
      // bisection
      return 0.5 * ( alpha_low + alpha_high );
    }
  }

  public:

  LineMinimizer(void) : sigma1_(0.1), alpha_start(1.0), first_use(true),
   debug_print(false) {}
  void reset(void) { first_use = true; }
  double sigma1(void) const { return sigma1_; }
  void set_sigma1(double s1) { sigma1_ = s1; }
  double newalpha(double alpha, double f, double fp)
  {
    // compute a new alpha
    // fp = d_f / d_alpha
    if ( first_use )
    {
      alpha_low = 0.0;
      f_low = f;
      fp_low = fp;
      f0 = f;
      fp0 = fp;
      first_use = false;
      bracket = false;
      return alpha_start;
    }

    if ( !bracket )
    {
      if ( debug_print )
      {
        std::cout << " LineMinimizer: not bracketing, alpha_low = "
                  << alpha_low << std::endl;
        std::cout << " LineMinimizer: not bracketing, f_low = "
                  << f_low << std::endl;
        std::cout << " LineMinimizer: not bracketing, fp_low = "
                  << fp_low  << std::endl;
      }
      if ( fp*fp_low < 0.0 || (f > f0 + fp0 * sigma1_ * alpha) || f > f_low )
      {
        alpha_high = alpha;
        f_high = f;
        fp_high = fp;
        bracket = true;
        if ( f_low > f_high )
        {
          double tmp;
          tmp = alpha_low; alpha_low = alpha_high; alpha_high = tmp;
          tmp = f_low; f_low = f_high; f_high = tmp;
          tmp = fp_low; fp_low = fp_high; fp_high = tmp;
        }
        return interpolate();
      }
      else
      {
        const double delta_alpha = 1.1 * ( alpha - alpha_low );
        // increase alpha by at least 0.1
        const double alpha_t = alpha + std::max(delta_alpha,0.1);
        if ( debug_print )
        {
          std::cout << " LineMinimizer: expanding, alpha = " << alpha
                    << " alpha_low = "
                    << alpha_low << " alpha_t = " << alpha_t  << std::endl;
        }
        alpha_low = alpha;
        f_low = f;
        fp_low = fp;
        return alpha_t;
      }
    }
    else
    {
      // we are in the bracketing phase
      if ( debug_print )
      {
        std::cout << " LineMinimizer: bracketing, alpha_low/high = "
                  << alpha_low << "  " << alpha_high << std::endl;
        std::cout << " LineMinimizer: bracketing, f_low/high = "
                  << f_low << " " << f_high << std::endl;
        std::cout << " LineMinimizer: bracketing, fp_low/high = "
                  << fp_low << " " << fp_high << std::endl;
      }
      if ( fp*fp_low < 0.0 )
      {
        // bracket(alpha, alpha_low)
        alpha_high = alpha;
        f_high = f;
        fp_high = fp;
        bracket = true;
        if ( f_low > f_high )
        {
          double tmp;
          tmp = alpha_low; alpha_low = alpha_high; alpha_high = tmp;
          tmp = f_low; f_low = f_high; f_high = tmp;
          tmp = fp_low; fp_low = fp_high; fp_high = tmp;
        }
        return interpolate();
      }
      else
      {
        // bracket(alpha, alpha_high);
        alpha_low = alpha;
        f_low = f;
        fp_low = fp;
        bracket = true;
        if ( f_low > f_high )
        {
          double tmp;
          tmp = alpha_low; alpha_low = alpha_high; alpha_high = tmp;
          tmp = f_low; f_low = f_high; f_high = tmp;
          tmp = fp_low; fp_low = fp_high; fp_high = tmp;
        }
        return interpolate();
      }
    }
  }
};
#endif
