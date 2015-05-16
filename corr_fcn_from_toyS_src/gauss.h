#ifndef GAUSS_H
#define GAUSS_H

//  file: gauss.h
//
//  This program determines points and weights for Gaussian quadrature
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      04-Jan-2004  original version, for 780.20 Computational Physics
//      07-Mar-2006  minor upgrades
//
//  Notes:  
//   * compile with:  "g++ -Wall -c gauss.cpp"
//   * adapted from: "Projects in Computational Physics" by Landau and Paez  
//             copyrighted by John Wiley and Sons, New York               
//             code copyrighted by RH Landau                           
// 
//************************************************************************

// include files
#include <cmath>

using namespace std;

void gauss (int npts, int job, double a, double b, double x[], double w[]);
void set_gaussian_points();

//************************************************************************
void
gauss (int npts, int job, double a, double b, double xpts[], double weights[])
{
  //     npts     number of points                                       
  //     job = 0  rescaling uniformly between (a,b)                      
  //           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
  //           2  for integral (a,inf) with 50% inside (a,b+2a)          
  //     xpts, weights     output grid points and weights.

  const double pi = M_PI;
  const double eps = 3.e-15;	// limit for accuracy 

  double t = 0, t1 = 0, p1 = 0, p2 = 0, p3 = 0, pp = 0;

  int m = (npts + 1) / 2;
  for (int i = 1; i <= m; i++)
    {
      t = cos (pi * (i - 0.25) / (npts + 0.5));
      t1 = 1;
      while ((fabs (t - t1)) >= eps)
	{
	  p1 = 1.0;
	  p2 = 0.0;
	  for (int j = 1; j <= npts; j++)
	    {
	      p3 = p2;
	      p2 = p1;
	      p1 = ((2 * j - 1) * t * p2 - (j - 1) * p3) / j;
	    }
	  pp = npts * (t * p1 - p2) / (t * t - 1);
	  t1 = t;
	  t = t1 - p1 / pp;
	}
      xpts[i - 1] = -t;
      xpts[npts - i] = t;
      weights[i - 1] = 2.0 / ((1 - t * t) * pp * pp);
      weights[npts - i] = weights[i - 1];
    }

  if (job == 0)		// rescaling uniformly between (a,b) 
    {
      for (int i = 0; i < npts; i++)
	{
	  xpts[i] = xpts[i] * (b - a) / 2.0 + (b + a) / 2.0;
	  weights[i] = weights[i] * (b - a) / 2.0;
	}
    }
  if (job == 1)		// integral (0,b) with 50% points inside (0, ab/(a+b))
    {
      for (int i = 0; i < npts; i++)
	{
	  t = (b + a) - (b - a) * xpts[i];
	  xpts[i] = a * b * (1 + xpts[i]) / t;
	  weights[i] = weights[i] * 2 * a * b * b / (t * t);
	}
    }
  if (job == 2)		// integral (a,inf) with 50% inside (a,b+2a)
    {
      for (int i = 0; i < npts; i++)
	{
	  t = 1.0 - xpts[i];
	  xpts[i] = (b * xpts[i] + b + a + a) / t;
	  weights[i] = weights[i] * 2 * (a + b) / (t * t);
	}
    }
}

void set_gaussian_points()
{
	//calculate zeroes of legendre polynomials (abscissae) and corresponding weights
	const double midpoint = 1.;
	double xpts[order];
	double wts[order];
	gauss(order,3,-1,1,xpts,wts);

	for(int i = 0; i < order; i++) {
		//finite
		(*xi_ptr)[i] = xpts[i];
		(*wi_ptr)[i] = wts[i];
//if (i==0) cout << (*xi_ptr)[i] << "\t" << (*wi_ptr)[i] << endl;

		//half-infinite
		(*xi_0pinf_ptr)[i] = midpoint * (1. + xpts[i]) / (1. - xpts[i]);
		(*wi_0pinf_ptr)[i] = 2. * midpoint * wts[i] / ( (1. - xpts[i]) * (1. - xpts[i]) );
//if (i==0) cout << (*xi_0pinf_ptr)[i] << "\t" << (*wi_0pinf_ptr)[i] << endl;

		//full-infinite
		(*xi_minfpinf_ptr)[i] = midpoint * xpts[i] / (1. - xpts[i] * xpts[i]);
		(*wi_minfpinf_ptr)[i] = midpoint * wts[i] * (1. + xpts[i] * xpts[i]) / ( (1. - xpts[i] * xpts[i]) * (1. - xpts[i] * xpts[i]) );
//if (i==0) cout << (*xi_minfpinf_ptr)[i] << "\t" << (*wi_minfpinf_ptr)[i] << endl;

//do it for q-integrations too
	const double midpointq = 1.;
		//half-infinite
		(*xiq_0pinf_ptr)[i] = midpointq * (1. + xpts[i]) / (1. - xpts[i]);
		(*wiq_0pinf_ptr)[i] = 2. * midpointq * wts[i] / ( (1. - xpts[i]) * (1. - xpts[i]) );
//if (i==0) cout << (*xi_0pinf_ptr)[i] << "\t" << (*wi_0pinf_ptr)[i] << endl;

		//full-infinite
		(*xiq_minfpinf_ptr)[i] = midpointq * xpts[i] / (1. - xpts[i] * xpts[i]);
		(*wiq_minfpinf_ptr)[i] = midpointq * wts[i] * (1. + xpts[i] * xpts[i]) / ( (1. - xpts[i] * xpts[i]) * (1. - xpts[i] * xpts[i]) );
//if (i==0) cout << (*xi_minfpinf_ptr)[i] << "\t" << (*wi_minfpinf_ptr)[i] << endl;
	}

	//cout<< "Computed Gaussian abscissae and corresponding weights..." << endl;
}

//End of file

#endif
