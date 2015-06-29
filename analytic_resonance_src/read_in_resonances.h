#ifndef READ_IN_RESONANCES_H
#define READ_IN_RESONANCES_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <string>
#include <time.h>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <limits>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_matrix.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_linalg.h>

using namespace std;

#include "parameters.h"

int read_in_resonances(resonance_info * resonances);

//End of file
#endif
