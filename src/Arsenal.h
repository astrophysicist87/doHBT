#ifndef ARSENAL_H
#define ARSENAL_H

using namespace std;

unsigned long int random_seed();

//miscellaneous functions needed for interpolation routines
long binarySearch(double * A, int length, double value, bool skip_out_of_range);
void get_1D_derivatives(double * x, double * f, double * derivatives, int length, double default_edge_fill = 0.0);
void get_2D_derivatives(double * x, double * y, double ** f, double ** f1, double ** f2, double ** f12, int x_length, int y_length, double default_edge_fill = 0.0);
void bcucof(double * y, double * y1, double * y2, double * y12, double d1, double d2, double ** c);
void bcuint(double * y, double * y1, double * y2, double * y12, double x1l, double x1u, double x2l, double x2u, double x1, double x2, double &ansy, double &ansy1, double &ansy2);
void polint(double xa[], double ya[], long n, double x, double *y, double *dy);
void polin2(double * x1a, double * x2a, double ** ya, long m, long n, double x1, double x2, double *y, double *dy);

//interpolation routines
double interpLinearDirect(double * x, double * y, double x0, long size);
double interpBiLinearDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size);
double interpCubicDirect(double * x, double * y, double x0, long size);
double interpBiCubicDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size);
double interpPolyDirect(double * x, double * y, double x0, long size);
double interpBiPolyDirect(double * x, double * y, double ** z, double x0, double y0, long x_size, long y_size);

#endif
