//Program: EBE_fluctuations_1.0
//Author: Christopher Plumberg
//Date: January 27, 2014

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
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

#include "corr_fcn_from_toyS.h"
#include "gauss.h"
#include "integration_routines.h"

//int read_corrfn(int length, vector<corr_fcn_data>* corrfn_ptr);
double S_function(double r, double phi, double eta, double tau, void * params_ptr);
double eta_t(double r, double phi, void * params_ptr);
phys_params get_parameters();
int Cal_correlationfunction_1D_from_toyS();
int Cal_correlationfunction_3D_from_toyS();
void Output_Correlationfunction_1D();
void Output_Correlationfunction_3D();
void initialize_corr_fcns();
double corr_fcn_integrand_re(vector<double>* xpts, void * params_ptr);
double corr_fcn_integrand_im(vector<double>* xpts, void * params_ptr);

bool output_to_screen = true;
bool output_to_file = false;

size_t n;
string input_filename;
string input_filepath;
//vector<corr_fcn_data> * corrfn_ptr;
//ostringstream corrfn_stream;

int main(int argc, char *argv[])
{
int success;
xi_ptr = new vector<double> (order);
wi_ptr = new vector<double> (order);
xi_0pinf_ptr = new vector<double> (order);
wi_0pinf_ptr = new vector<double> (order);
xi_minfpinf_ptr = new vector<double> (order);
wi_minfpinf_ptr = new vector<double> (order);
xiq_0pinf_ptr = new vector<double> (order);
wiq_0pinf_ptr = new vector<double> (order);
xiq_minfpinf_ptr = new vector<double> (order);
wiq_minfpinf_ptr = new vector<double> (order);

lower_limits_vec = new vector<double> (4);
upper_limits_vec = new vector<double> (4);
integration_intervals_vec = new vector<int> (4);

//set limits and integration intervals
(*lower_limits_vec)[0] = 0.;
(*upper_limits_vec)[0] = 0.;
(*integration_intervals_vec)[0] = 1;	//r
(*lower_limits_vec)[1] = -PI;
(*upper_limits_vec)[1] = PI;
(*integration_intervals_vec)[1] = 0;	//phi
(*lower_limits_vec)[2] = 0.;
(*upper_limits_vec)[2] = 0.;
(*integration_intervals_vec)[2] = 2;	//eta
(*lower_limits_vec)[3] = 0.;
(*upper_limits_vec)[3] = 0.;
(*integration_intervals_vec)[3] = 1;	//tau

//take arguments from command-line
if (argc == 1) { //if no arguments included in command-line, run with defaults
	eps_2_bar_cmd = eps_2_bar_default;
	eps_3_bar_cmd = eps_3_bar_default;
	v_2_bar_cmd = v_2_bar_default;
	v_3_bar_cmd = v_3_bar_default;
	psi_2_bar_cmd = psi_2_bar_default;
	psi_3_bar_cmd = psi_3_bar_default;
	tau_f_cmd = tau_f_default;
	Delta_tau_cmd = Delta_tau_default;
	K_perp_cmd = K_perp_default;
	eta_0_cmd = eta_0_default;
	Delta_eta_cmd = Delta_eta_default;
	eta_f_cmd = eta_f_default;
	Rad_cmd = Rad_default;
	Phi_K_cmd = Phi_K_default;
	T0_cmd = T0_default;
	flow_angle_cmd = flow_angle_default;
}
else if (is_multiple(argc , 2)) {
	//if incorrect number of arguments, end program and return (1)
	cerr << "Incorrect number of arguments: expected (function_call) (int param_key) (double param_val)" << endl
	     << "param_key: 1 - eps_2_bar, 2 - eps_3_bar, 3 - v_2_bar, 4 - v_3_bar" << endl
	     << "5 - psi_2_bar, 6 - psi_3_bar, 7 - tau_f, 8 - Delta_tau, 9 - K_perp" << endl
	     << "10 - eta_0, 11 - Delta_eta, 12 - eta_f, 13 - Rad, 14 - Phi_K, 15 - T0, 16 - flow_angle" << endl;

	return (1);
}
else {
	eps_2_bar_cmd = eps_2_bar_default;
	eps_3_bar_cmd = eps_3_bar_default;
	v_2_bar_cmd = v_2_bar_default;
	v_3_bar_cmd = v_3_bar_default;
	psi_2_bar_cmd = psi_2_bar_default;
	psi_3_bar_cmd = psi_3_bar_default;
	tau_f_cmd = tau_f_default;
	Delta_tau_cmd = Delta_tau_default;
	K_perp_cmd = K_perp_default;
	eta_0_cmd = eta_0_default;
	Delta_eta_cmd = Delta_eta_default;
	eta_f_cmd = eta_f_default;
	Rad_cmd = Rad_default;
	Phi_K_cmd = Phi_K_default;
	T0_cmd = T0_default;
	flow_angle_cmd = flow_angle_default;

int arg_index = 0;
do {
	switch (atoi(argv[arg_index+1]))
	{
		case 1:
			eps_2_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 2:
			eps_3_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 3:
			v_2_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 4:
			v_3_bar_cmd = atof(argv[arg_index+2]);
			break;
		case 5:
			psi_2_bar_cmd = atof(argv[arg_index+2])*PI;
			break;
		case 6:
			psi_3_bar_cmd = atof(argv[arg_index+2])*PI;
			break;
		case 7:
			tau_f_cmd = atof(argv[arg_index+2]);
			break;
		case 8:
			Delta_tau_cmd = atof(argv[arg_index+2]);
			break;
		case 9:
			K_perp_cmd = atof(argv[arg_index+2]);
			break;
		case 10:
			eta_0_cmd = atof(argv[arg_index+2]);
			break;
		case 11:
			Delta_eta_cmd = atof(argv[arg_index+2]);
			break;
		case 12:
			eta_f_cmd = atof(argv[arg_index+2]);
			break;
		case 13:
			Rad_cmd = atof(argv[arg_index+2]);
			break;
		case 14:
			Phi_K_cmd = atof(argv[arg_index+2]);
			break;
		case 15:
			T0_cmd = atof(argv[arg_index+2]);
			break;
		case 16:
			flow_angle_cmd = atof(argv[arg_index+2]);
			break;
	}

arg_index += 2;  //advance to the next two input arguments
  } while (arg_index < argc-1);
}
if (!output_to_file && !output_to_screen)
{
	cerr << "No output mode selected.  Exiting..." << endl;
	return (1);
}

if (output_to_screen)
{
	cout << "Initial parameters are:" << endl
		<< "eps_2_bar = " << eps_2_bar_cmd << endl
		<< "eps_3_bar = " << eps_3_bar_cmd << endl
		<< "v_2_bar = " << v_2_bar_cmd << endl
		<< "v_3_bar = " << v_3_bar_cmd << endl
		<< "psi_2_bar = " << psi_2_bar_cmd << endl
		<< "psi_3_bar = " << psi_3_bar_cmd << endl
		<< "tau_f = " << tau_f_cmd << endl
		<< "Delta_tau = " << Delta_tau_cmd << endl
		<< "eta_0 = " << eta_0_cmd << endl
		<< "Delta_eta = " << Delta_eta_cmd << endl
		<< "eta_f = " << eta_f_cmd << endl
		<< "Rad = " << Rad_cmd << endl
		<< "T_0 = " << T0_cmd << endl
		<< "K_perp = " << K_perp_cmd << endl
		<< "Phi_K = " << Phi_K_cmd << endl << endl;
}

/*
ostringstream test_stream;
test_stream << path << "/test.dat";
if (fexists( test_stream.str().c_str() ) )
{
	cout << "exists!!!" << endl;
	int fileindex = 2;
	do
	{
		test_stream.str("");
		test_stream.clear();
		test_stream << path << "/test" << fileindex << ".dat";
		fileindex++;
	} while (fexists( test_stream.str().c_str() ));
}
cout << "final stringstream was " << test_stream.str() << endl;
   ofstream testout;
   testout.open(test_stream.str().c_str());
if(1) return (1);
*/

/**************************************************************************/
/*************************MAIN PART OF THE PROGRAM*************************/
/**************************************************************************/

	//initialize model parameters
	get_parameters();

	//set the points needed to perform the integrations
	set_gaussian_points();

	//set up correlation functions
	initialize_corr_fcns();

	//do actual integration to get and fill correlation functions
	success = Cal_correlationfunction_1D_from_toyS();
	//success = Cal_correlationfunction_3D_from_toyS();

	//output results
	Output_Correlationfunction_1D();
	//Output_Correlationfunction_3D();

/**************************************************************************/
/***********************END MAIN PART OF THE PROGRAM***********************/
/**************************************************************************/
	return (0);
}

/**************************************************************************/

phys_params get_parameters()
{
	//struct phys_params inits;

	//Define parameters for function to be integrated here
	inits.T_0 = T0_cmd;
	inits.Y_rapidity = 0.;
	inits.Rad = Rad_cmd;
	inits.eta_f = eta_f_cmd;
	inits.tau_f = tau_f_cmd;
	inits.K_perp = K_perp_cmd;
	inits.M_perp = sqrt(m_pion*m_pion + inits.K_perp*inits.K_perp);
	inits.Phi_K = Phi_K_cmd;
	inits.eta_0 = eta_0_cmd;
	inits.Delta_tau = Delta_tau_cmd;
	inits.Delta_eta = Delta_eta_cmd;
	inits.beta_perp = inits.K_perp/inits.M_perp;
	inits.beta_long = 0.;

	inits.eps_2_bar = eps_2_bar_cmd;
	inits.eps_3_bar = eps_3_bar_cmd;
	inits.v_2_bar = v_2_bar_cmd;
	inits.v_3_bar = v_3_bar_cmd;
	inits.psi_2_bar = psi_2_bar_cmd;
	inits.psi_3_bar = psi_3_bar_cmd;

	return (inits);
}


/*double S_function(double r, double phi, double eta, double tau, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double M_perp = params.M_perp;
	double T_0 = params.T_0;
	double eta_0 = params.eta_0;
	double Y_rapidity = params.Y_rapidity;
	double Rad = params.Rad;
	double tau_f = params.tau_f;
	double K_perp = params.K_perp;
	double Delta_tau = params.Delta_tau;
	double Delta_eta = params.Delta_eta;
	double Phi_K = params.Phi_K;

	double eps_2_bar = params.eps_2_bar;
	double psi_2_bar = params.psi_2_bar;
	double eps_3_bar = params.eps_3_bar;
	double psi_3_bar = params.psi_3_bar;

	//Expression for function is long and complicated,
	//so break it up into constituent terms
	double term0, term1, term2, term3, term4;
	term0 = (tau-tau_f)*(tau-tau_f) / (2.*Delta_tau*Delta_tau);
	//term0 = 0.;
	term1 = (eta-eta_0)*(eta-eta_0) / (2.*Delta_eta*Delta_eta);
	//term1 = 0.;
	term2 = (M_perp/T_0) * cosh(eta - Y_rapidity) * cosh(eta_t(r, phi, &params));
	//term2 = 0.;
	term3 = (K_perp/T_0) * (cos(phi) * cos(Phi_K) + sin(phi) * sin(Phi_K)) * sinh(eta_t(r, phi, &params));
	//term3 = 0.;
	term4 = (r*r)/(2.*Rad*Rad) * (1.+2.*eps_3_bar*( cos(3.*phi) * cos(3.*psi_3_bar) + sin(3.*phi) * sin(3.*psi_3_bar) )    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					+2.*eps_2_bar*( cos(2.*phi) * cos(2.*psi_2_bar) + sin(2.*phi) * sin(2.*psi_2_bar) ));  //changed minus signs to plus signs //5/16/2013
	//term4 = 0.;

	return (exp(-term0 - term1 - term2 + term3 - term4));
	//return (1.);
}*/

double S_function(double r, double phi, double eta, double tau, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double Rad = params.Rad;

	return (heaviside_theta(Rad - r) * heaviside_theta(r));
}

double eta_t(double r, double phi, void * params_ptr)
{
	phys_params params = * (struct phys_params *) params_ptr;

	double v_2_bar = params.v_2_bar;
	double psi_2_bar = params.psi_2_bar;
	double eta_f = params.eta_f;
	double Rad = params.Rad;
	double v_3_bar = params.v_3_bar;
	double psi_3_bar = params.psi_3_bar;

	return (eta_f * (r/Rad) * (1.+2.*v_3_bar*( cos(3.*phi) * cos(3.*psi_3_bar) + sin(3.*phi) * sin(3.*psi_3_bar) )
					+2.*v_2_bar*( cos(2.*phi) * cos(2.*psi_2_bar) + sin(2.*phi) * sin(2.*psi_2_bar) )));
}

int Cal_correlationfunction_1D_from_toyS()
{
cout << "K_T = "  << inits.K_perp << "\t" << " and K_phi = " << inits.Phi_K << endl;
double integ1 = 0., integ2 = 0., sum = 0., error = 1.e-4;

my_q_out = 0.;
my_q_side = 0.;
my_q_long = 0.;
double spectra = integrate(&corr_fcn_integrand_re, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec)/(hbarC*hbarC*hbarC);

   cout << "generating the 1d slices of the correlation function along q_out, q_side, and q_long direction... " << endl;
		  
   for(int i = 0; i < qnpts; i++)
   {
      cout << "calculating q_mu = " << q_out[i] << endl;
      double values[3];
      for (int ops = 0; ops < 3; ops++)
         values[ops] = 0.0;
      for (int l = 0; l < 3; l++)
      {
         switch (l)
         {
            case 0:
            {
               my_q_out  = q_out[i];
               my_q_side = 0.0e0;
               my_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               my_q_out  = 0.0e0;
               my_q_side = q_side[i];
               my_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               my_q_out  = 0.0e0;
               my_q_side = 0.0e0;
               my_q_long = q_long[i];
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

         integ1 = integrate(&corr_fcn_integrand_re, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec);
         integ2 = integrate(&corr_fcn_integrand_im, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec);
         integ1 = integ1/hbarC/hbarC/hbarC/spectra;
         integ2 = integ2/hbarC/hbarC/hbarC/spectra;
//cerr << integ1 << endl;
//cerr << integ2 << endl;
//cerr << spectra << endl;
         double localvalue = integ1*integ1+integ2*integ2;
         values[l] = localvalue;
      }
      Correl_1D_out[i]  = values[0];
      Correl_1D_side[i] = values[1];
      Correl_1D_long[i] = values[2];
      Correl_1D_out_err[i] = error;
      Correl_1D_side_err[i] = error;
      Correl_1D_long_err[i] = error;
cout << "RESULTS: " << scientific << setprecision(7) << setw(15)
                   << q_out[i] << "  " << Correl_1D_out[i] << "  " << Correl_1D_out_err[i] << "  "
                   << q_side[i] << "  " << Correl_1D_side[i] << "  " << Correl_1D_side_err[i] << "  "
                   << q_long[i] << "  " << Correl_1D_long[i] << "  " << Correl_1D_long_err[i] 
                   << endl;
   }
   return(0);
}

int Cal_correlationfunction_3D_from_toyS()
{
cout << "K_T = "  << inits.K_perp << "\t" << " and K_phi = " << inits.Phi_K << endl;
double integ1 = 0., integ2 = 0., sum = 0., error = 1.e-4;

my_q_out = 0.;
my_q_side = 0.;
my_q_long = 0.;
double spectra = integrate(&corr_fcn_integrand_re, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec)/hbarC/hbarC/hbarC;

   cout << "generating correlation function in 3D... " << endl;
   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      //double local_q_out = q_out[i];
      double my_q_out = q_out[i];
      //cout << "q_out = " << local_q_out << endl;
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         //double local_q_side = q_side[j];
         double my_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            //double local_q_long = q_long[k];
            double my_q_long = q_long[k];

            integ1 = integrate(&corr_fcn_integrand_re, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec);
            integ2 = integrate(&corr_fcn_integrand_im, &inits, lower_limits_vec, upper_limits_vec, integration_intervals_vec);
            integ1 = integ1/hbarC/hbarC/hbarC/spectra;
            integ2 = integ2/hbarC/hbarC/hbarC/spectra;
            sum = integ1*integ1+integ2*integ2;
            Correl_3D[i][j][k] = sum;
            Correl_3D_err[i][j][k] = error;
         }
      }
   }
   return (0);
}

void initialize_corr_fcns()
{
   q_out = new double [qnpts];
   q_side = new double [qnpts];
   q_long = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      q_out[i] = init_q + (double)i * delta_q;
      q_side[i] = init_q + (double)i * delta_q;
      q_long[i] = init_q + (double)i * delta_q;
   }

   Correl_1D_out = new double [qnpts];
   Correl_1D_side = new double [qnpts];
   Correl_1D_long = new double [qnpts];
   Correl_1D_out_err = new double [qnpts];
   Correl_1D_side_err = new double [qnpts];
   Correl_1D_long_err = new double [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      Correl_1D_out[i] = 0.0;
      Correl_1D_side[i] = 0.0;
      Correl_1D_long[i] = 0.0;
      Correl_1D_out_err[i] = 0.0;
      Correl_1D_side_err[i] = 0.0;
      Correl_1D_long_err[i] = 0.0;
   }
   
   Correl_3D = new double** [qnpts];
   Correl_3D_err = new double** [qnpts];
   for(int i=0; i<qnpts; i++)
   {
      Correl_3D[i] = new double* [qnpts];
      Correl_3D_err[i] = new double* [qnpts];
      for(int j=0; j<qnpts; j++)
         Correl_3D[i][j] = new double [qnpts];
      for(int j=0; j<qnpts; j++)
         Correl_3D_err[i][j] = new double [qnpts];
   }
   for(int i=0; i<qnpts; i++)
      for(int j=0; j<qnpts; j++)
         for(int k=0; k<qnpts; k++)
         {
            Correl_3D[i][j][k] = 0.0;
            Correl_3D_err[i][j][k] = 0.0;
         }
return;
}

void Output_Correlationfunction_1D()
{
if (output_to_file)
   {
   ostringstream oCorrelfun_1D_stream;
   oCorrelfun_1D_stream << path << "/correlfunct1D" << "_" << particle_name << "_kt_" << inits.K_perp << "_phi_" << inits.Phi_K << ".dat";
if (fexists( oCorrelfun_1D_stream.str().c_str() ) )
{
	int fileindex = 2;
	do
	{
		oCorrelfun_1D_stream.str("");
		oCorrelfun_1D_stream.clear();
		oCorrelfun_1D_stream << path << "/correlfunct1D" << "_" << particle_name << "_kt_" << inits.K_perp << "_phi_" << inits.Phi_K << "_" << fileindex << ".dat";
		fileindex++;
	} while (fexists( oCorrelfun_1D_stream.str().c_str() ));
}
   ofstream oCorrelfun_1D;
   oCorrelfun_1D.open(oCorrelfun_1D_stream.str().c_str());
   for(int i=0; i < qnpts; i++)
     oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                   << q_out[i] << "  " << Correl_1D_out[i] << "  " << Correl_1D_out_err[i] << "  "
                   << q_side[i] << "  " << Correl_1D_side[i] << "  " << Correl_1D_side_err[i] << "  "
                   << q_long[i] << "  " << Correl_1D_long[i] << "  " << Correl_1D_long_err[i] 
                   << endl;
   }
if (output_to_screen)
   {
   for(int i=0; i < qnpts; i++)
   	cout << scientific << setprecision(7) << setw(15)
                   << q_out[i] << "  " << Correl_1D_out[i] << "  " << Correl_1D_out_err[i] << "  "
                   << q_side[i] << "  " << Correl_1D_side[i] << "  " << Correl_1D_side_err[i] << "  "
                   << q_long[i] << "  " << Correl_1D_long[i] << "  " << Correl_1D_long_err[i] 
                   << endl;
   }
   return;
}

void Output_Correlationfunction_3D()
{
if (output_to_file)
   {
   ostringstream oCorrelfun_3D_stream;
   oCorrelfun_3D_stream << path << "/correlfunct3D" << "_" << particle_name << "_kt_" << inits.K_perp << "_phi_" << inits.Phi_K << ".dat";
   ofstream oCorrelfun_3D;
   oCorrelfun_3D.open(oCorrelfun_3D_stream.str().c_str());
   for(int i=0; i < qnpts; i++)
      for(int j=0; j < qnpts; j++)
         for(int k=0; k < qnpts; k++)
              oCorrelfun_3D << scientific << setprecision(7) << setw(15)
                            << q_out[i] << "  " << q_side[j] << "  " 
                            << q_long[k] << "  " << Correl_3D[i][j][k] << "  "
                            << Correl_3D_err[i][j][k] << endl;
   }
if (output_to_screen)
   {
   for(int i=0; i < qnpts; i++)
      for(int j=0; j < qnpts; j++)
         for(int k=0; k < qnpts; k++)
         	cout << scientific << setprecision(7) << setw(15)
                            << q_out[i] << "  " << q_side[j] << "  " 
                            << q_long[k] << "  " << Correl_3D[i][j][k] << "  "
                            << Correl_3D_err[i][j][k] << endl;
   }
   return;
}

double corr_fcn_integrand_re(vector<double>* xpts, void * params_ptr)
{
	double rpt = (*xpts)[0];
	double phipt = (*xpts)[1];
	double etapt = (*xpts)[2];
	double taupt = (*xpts)[3];

	phys_params params = * (struct phys_params *) params_ptr;
	//double K_phi = params.Phi_K;
	double my_K_T = params.K_perp;
	//double K_y = 0.;
	double cosK_phi = cos(params.Phi_K);
	double sinK_phi = sin(params.Phi_K);

	double xsi = my_K_T*my_K_T + m_pion*m_pion + (my_q_out*my_q_out + my_q_side*my_q_side + my_q_long*my_q_long)/4.0;  //Set Xsi
	double E1sq = xsi + my_K_T*my_q_out;
	double E2sq = xsi - my_K_T*my_q_out;
	double qt = sqrt(E1sq) - sqrt(E2sq);
	double qx = my_q_out*cosK_phi - my_q_side*sinK_phi;
	double qy = my_q_side*cosK_phi + my_q_out*sinK_phi;
	double qz = my_q_long;

	//double tpt = taupt*cosh(etapt);
	double xpt = rpt*cos(phipt);
	double ypt = rpt*sin(phipt);
	//double zpt = taupt*sinh(etapt);
	double tpt = 0.;  //just using this for check
	double zpt = 0.;  //just using this for check
	double ss = S_function(rpt, phipt, etapt, taupt, params_ptr);
	double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;

/*cout << endl << endl;
cout << rpt << endl;
cout << phipt << endl;
cout << etapt << endl;
cout << taupt << endl;
cout << xsi << endl;
cout << qt << endl;
cout << qx << endl;
cout << qy << endl;
cout << qz << endl;
cout << rpt << endl;
cout << xpt << endl;
cout << ypt << endl;
cout << zpt << endl;
cout << ss << endl;
cout << arg << endl;
cout << endl << endl;*/

return (ss * cos(arg));
}

double corr_fcn_integrand_im(vector<double>* xpts, void * params_ptr)
{
	double rpt = (*xpts)[0];
	double phipt = (*xpts)[1];
	double etapt = (*xpts)[2];
	double taupt = (*xpts)[3];

	phys_params params = * (struct phys_params *) params_ptr;
	//double K_phi = params.Phi_K;
	double my_K_T = params.K_perp;
	//double K_y = 0.;
	double cosK_phi = cos(params.Phi_K);
	double sinK_phi = sin(params.Phi_K);

	double xsi = my_K_T*my_K_T + m_pion*m_pion + (my_q_out*my_q_out + my_q_side*my_q_side + my_q_long*my_q_long)/4.0;  //Set Xsi
	double E1sq = xsi + my_K_T*my_q_out;
	double E2sq = xsi - my_K_T*my_q_out;
	double qt = sqrt(E1sq) - sqrt(E2sq);
	double qx = my_q_out*cosK_phi - my_q_side*sinK_phi;
	double qy = my_q_side*cosK_phi + my_q_out*sinK_phi;
	double qz = my_q_long;

	//double tpt = taupt*cosh(etapt);
	double xpt = rpt*cos(phipt);
	double ypt = rpt*sin(phipt);
	//double zpt = taupt*sinh(etapt);
	double tpt = 0.;  //just using this for check
	double zpt = 0.;  //just using this for check
	double ss = S_function(rpt, phipt, etapt, taupt, params_ptr);
	double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;

return (ss * sin(arg));
}

//End of file
