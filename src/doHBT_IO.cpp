#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "doHBT.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

void doHBT::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << path << folderindex << "/HBTradii_ev" << folderindex << ".dat";
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << ".dat";
	ofstream outputHBT(filename_stream_HBT.str().c_str(), ios::app);
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << path << folderindex << "/HBTradii_cfs_ev" << folderindex << ".dat";
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str(), ios::app);
	ostringstream filename_stream_S;
	//filename_stream_S << path << folderindex << "/Sourcefunction_variances.dat";
	filename_stream_S << global_path << "/Sourcefunction_variances.dat";
	ofstream output_Svars(filename_stream_S.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		outputHBTcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
			<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
			<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
			<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		outputHBT << folderindex << "  " << K_T[iKT] << "  " << K_phi[iKphi]
			<< "  " << R2_side[iKT][iKphi] << "  " << R2_out[iKT][iKphi]
			<< "  " << R2_outside[iKT][iKphi] << "  " << R2_long[iKT][iKphi]
			<< "  " << R2_sidelong[iKT][iKphi] << "  " << R2_outlong[iKT][iKphi] << endl;
         	output_Svars << scientific << setprecision(8) << setw(15) 
			<< K_T[iKT] << "   " << K_phi[iKphi] << "   " << S_func[iKT][iKphi] << "   "
			<< xs_S[iKT][iKphi]  << "   " << xo_S[iKT][iKphi]  << "   " << xl_S[iKT][iKphi]  << "   " << t_S[iKT][iKphi]  << "   "
			<< xs_t_S[iKT][iKphi]  << "   " << xo_t_S[iKT][iKphi]  << "   " << xl_t_S[iKT][iKphi]  << "   " << xo_xs_S[iKT][iKphi]  << "   "
			<< xl_xs_S[iKT][iKphi]  << "   " << xo_xl_S[iKT][iKphi]  << "   " << xs2_S[iKT][iKphi]  << "   " << xo2_S[iKT][iKphi]  << "   "
			<< xl2_S[iKT][iKphi]  << "   " << t2_S[iKT][iKphi] << endl;
	}
}

	outputHBT.close();
	outputHBTcoeffs.close();
	output_Svars.close();

	return;
}

void doHBT::Output_Svars_results(int folderindex)
{
	ostringstream filename_stream_Scos;
	//filename_stream_Scos << path << folderindex << "/Sourcefunction_variances_cfs_COS.dat";
	filename_stream_Scos << global_path << "/Sourcefunction_variances_cfs_COS.dat";
	ofstream output_SvarscoeffsCOS(filename_stream_Scos.str().c_str(), ios::app);
	ostringstream filename_stream_Ssin;
	//filename_stream_Ssin << path << folderindex << "/Sourcefunction_variances_cfs_SIN.dat";
	filename_stream_Ssin << global_path << "/Sourcefunction_variances_cfs_SIN.dat";
	ofstream output_SvarscoeffsSIN(filename_stream_Ssin.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
//only outputting shifted, i.e., "xs2_cos" is shorthand
//for the cosine terms in <xs^2>-<xs>^2, etc.
         	output_SvarscoeffsCOS << scientific << setprecision(8)
			<< folderindex << "   " << K_T[iKT] << "   " << Morder << "   "
			<< xs_t_cos[iKT][Morder]  << "   " << xo_t_cos[iKT][Morder]  << "   " << xl_t_cos[iKT][Morder]  << "   " << xo_xs_cos[iKT][Morder]  << "   "
			<< xl_xs_cos[iKT][Morder]  << "   " << xo_xl_cos[iKT][Morder]  << "   " << xs2_cos[iKT][Morder]  << "   " << xo2_cos[iKT][Morder]  << "   "
			<< xl2_cos[iKT][Morder]  << "   " << t2_cos[iKT][Morder] << endl;
         	output_SvarscoeffsSIN << scientific << setprecision(8)
			<< folderindex << "   " << K_T[iKT] << "   " << Morder << "   "
			<< xs_t_sin[iKT][Morder]  << "   " << xo_t_sin[iKT][Morder]  << "   " << xl_t_sin[iKT][Morder]  << "   " << xo_xs_sin[iKT][Morder]  << "   "
			<< xl_xs_sin[iKT][Morder]  << "   " << xo_xl_sin[iKT][Morder]  << "   " << xs2_sin[iKT][Morder]  << "   " << xo2_sin[iKT][Morder]  << "   "
			<< xl2_sin[iKT][Morder]  << "   " << t2_sin[iKT][Morder] << endl;
	}
}

	output_SvarscoeffsCOS.close();
	output_SvarscoeffsSIN.close();

	return;
}

void doHBT::Output_AVG_results()
{
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << runfolder << "/avgHBTradii.dat";
	filename_stream_HBT << global_runfolder << "/avgHBTradii.dat";
	ofstream outputHBT(filename_stream_HBT.str().c_str(), ios::app);
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << runfolder << "/avgHBTradii_cfs.dat";
	filename_stream_HBTcfs << global_runfolder << "/avgHBTradii_cfs.dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str(), ios::app);
	ostringstream filename_stream_S;
	//filename_stream_S << runfolder << "/avgSourcefunction_variances.dat";
	filename_stream_S << global_runfolder << "/avgSourcefunction_variances.dat";
	ofstream output_Svars(filename_stream_S.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		outputHBTcoeffs << "avg  " << K_T[iKT] << "  " << Morder
			<< "  " << avgR2_side_C[iKT][Morder] << "   " << avgR2_side_S[iKT][Morder] << "  " << avgR2_out_C[iKT][Morder] << "  " << avgR2_out_S[iKT][Morder]
			<< "  " << avgR2_outside_C[iKT][Morder] << "   " << avgR2_outside_S[iKT][Morder] << "  " << avgR2_long_C[iKT][Morder] << "  " << avgR2_long_S[iKT][Morder]
			<< "  " << avgR2_sidelong_C[iKT][Morder] << "   " << avgR2_sidelong_S[iKT][Morder] << "  " << avgR2_outlong_C[iKT][Morder] << "  " << avgR2_outlong_S[iKT][Morder] << endl;
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		outputHBT << "avg  " << K_T[iKT] << "  " << K_phi[iKphi]
			<< "  " << avgR2_side[iKT][iKphi] << "  " << avgR2_out[iKT][iKphi]
			<< "  " << avgR2_outside[iKT][iKphi] << "  " << avgR2_long[iKT][iKphi]
			<< "  " << avgR2_sidelong[iKT][iKphi] << "  " << avgR2_outlong[iKT][iKphi] << endl;
         	output_Svars << scientific << setprecision(8) << setw(15) 
			<< K_T[iKT] << "   " << K_phi[iKphi] << "   " << avgS_func[iKT][iKphi] << "   "
			<< avgxs_S[iKT][iKphi]  << "   " <<avgxo_S[iKT][iKphi]  << "   " << avgxl_S[iKT][iKphi]  << "   " << avgt_S[iKT][iKphi]  << "   "
			<< avgxs_t_S[iKT][iKphi]  << "   " << avgxo_t_S[iKT][iKphi]  << "   " << avgxl_t_S[iKT][iKphi]  << "   " << avgxo_xs_S[iKT][iKphi]  << "   "
			<< avgxl_xs_S[iKT][iKphi]  << "   " << avgxo_xl_S[iKT][iKphi]  << "   " << avgxs2_S[iKT][iKphi]  << "   " << avgxo2_S[iKT][iKphi]  << "   "
			<< avgxl2_S[iKT][iKphi]  << "   " << avgt2_S[iKT][iKphi] << endl;
	}
}

	outputHBT.close();
	outputHBTcoeffs.close();
	output_Svars.close();

	return;
}

void doHBT::Output_CAVG_results()
{
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << runfolder << "/CavgHBTradii.dat";
	filename_stream_HBT << global_runfolder << "/CavgHBTradii.dat";
	ofstream outputHBT(filename_stream_HBT.str().c_str(), ios::app);
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << runfolder << "/CavgHBTradii_cfs.dat";
	filename_stream_HBTcfs << global_runfolder << "/CavgHBTradii_cfs.dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		outputHBTcoeffs << "Cavg  " << K_T[iKT] << "  " << Morder
			<< "  " << CavgR2_side_C[iKT][Morder] << "   " << CavgR2_side_S[iKT][Morder] << "  " << CavgR2_out_C[iKT][Morder] << "  " << CavgR2_out_S[iKT][Morder]
			<< "  " << CavgR2_outside_C[iKT][Morder] << "   " << CavgR2_outside_S[iKT][Morder] << "  " << CavgR2_long_C[iKT][Morder] << "  " << CavgR2_long_S[iKT][Morder]
			<< "  " << CavgR2_sidelong_C[iKT][Morder] << "   " << CavgR2_sidelong_S[iKT][Morder] << "  " << CavgR2_outlong_C[iKT][Morder] << "  " << CavgR2_outlong_S[iKT][Morder] << endl;
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		outputHBT << "Cavg  " << K_T[iKT] << "  " << K_phi[iKphi]
			<< "  " << CavgR2_side[iKT][iKphi] << "  " << CavgR2_out[iKT][iKphi]
			<< "  " << CavgR2_outside[iKT][iKphi] << "  " << CavgR2_long[iKT][iKphi]
			<< "  " << CavgR2_sidelong[iKT][iKphi] << "  " << CavgR2_outlong[iKT][iKphi] << endl;
	}
}

outputHBT.close();
outputHBTcoeffs.close();
	return;
}

void doHBT::Readin_results(int folderindex)
{
double dummy;
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << path << folderindex << "/HBTradii_ev" << folderindex << ".dat";
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << path << folderindex << "/HBTradii_cfs_ev" << folderindex << ".dat";
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << ".dat";
	ifstream inputHBTcoeffs(filename_stream_HBTcfs.str().c_str());
	ostringstream filename_stream_S;
	//filename_stream_S << path << folderindex << "/Sourcefunction_variances.dat";
	filename_stream_S << global_path << "/Sourcefunction_variances.dat";
	ifstream input_Svars(filename_stream_S.str().c_str());
//cout << "Reading in from files " << filename_stream_HBT.str() << ", " << filename_stream_HBTcfs.str() << " and " << filename_stream_S.str() << endl;

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		inputHBTcoeffs >> dummy;
		inputHBTcoeffs >> dummy;
		inputHBTcoeffs >> dummy;
		inputHBTcoeffs >> R2_side_C[iKT][Morder];
		inputHBTcoeffs >> R2_side_S[iKT][Morder];
		inputHBTcoeffs >> R2_out_C[iKT][Morder];
		inputHBTcoeffs >> R2_out_S[iKT][Morder];
		inputHBTcoeffs >> R2_outside_C[iKT][Morder];
		inputHBTcoeffs >> R2_outside_S[iKT][Morder];
		inputHBTcoeffs >> R2_long_C[iKT][Morder];
		inputHBTcoeffs >> R2_long_S[iKT][Morder];
		inputHBTcoeffs >> R2_sidelong_C[iKT][Morder];
		inputHBTcoeffs >> R2_sidelong_S[iKT][Morder];
		inputHBTcoeffs >> R2_outlong_C[iKT][Morder];
		inputHBTcoeffs >> R2_outlong_S[iKT][Morder];
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		inputHBT >> dummy;
		inputHBT >> dummy;
        	inputHBT >> dummy;
		inputHBT >> R2_side[iKT][iKphi];
		inputHBT >> R2_out[iKT][iKphi];
		inputHBT >> R2_outside[iKT][iKphi];
		inputHBT >> R2_long[iKT][iKphi];
		inputHBT >> R2_sidelong[iKT][iKphi];
		inputHBT >> R2_outlong[iKT][iKphi];
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
         	input_Svars >> dummy;
        	input_Svars >> dummy;
        	input_Svars >> S_func[iKT][iKphi];
        	input_Svars >> xs_S[iKT][iKphi];
        	input_Svars >> xo_S[iKT][iKphi];
        	input_Svars >> xl_S[iKT][iKphi];
        	input_Svars >> t_S[iKT][iKphi];
        	input_Svars >> xs_t_S[iKT][iKphi];
        	input_Svars >> xo_t_S[iKT][iKphi];
        	input_Svars >> xl_t_S[iKT][iKphi];
        	input_Svars >> xo_xs_S[iKT][iKphi];
        	input_Svars >> xl_xs_S[iKT][iKphi];
        	input_Svars >> xo_xl_S[iKT][iKphi];
        	input_Svars >> xs2_S[iKT][iKphi];
        	input_Svars >> xo2_S[iKT][iKphi];
        	input_Svars >> xl2_S[iKT][iKphi];
        	input_Svars >> t2_S[iKT][iKphi];
	}
}

	inputHBT.close();
	inputHBTcoeffs.close();
	input_Svars.close();

	return;
}

void doHBT::Readin_HBTev_results_only(int folderindex)
{
double dummy;
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << path << folderindex << "/HBTradii_ev" << folderindex << ".dat";
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		inputHBT >> dummy;
		inputHBT >> dummy;
        	inputHBT >> dummy;
		inputHBT >> R2_side[iKT][iKphi];
		inputHBT >> R2_out[iKT][iKphi];
		inputHBT >> R2_outside[iKT][iKphi];
		inputHBT >> R2_long[iKT][iKphi];
		inputHBT >> R2_sidelong[iKT][iKphi];
		inputHBT >> R2_outlong[iKT][iKphi];
	}
}

	inputHBT.close();

	return;
}

void doHBT::Output_HBTcfsev_results_only(int folderindex)
{
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << path << folderindex << "/HBTradii_cfs_ev" << folderindex << ".dat";
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{Readin_AVG_results();
	for(int Morder=0; Morder<n_order; Morder++)
	{
		outputHBTcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
			<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
			<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
			<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
	}
}

	outputHBTcoeffs.close();

	return;
}

void doHBT::Readin_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << ".dat";
	ifstream inputplanepsi(filename_stream_planepsi.str().c_str());

	inputplanepsi >> global_plane_psi;

	inputplanepsi.close();

	return;
}

void doHBT::Output_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << ".dat";
	ofstream outputplanepsi(filename_stream_planepsi.str().c_str());

	outputplanepsi << global_plane_psi << endl;

	outputplanepsi.close();

	return;
}

void doHBT::Output_ev_plane_psis(int folderindex)
{
	ostringstream filename_stream_planepsis;
	//filename_stream_planepsis << path << folderindex << "/plane_psis_ev" << folderindex << ".dat";
	filename_stream_planepsis << global_path << "/plane_psis_ev" << folderindex << ".dat";
	ofstream outputplanepsis(filename_stream_planepsis.str().c_str());

	for (int i = 0; i < n_order; i++)
		outputplanepsis << i << "   " << plane_angle[i] << endl;

	outputplanepsis.close();

	return;
}

void doHBT::Output_ev_anisotropic_flows(int folderindex)
{
	ostringstream filename_stream_anisotropicflows;
	//filename_stream_anisotropicflows << path << folderindex << "/anisotropic_flows_ev" << folderindex << ".dat";
	filename_stream_anisotropicflows << global_path << "/anisotropic_flows_ev" << folderindex << ".dat";
	ofstream outputanisotropicflows(filename_stream_anisotropicflows.str().c_str());

	for (int i = 0; i < n_order; i++)
		outputanisotropicflows << i << "   " << anisotropic_flows[i] << endl;

	outputanisotropicflows.close();

	return;
}

void doHBT::Output_ev_anisotropic_flows_pTdiff(int folderindex)
{
	ostringstream filename_stream_anisotropicflowspTdiff;
	//filename_stream_anisotropicflowspTdiff << path << folderindex << "/anisotropic_flows_pTdiff_ev" << folderindex << ".dat";
	filename_stream_anisotropicflowspTdiff << global_path << "/anisotropic_flows_pTdiff_ev" << folderindex << ".dat";
	ofstream outputanisotropicflowspTdiff(filename_stream_anisotropicflowspTdiff.str().c_str());

	for (int i = 0; i < n_SP_pT; i++)
	for (int j = 0; j < n_order; j++)
		outputanisotropicflowspTdiff << j << "   " << SP_pT[i] << "   " << anisotropic_flows_pTdiff[i][j] << "   " << anisotropic_flows_pTdiff_psin[i][j] << endl;

	outputanisotropicflowspTdiff.close();

	return;
}

void doHBT::Output_ev_mean_pT(int folderindex)
{
	ostringstream omeanpT_stream;
	//omeanpT_stream << path << folderindex << "/meanpT_ev" << folderindex << ".dat";
	omeanpT_stream << global_path << "/meanpT_ev" << folderindex << ".dat";
	ofstream omeanpT;
	omeanpT.open(omeanpT_stream.str().c_str());

	omeanpT << mean_pT << endl;

	omeanpT.close();

	return;
}

void doHBT::Readin_AVG_results()
{
	ostringstream filename_stream_HBT;
	//filename_stream_HBT << runfolder << "/avgHBTradii.dat";
	filename_stream_HBT << global_runfolder << "/avgHBTradii.dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	//filename_stream_HBTcfs << runfolder << "/avgHBTradii_cfs.dat";
	filename_stream_HBTcfs << global_runfolder << "/avgHBTradii_cfs.dat";
	ifstream inputHBTcoeffs(filename_stream_HBTcfs.str().c_str());
	ostringstream filename_stream_S;
	//filename_stream_S << runfolder << "/avgSourcefunction_variances.dat";
	filename_stream_S << global_runfolder << "/avgSourcefunction_variances.dat";
	ifstream input_Svars(filename_stream_S.str().c_str());
//cout << "Reading in from files " << filename_stream_HBT.str() << " and " << filename_stream_HBTcfs.str() << endl;
double dummy;
string dummystring = "";

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		inputHBTcoeffs >> dummystring;
		inputHBTcoeffs >> dummy;
		inputHBTcoeffs >> dummy;
		inputHBTcoeffs >> avgR2_side_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_side_S[iKT][Morder];
		inputHBTcoeffs >> avgR2_out_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_out_S[iKT][Morder];
		inputHBTcoeffs >> avgR2_outside_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_outside_S[iKT][Morder];
		inputHBTcoeffs >> avgR2_long_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_long_S[iKT][Morder];
		inputHBTcoeffs >> avgR2_sidelong_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_sidelong_S[iKT][Morder];
		inputHBTcoeffs >> avgR2_outlong_C[iKT][Morder];
		inputHBTcoeffs >> avgR2_outlong_S[iKT][Morder];
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		inputHBT >> dummystring;
		inputHBT >> dummy;
		inputHBT >> dummy;
		inputHBT >> avgR2_side[iKT][iKphi];
		inputHBT >> avgR2_out[iKT][iKphi];
		inputHBT >> avgR2_outside[iKT][iKphi];
		inputHBT >> avgR2_long[iKT][iKphi];
		inputHBT >> avgR2_sidelong[iKT][iKphi];
		inputHBT >> avgR2_outlong[iKT][iKphi];

		input_Svars >> dummy;
		input_Svars >> dummy;
		input_Svars >> avgS_func[iKT][iKphi];
		input_Svars >> avgxs_S[iKT][iKphi];
		input_Svars >> avgxo_S[iKT][iKphi];
		input_Svars >> avgxl_S[iKT][iKphi];
		input_Svars >> avgt_S[iKT][iKphi];
		input_Svars >> avgxs_t_S[iKT][iKphi];
		input_Svars >> avgxo_t_S[iKT][iKphi];
		input_Svars >> avgxl_t_S[iKT][iKphi];
		input_Svars >> avgxo_xs_S[iKT][iKphi];
		input_Svars >> avgxl_xs_S[iKT][iKphi];
		input_Svars >> avgxo_xl_S[iKT][iKphi];
		input_Svars >> avgxs2_S[iKT][iKphi];
		input_Svars >> avgxo2_S[iKT][iKphi];
		input_Svars >> avgxl2_S[iKT][iKphi];
		input_Svars >> avgt2_S[iKT][iKphi];
	}
}

	inputHBT.close();
	inputHBTcoeffs.close();
	input_Svars.close();

	return;
}

void doHBT::Output_Correlationfunction_1D(int iKT, int iKphi, int folderindex)
{
   ostringstream oCorrelfun_1D_stream;
   //oCorrelfun_1D_stream << path << folderindex << "/correlfunct1D" << "_" << particle_name << "_kt_" << K_T[iKT] << "_phi_" << K_phi[iKphi] << ".dat";
   oCorrelfun_1D_stream << global_path << "/correlfunct1D" << "_" << particle_name << "_kt_" << K_T[iKT] << "_phi_" << K_phi[iKphi] << ".dat";
   ofstream oCorrelfun_1D;
   oCorrelfun_1D.open(oCorrelfun_1D_stream.str().c_str(), ios::app);
   for(int i=0; i < qnpts; i++)
     oCorrelfun_1D << scientific << setprecision(7) << setw(15)
                   << q_out[i] << "  " << Correl_1D_out[i] << "  " << Correl_1D_out_err[i] << "  "
                   << q_side[i] << "  " << Correl_1D_side[i] << "  " << Correl_1D_side_err[i] << "  "
                   << q_long[i] << "  " << Correl_1D_long[i] << "  " << Correl_1D_long_err[i] 
                   << endl;
   return;
}

void doHBT::Output_Correlationfunction_3D(int iKT, int iKphi, int folderindex)
{
   ostringstream oCorrelfun_3D_stream;
   //oCorrelfun_3D_stream << path << folderindex << "/correlfunct3D" << "_" << particle_name << "_kt_" << K_T[iKT] << "_phi_" << K_phi[iKphi] << ".dat";
   oCorrelfun_3D_stream << global_path << "/correlfunct3D" << "_" << particle_name << "_kt_" << K_T[iKT] << "_phi_" << K_phi[iKphi] << ".dat";
   ofstream oCorrelfun_3D;
   oCorrelfun_3D.open(oCorrelfun_3D_stream.str().c_str(), ios::app);
   for(int i=0; i < qnpts; i++)
      for(int j=0; j < qnpts; j++)
         for(int k=0; k < qnpts; k++)
              oCorrelfun_3D << scientific << setprecision(7) << setw(15)
                            << q_out[i] << "  " << q_side[j] << "  " 
                            << q_long[k] << "  " << Correl_3D[i][j][k] << "  "
                            << Correl_3D_err[i][j][k] << endl;
   return;
}

void doHBT::Output_GF_results(int folderindex)
{
	ostringstream filename_stream_HBTGF;
	//filename_stream_HBTGF << path << folderindex << "/HBTradii_GF_ev" << folderindex << ".dat";
	filename_stream_HBTGF << global_path << "/HBTradii_GF_ev" << folderindex << ".dat";
	ofstream outputHBTGF(filename_stream_HBTGF.str().c_str(), ios::app);
	ostringstream filename_stream_HBTGFcfs;
	//filename_stream_HBTGFcfs << path << folderindex << "/HBTradii_GF_cfs_ev" << folderindex << ".dat";
	filename_stream_HBTGFcfs << global_path << "/HBTradii_GF_cfs_ev" << folderindex << ".dat";
	ofstream outputHBTGFcoeffs(filename_stream_HBTGFcfs.str().c_str(), ios::app);

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int Morder=0; Morder<n_order; Morder++)
	{
		outputHBTGFcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
			<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
			<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
			<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
	}
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		outputHBTGF << folderindex << "  " << K_T[iKT] << "  " << K_phi[iKphi]
			<< "  " << R2_side[iKT][iKphi] << "  " << R2_out[iKT][iKphi]
			<< "  " << R2_outside[iKT][iKphi] << "  " << R2_long[iKT][iKphi]
			<< "  " << R2_sidelong[iKT][iKphi] << "  " << R2_outlong[iKT][iKphi] << endl;
	}
}

	outputHBTGF.close();
	outputHBTGFcoeffs.close();

	return;
}

void doHBT::Output_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpTdphi;
	//filename_stream_dN_dypTdpTdphi << path << folderindex << "/dN_dypTdpTdphi_ev" << folderindex << ".dat";
	filename_stream_dN_dypTdpTdphi << global_path << "/dN_dypTdpTdphi_ev" << folderindex << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpTdphi << SP_pT[ipt] << "   " << SP_pphi[iphi] << "   " << dN_dypTdpTdphi[ipt][iphi] << endl;

	output_dN_dypTdpTdphi.close();

	return;
}

void doHBT::Output_dN_dypTdpT(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpT;
	//filename_stream_dN_dypTdpT << path << folderindex << "/dN_dypTdpT_ev" << folderindex << ".dat";
	filename_stream_dN_dypTdpT << global_path << "/dN_dypTdpT_ev" << folderindex << ".dat";
	ofstream output_dN_dypTdpT(filename_stream_dN_dypTdpT.str().c_str());

	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpT << SP_pT[ipt] << "   " << dN_dypTdpT[ipt] << endl;

	output_dN_dypTdpT.close();

	return;
}

void doHBT::Output_avgEmission_Function_on_FOsurface(int folderindex)
{
	ostringstream filename_stream_S;

	//filename_stream_S << path << folderindex << "/averaged_S_on_FOsurface_ev" << folderindex << ".dat";
	filename_stream_S << global_path << "/averaged_S_on_FOsurface_ev" << folderindex << ".dat";
	*global_out_stream_ptr << "Output averaged S to " << filename_stream_S.str() << endl;

	ofstream output_S(filename_stream_S.str().c_str());

	int idx = 0;

	for (int iKT = 0; iKT < n_localp_T; iKT++)
	{
		for (int iFOcell = 0; iFOcell < FO_length; iFOcell++)
		{
			output_S << scientific << setprecision(8) << setw(15)
				<< K_T[iKT] << "   "
				<< (*avgFOsurf_ptr)[idx].tau << "   "
				<< (*avgFOsurf_ptr)[idx].x << "   "
				<< (*avgFOsurf_ptr)[idx].y << "   "
				<< (*avgFOsurf_ptr)[idx].data << endl;
			idx++;
		}
	}

	output_S.close();

	return;
}

void doHBT::Output_Emission_Function(int iKT, int iKphi, int folderindex)
{
int coords = 2;
//coords:	0 - output in (x,y,z,t,data) format
//		1 - output in (r,phi,eta,tau,data) format
//		2 - output in (x,y,eta,tau,data) format

	ostringstream filename_stream_S;

//	filename_stream_S << path << folderindex << "/Emissionfunction_S_spacetime_kT_" 
//			<< fixed << setprecision(2) << setw(3) << K_T[iKT] << "_kphi_"
//			<< fixed << setprecision(5) << setw(6) << K_phi[iKphi] << ".dat";
	filename_stream_S << global_path << "/Emissionfunction_S_spacetime_kT_" 
			<< fixed << setprecision(2) << setw(3) << K_T[iKT] << "_kphi_"
			<< fixed << setprecision(5) << setw(6) << K_phi[iKphi] << ".dat";
	cout << "Output Emissionfunction to " << filename_stream_S.str() << endl;

	ofstream output_S(filename_stream_S.str().c_str());

cerr << "Emissionfunction_length in Output_Emission_Function is " << Emissionfunction_length << endl;

	for(int i=0; i<Emissionfunction_length; i++)
	{
		if (coords == 2) {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].eta << "   "
				<< (*Emissionfunction_ptr)[i].tau << "   "
				<< (*Emissionfunction_ptr)[i].x << "   "
				<< (*Emissionfunction_ptr)[i].y << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		else if (coords == 1) {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].r << "   "
				<< (*Emissionfunction_ptr)[i].phi << "   "
				<< (*Emissionfunction_ptr)[i].eta << "   "
				<< (*Emissionfunction_ptr)[i].tau << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		else {
		output_S << scientific << setprecision(8) << setw(15)
				<< (*Emissionfunction_ptr)[i].x << "   "
				<< (*Emissionfunction_ptr)[i].y << "   "
				<< (*Emissionfunction_ptr)[i].z << "   "
				<< (*Emissionfunction_ptr)[i].t << "   "
				<< (*Emissionfunction_ptr)[i].data << "   ";
		}
		output_S << endl;
	}

	output_S.close();

	return;
}

//print output to output filestream, one line at a time
void doHBT::Set_ofstream(ofstream& myout)
{
	global_out_stream_ptr = &myout;

	return;
}

/*void doHBT::Set_folderindex(int folderindex)
{
	global_folderindex = folderindex;

	return;
}*/

//print output to output filestream, one line at a time
void doHBT::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void doHBT::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

/*void doHBT::Set_n_order(int user_def_n_order)
{
	n_order = user_def_n_order;

	return;
}*/

//End of file
