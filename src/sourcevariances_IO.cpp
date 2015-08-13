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

#include "sourcevariances.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

void SourceVariances::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str());
	ostringstream filename_stream_S;
	filename_stream_S << global_path << "/Sourcefunction_variances" << no_df_stem << ".dat";
	ofstream output_Svars(filename_stream_S.str().c_str());

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
	         	//output_Svars << scientific << setprecision(8) << setw(15) 
	         	output_Svars << setprecision(8) << setw(15) 
				<< K_T[iKT] << "   " << K_phi[iKphi] << "   " << S_func[iKT][iKphi] << "   "
				<< xs_S[iKT][iKphi] << "   " << xo_S[iKT][iKphi] << "   " << xl_S[iKT][iKphi] << "   "
				<< t_S[iKT][iKphi]  << "   " << xs_t_S[iKT][iKphi] << "   "
				<< xo_t_S[iKT][iKphi] << "   " << xl_t_S[iKT][iKphi] << "   "
				<< xo_xs_S[iKT][iKphi] << "   " << xl_xs_S[iKT][iKphi] << "   "
				<< xo_xl_S[iKT][iKphi] << "   " << xs2_S[iKT][iKphi] << "   " << xo2_S[iKT][iKphi] << "   "
				<< xl2_S[iKT][iKphi] << "   " << t2_S[iKT][iKphi] << endl;
		}
	}

	outputHBT.close();
	output_Svars.close();

	return;
}

void SourceVariances::Readin_results(int folderindex)
{
double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_S;
	filename_stream_S << global_path << "/Sourcefunction_variances" << no_df_stem << ".dat";
	ifstream input_Svars(filename_stream_S.str().c_str());

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
	input_Svars.close();

	return;
}

/*void SourceVariances::Readin_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputplanepsi(filename_stream_planepsi.str().c_str());

	inputplanepsi >> global_plane_psi;

	inputplanepsi.close();

	return;
}

void SourceVariances::Output_ev_plane_psi(int folderindex)
{
	ostringstream filename_stream_planepsi;
	//filename_stream_planepsi << path << folderindex << "/plane_psi_ev" << folderindex << ".dat";
	filename_stream_planepsi << global_path << "/plane_psi_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputplanepsi(filename_stream_planepsi.str().c_str());

	outputplanepsi << global_plane_psi << endl;

	outputplanepsi.close();

	return;
}

void SourceVariances::Output_ev_plane_psis(int folderindex)
{
	ostringstream filename_stream_planepsis;
	//filename_stream_planepsis << path << folderindex << "/plane_psis_ev" << folderindex << ".dat";
	filename_stream_planepsis << global_path << "/plane_psis_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputplanepsis(filename_stream_planepsis.str().c_str());

	for (int i = 0; i < n_order; i++)
		outputplanepsis << i << "   " << plane_angle[i] << endl;

	outputplanepsis.close();

	return;
}*/

void SourceVariances::Output_SVdN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << global_path << "/SV_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpTdphi << SP_pT[ipt] << "   " << SP_pphi[iphi] << "   " << SV_dN_dypTdpTdphi[ipt][iphi] << endl;

	output_dN_dypTdpTdphi.close();

	return;
}

void SourceVariances::Output_SVdN_dypTdpT(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpT;
	filename_stream_dN_dypTdpT << global_path << "/SV_dN_dypTdpT_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpT(filename_stream_dN_dypTdpT.str().c_str());

	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpT << SP_pT[ipt] << "   " << SV_dN_dypTdpT[ipt] << endl;

	output_dN_dypTdpT.close();

	return;
}

void SourceVariances::Output_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << global_path << "/dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int iphi=0; iphi<n_SP_pphi; iphi++)
	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpTdphi << SP_pT[ipt] << "   " << SP_pphi[iphi] << "   " << dN_dypTdpTdphi[ipt][iphi] << endl;

	output_dN_dypTdpTdphi.close();

	return;
}

void SourceVariances::Output_dN_dypTdpTdphi_grid(int folderindex, int dc_idx)
{//assuming polar grid for interpolation for now
	ostringstream filename_stream_dN_dypTdpTdphi;
	filename_stream_dN_dypTdpTdphi << global_path << "/dN_dypTdpTdphi_polar_grid_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpTdphi(filename_stream_dN_dypTdpTdphi.str().c_str());

	for(int imom=0; imom<n_weighting_functions; imom++)
	{
		for(int iphi=0; iphi<n_interp2_pphi_pts; iphi++)
		for(int ipt=0; ipt<n_interp2_pT_pts; ipt++)
			output_dN_dypTdpTdphi << imom << "   " << SPinterp2_pT[ipt] << "   " << SPinterp2_pphi[iphi]
						<< "   " << dN_dypTdpTdphi_moments[dc_idx][imom][ipt][iphi] << endl;
		output_dN_dypTdpTdphi << endl;
	}

	output_dN_dypTdpTdphi.close();

	return;
}

void SourceVariances::Output_dN_dypTdpT(int folderindex)
{
	ostringstream filename_stream_dN_dypTdpT;
	filename_stream_dN_dypTdpT << global_path << "/dN_dypTdpT_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_dN_dypTdpT(filename_stream_dN_dypTdpT.str().c_str());

	for(int ipt=0; ipt<n_SP_pT; ipt++)
		output_dN_dypTdpT << SP_pT[ipt] << "   " << dN_dypTdpT[ipt] << endl;

	output_dN_dypTdpT.close();

	return;
}

void SourceVariances::Output_all_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_all_dN_dypTdpTdphi;
	filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());

	for(int ii = 0; ii < Nparticle; ii++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	{
		for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
			output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ii][0][ipt][iphi] << "   ";
		output_all_dN_dypTdpTdphi << endl;
	}

	output_all_dN_dypTdpTdphi.close();

	return;
}

/*void SourceVariances::Output_all_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_all_dN_dypTdpTdphi;
	filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());

	//for(int idc = 1; idc <= n_decay_channels; idc++)
	//{
		//for(int iKT = 0; iKT < n_localp_T; iKT++)
		//for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		//	output_all_dN_dypTdpTdphi << K_T[iKT] << "   " << K_phi[iKphi] << "   " << integrated_spacetime_moments[idc-1][0][iKT][iKphi] << endl;
		//output_all_dN_dypTdpTdphi << endl;
	//}

	//for(int ir = 0; ir <= n_resonance; ir++)
	for(int ii = 0; ii < Nparticle; ii++)
	{
		//thermal only for timebeing...
		std::vector<int>::iterator it = find (chosen_resonances.begin(), chosen_resonances.end(), ii);
		//if (VERBOSE > 0) *global_out_stream_ptr << "particle_id = " << particle_id << endl;
		if (ii == particle_id)
		{
			for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
			{
				for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
					output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[0][0][ipt][iphi] << "   ";
				output_all_dN_dypTdpTdphi << endl;
			}

		}
		else if (it != chosen_resonances.end())	// if particle was one of chosen resonances
		{
			int ir = it - chosen_resonances.begin() + 1;
			//if (VERBOSE > 0) *global_out_stream_ptr << "ir = " << ir << endl;

			for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
			{
				for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
					output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ir][0][ipt][iphi] << "   ";
				output_all_dN_dypTdpTdphi << endl;
			}
		}
		else
		{
			for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
			{
				for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
					output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << 0.0 << "   ";
				output_all_dN_dypTdpTdphi << endl;
			}
		}

	}

	output_all_dN_dypTdpTdphi.close();

	return;
}*/

void SourceVariances::Read_in_all_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_all_dN_dypTdpTdphi;
	filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ifstream input_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());

	int local_filelength = get_filelength(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	int local_filewidth = get_filewidth(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	if (VERBOSE > 0) *global_out_stream_ptr << "Read_in_all_dN_dypTdpTdphi(): nrows = " << local_filelength << " and ncols = " << local_filewidth << endl;
	if ((Nparticle * n_interp2_pphi_pts != local_filelength) || (n_interp2_pT_pts != local_filewidth))
	{
		cerr << "Read_in_all_dN_dypTdpTdphi(): Mismatch in dimensions!" << endl;
		exit(1);
	}

	for(int ii = 0; ii < Nparticle; ii++)
	for(int iphi = 0; iphi < n_interp2_pphi_pts; iphi++)
	for(int ipt = 0; ipt < n_interp2_pT_pts; ipt++)
		input_all_dN_dypTdpTdphi >> dN_dypTdpTdphi_moments[ii][0][ipt][iphi];

	input_all_dN_dypTdpTdphi.close();

	return;
}

//End of file
