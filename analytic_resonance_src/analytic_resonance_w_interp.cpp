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

#include "analytic_resonance.h"
#include "gauss.h"
//#include "CMDreader.h"

bool output_to_screen = true;
bool output_to_file = true;

size_t n;
string input_filename;
string input_filepath;

double m, Mres, Gamma, m2, m3, br, sign, Qfunc;
int n_body;

int main(int argc, char *argv[])
{
	/*int success = Read_in_parameters_from_CMD(argc, argv);
	if (success > 0)
	{
		cerr << "Problem reading in command-line parameters!" << endl;
		exit(1);
	}
	if (output_to_screen)
	{
		cout << "Initial parameters are:" << endl
			<< "eps_2_bar = " << eps_2_bar << endl
			<< "v_2_bar = " << v_2_bar << endl
			<< "psi_2_bar = " << psi_2_bar << endl
			<< "eta_f = " << eta_f << endl
			<< "Rad = " << Rad << endl
			<< "T0 = " << T0 << endl;
	}*/
	
/**************************************************************************/
/*************************MAIN PART OF THE PROGRAM*************************/
/**************************************************************************/

	//set up stuff for integrations
	Create_integrations_points_and_weights();

	//get spectra
	//Compute_direct_pion_spectra_OLD();
	for (int iKT = 0; iKT < n_K_T; iKT++)
		Direct_contributions_to_pion_spectra(K_T[iKT], 0.0, 0.0);
	Compute_direct_resonance_spectra();
	
	for (int iKT = 0; iKT < n_K_T; iKT++)
		Resonance_decay_contributions_to_pion_spectra(K_T[iKT], 0.0, 0.0, 1);

/**************************************************************************/
/***********************END MAIN PART OF THE PROGRAM***********************/
/**************************************************************************/
	return (0);
}

/**************************************************************************/

double S_thermal(double r, double phi, double eta, double tau, double PT, double Y, double Phi, int reso_idx /*= 0*/)
{
	double M = m_pion;
	if (reso_idx != 0)
		M = resonances.resonance_mass[reso_idx-1];
	double MT = sqrt(M*M + PT*PT);
	double rt = r / Rad;
	double He1 = 0.5*rt*rt, He2 = 0.5*(eta-eta0)*(eta-eta0)/(Deleta*Deleta), He3 = 0.5*(tau-tau0)*(tau-tau0)/(Deltau*Deltau);
	//He2 = 0.0;
	double H = exp(-He1 - He2 - He3) / (M_PI*Deltau);
	double eta_t = eta_f * rt;
	double prefactor = MT*cosh(Y-eta)/(8.*M_PI*M_PI*M_PI*hbarC3);	//2Jr+1 == 1 for pions, lumped into br for all other resonances...
	double term1 = (MT/Tdec)*cosh(Y-eta)*cosh(eta_t), term2 = (PT/Tdec)*sinh(eta_t)*cos(phi-Phi);
	return (r * tau * prefactor * H * exp(-term1 + term2));	//r*tau factor is jacobian factor
	//return (r*tau);
}

double tauintegrated_S_thermal(double r, double phi, double eta, double PT, double Y, double Phi, int reso_idx /*= 0*/)
{
	double M = m_pion;
	if (reso_idx != 0)
		M = resonances.resonance_mass[reso_idx-1];
	double MT = sqrt(M*M + PT*PT);
	double rt = r / Rad;
	double He1 = 0.5*rt*rt, He2 = 0.5*(eta-eta0)*(eta-eta0)/(Deleta*Deleta);
	//He2 = 0.0;
	double H = exp(-He1 - He2) / (M_PI*Deltau);
	double eta_t = eta_f * rt;
	double prefactor = MT*cosh(Y-eta)/(8.*M_PI*M_PI*M_PI*hbarC3);	//2Jr+1 == 1 for pions, lumped into br for all other resonances...
	double term1 = (MT/Tdec)*cosh(Y-eta)*cosh(eta_t), term2 = (PT/Tdec)*sinh(eta_t)*cos(phi-Phi);
	return (r * prefactor * H * exp(-term1 + term2));	//r factor is jacobian factor
}

void Create_integrations_points_and_weights()
{
	r_pts = new double [nrpts];
	r_wts = new double [nrpts];
	phi_pts = new double [nphipts];
	phi_wts = new double [nphipts];
	eta_pts = new double [netapts];
	eta_wts = new double [netapts];
	tau_pts = new double [ntaupts];
	tau_wts = new double [ntaupts];
	gauss (nrpts, 0, r_min, r_max, r_pts, r_wts);
	gauss (nphipts, 0, phi_min, phi_max, phi_pts, phi_wts);
	gauss (netapts, 0, eta_min, eta_max, eta_pts, eta_wts);
	gauss (ntaupts, 0, tau_min, tau_max, tau_pts, tau_wts);
	
	tau_integrated_S_prefactor = Deltau*Deltau*exp(-0.5*tau0*tau0/(Deltau*Deltau))
					+ sqrt(0.5*M_PI)*Deltau*tau0*( 1. + erf( tau0/( sqrt(2.)*Deltau ) ) );
	tau_integrated_Stau_prefactor = Deltau*Deltau*tau0*exp(-0.5*tau0*tau0/(Deltau*Deltau))
					+ sqrt(0.5*M_PI)*Deltau*(Deltau*Deltau + tau0*tau0)*( 1. + erf( tau0/( sqrt(2.)*Deltau ) ) );
	tau_integrated_Stau2_prefactor = Deltau*Deltau*(2.*Deltau*Deltau + tau0*tau0)*exp(-0.5*tau0*tau0/(Deltau*Deltau))
					+ sqrt(0.5*M_PI)*Deltau*tau0*(3.*Deltau*Deltau + tau0*tau0)*( 1. + erf( tau0/( sqrt(2.)*Deltau ) ) );

	//tau_integrated_S_prefactor /= M_PI*Deltau;
	//tau_integrated_Stau_prefactor /= M_PI*Deltau;
	//tau_integrated_Stau2_prefactor /= M_PI*Deltau;
	
	//load resonances
	n_resonance = read_in_resonances(&resonances);
	for (int ires = 0; ires < n_resonance; ires++)
		cout << resonances.resonance_mass[ires] << endl;

	//set pair-momentum info
	K_T = new double [n_K_T];
	K_phi = new double [n_K_phi];
	K_phi_weights = new double [n_K_phi];
	double dK_T = (K_T_max - K_T_min)/((double)n_K_T - 1. + 1.e-100);
	for (int iKT = 0; iKT < n_K_T; iKT++)
		K_T[iKT] = K_T_min + iKT*dK_T;
	gauss (n_K_phi, 0, K_phi_min, K_phi_max, K_phi, K_phi_weights);
	
	//set points for interpolation grid
	PTpts = new double [nPTpts];
	double dPT = (PTmax - PTmin)/((double)nPTpts - 1. + 1.e-100);
	for (int ipt = 0; ipt < nPTpts; ipt++)
		PTpts[ipt] = PTmin + ipt*dPT;
	PPhipts = new double [nPPhipts];
	double dPPhi = (PPhimax - PPhimin)/((double)nPPhipts - 1. + 1.e-100);	//since 0 ~ 2*pi
	for (int ipphi = 0; ipphi < nPPhipts; ipphi++)
		PPhipts[ipphi] = PPhimin + ipphi*dPPhi;
	PYpts = new double [nPYpts];
	double dPY = (PYmax - PYmin)/((double)nPYpts - 1. + 1.e-100);
	for (int ipy = 0; ipy < nPYpts; ipy++)
		PYpts[ipy] = PYmin + ipy*dPY;
	
	//set up array for all resonance spectra at interpolation grid points
	resonance_spectra = new double *** [n_resonance];
	for (int ires = 0; ires < n_resonance; ires++)
	{
		resonance_spectra[ires] = new double ** [nPTpts];
		for (int ipt = 0; ipt < nPTpts; ipt++)
		{
			resonance_spectra[ires][ipt] = new double * [nPPhipts];
			for (int ipphi = 0; ipphi < nPPhipts; ipphi++)
			{
				resonance_spectra[ires][ipt][ipphi] = new double [nPYpts];
				for (int ipy = 0; ipy < nPYpts; ipy++)
					resonance_spectra[ires][ipt][ipphi][ipy] = 0.0;
			}
		}
	}
	
	return;
}

double Direct_contributions_to_pion_spectra(double pT, double y, double pphi)
{
	double sum = 0.0;
	//double m = m_pion;
	double symmetry_factor = 1.0;	//for eta integration
	if (ASSUME_ETA_SYMMETRIC)
		symmetry_factor = 2.0;
	//double mT = sqrt(m*m + pT*pT);
	for (int ir = 0; ir < nrpts; ir++)
	for (int iphi = 0; iphi < nphipts; iphi++)
	for (int ieta = 0; ieta < netapts; ieta++)
		sum += symmetry_factor * r_wts[ir] * phi_wts[iphi] * eta_wts[ieta]
			* tau_integrated_S_prefactor * tauintegrated_S_thermal(r_pts[ir], phi_pts[iphi], eta_pts[ieta], pT, y, pphi);
	cout << "dN_dypTdpTdpphi(pT = " << pT << ", pphi = " << pphi << ") = " << sum << endl;
	return sum;
}

void Compute_direct_resonance_spectra()
{
	double symmetry_factor = 1.0;	//for eta integration
	if (ASSUME_ETA_SYMMETRIC)
		symmetry_factor = 2.0;
	/*double sum = 0.0, KTloc = 0.0, Kphiloc = 0.0;
	double m = m_pion, Yloc = 0.0;
	for (int iKT = 0; iKT < n_K_T; iKT++)
	//for (int iKphi = 0; iKphi < n_K_phi; iKphi++)
	{
		sum = 0.0;
		KTloc = K_T[iKT];
		//Kphiloc = K_phi[iKphi];
		for (int ir = 0; ir < nrpts; ir++)
		for (int iphi = 0; iphi < nphipts; iphi++)
		for (int ieta = 0; ieta < netapts; ieta++)
		for (int itau = 0; itau < ntaupts; itau++)
			sum += symmetry_factor * r_wts[ir] * phi_wts[iphi] * eta_wts[ieta] * tau_wts[itau]
				* S_thermal(r_pts[ir], phi_pts[iphi], eta_pts[ieta], tau_pts[itau], KTloc, Yloc, Kphiloc);
		cout << "dN_dypTdpTdpphi(pT = " << KTloc << ", pphi = " << Kphiloc << ") = " << sum << endl;
	}*/
	for (int ires = 0; ires < 1; ires++)
	for (int ipt = 0; ipt < nPTpts; ipt++)
	for (int ipphi = 0; ipphi < nPPhipts; ipphi++)
	for (int ipy = 0; ipy < nPYpts; ipy++)
	{
	for (int ir = 0; ir < nrpts; ir++)
	for (int iphi = 0; iphi < nphipts; iphi++)
	for (int ieta = 0; ieta < netapts; ieta++)
		resonance_spectra[ires][ipt][ipphi][ipy] += symmetry_factor * tau_integrated_S_prefactor * r_wts[ir] * phi_wts[iphi] * eta_wts[ieta]
				* tauintegrated_S_thermal(r_pts[ir], phi_pts[iphi], eta_pts[ieta], PTpts[ipt], PYpts[ipy], PPhipts[ipphi], ires + 1);
	//cerr << PTpts[ipt] << "   " << PPhipts[ipphi] << "   " << PYpts[ipy] << "   " << resonance_spectra[ires][ipt][ipphi][ipy] << endl;
	}
}

void Compute_direct_pion_spectra_OLD()
{
	double sum = 0.0, KTloc = 0.0, Kphiloc = 0.0;
	double symmetry_factor = 1.0;	//for eta integration
	double Yloc = 0.0;
	for (int iKT = 0; iKT < n_K_T; iKT++)
	//for (int iKphi = 0; iKphi < n_K_phi; iKphi++)
	{
		sum = 0.0;
		KTloc = K_T[iKT];
		//Kphiloc = K_phi[iKphi];
		for (int ir = 0; ir < nrpts; ir++)
		for (int iphi = 0; iphi < nphipts; iphi++)
		for (int ieta = 0; ieta < netapts; ieta++)
		for (int itau = 0; itau < ntaupts; itau++)
			sum += symmetry_factor * r_wts[ir] * phi_wts[iphi] * eta_wts[ieta] * tau_wts[itau]
				* S_thermal(r_pts[ir], phi_pts[iphi], eta_pts[ieta], tau_pts[itau], KTloc, Yloc, Kphiloc);
		cout << setw(8) << setprecision(15) << "dN_dypTdpTdpphi(pT = " << KTloc << ", pphi = " << Kphiloc << ") = " << sum << endl;
	}
}

int read_in_resonances(resonance_info * resonances)
{
	int number_of_resonances = 0;
	string resonancefilename = "/home/plumberg.1/HBTPlumberg/EOS/temporary_resonance_data.dat";
	ifstream resonanceinput (resonancefilename.c_str());
	resonanceinput >> number_of_resonances;
	resonances->resonance_mass = new double [number_of_resonances];
	resonances->resonance_Gamma = new double [number_of_resonances];
	resonances->resonance_total_br = new double [number_of_resonances];
	resonances->resonance_mu = new double [number_of_resonances];
	resonances->resonance_gspin = new double [number_of_resonances];
	resonances->resonance_sign = new int [number_of_resonances];
	resonances->resonance_decay_masses = new double* [number_of_resonances];
	for (int ir=0; ir<number_of_resonances; ir++)
	{
		resonances->resonance_decay_masses[ir] = new double [2];
		resonances->resonance_decay_masses[ir][0] = 0.0;
		resonances->resonance_decay_masses[ir][1] = 0.0;
		resonances->resonance_mu[ir] = 0.0;
		resonances->resonance_gspin[ir] = 1.0;	//actual g's have been absorbed into definitions of br
		//resonances->resonance_sign[ir] = 1;	//not quite right
	}
	int row_index = 0;
	resonanceinput >> row_index;
	while (!resonanceinput.eof() && row_index != 0)
	{
		//note that we have to convert given table values to GeV
		resonanceinput >> resonances->resonance_mass[row_index-1];
		resonanceinput >> resonances->resonance_decay_masses[row_index-1][0];
		resonanceinput >> resonances->resonance_decay_masses[row_index-1][1];
		resonanceinput >> resonances->resonance_Gamma[row_index-1];
		resonanceinput >> resonances->resonance_total_br[row_index-1];
		resonanceinput >> resonances->resonance_sign[row_index-1];
		//if (DEBUG)
		//	cerr << "Made it through row_index = " << row_index << endl;
		resonanceinput >> row_index;
	}
	resonanceinput.close();
	for (int ir=0; ir<number_of_resonances; ir++)
	{
		resonances->resonance_mass[ir] *= MeVToGeV;
		resonances->resonance_decay_masses[ir][0] *= MeVToGeV;
		resonances->resonance_decay_masses[ir][1] *= MeVToGeV;
		resonances->resonance_Gamma[ir] *= MeVToGeV;
	}
	
	return (number_of_resonances);
}

void Set_resonance_integration_points(double smin, double smax)
{
	s_pts = new double [nspts];
	s_wts = new double [nspts];
	v_pts = new double [nvpts];
	v_wts = new double [nvpts];
	zeta_pts = new double [nzetapts];
	zeta_wts = new double [nzetapts];
	ptau_pts = new double [nptaupts];
	ptau_wts = new double [nptaupts];
	gauss (nspts, 0, smin, smax, s_pts, s_wts);
	gauss (nvpts, 0, -1.0, 1.0, v_pts, v_wts);
	gauss (nzetapts, 0, 0.0, M_PI, zeta_pts, zeta_wts);
	gauss (nptaupts, 0, 0.0, 20.0, ptau_pts, ptau_wts);
	
	return;
}

double Resonance_decay_contributions_to_pion_spectra(double pT, double y, double pphi, int reso_idx)
{
	double PkT, PkPhi, PkY;
	if (fabs(pT) < 1e-12) return (0.0);
	double ssum = 0.0;
	m = m_pion;
	double mT = sqrt(m*m+pT*pT);
	Mres = resonances.resonance_mass[reso_idx-1];
	Gamma = resonances.resonance_Gamma[reso_idx-1];
	m2 = resonances.resonance_decay_masses[reso_idx-1][0];
	m3 = resonances.resonance_decay_masses[reso_idx-1][1];
	br = resonances.resonance_total_br[reso_idx-1];
	sign = (double)resonances.resonance_sign[reso_idx-1];
	double s_min = (m2 + m3)*(m2 + m3);
	double s_max = (Mres - m)*(Mres - m);
	n_body = 2;
	if (m2 > 1.e-6 && m3 > 1.e-6) n_body = 3;
	
	Set_resonance_integration_points(s_min, s_max);
	Qfunc = get_Q();
	
	if (n_body == 2)
	{
		//then g(s) is delta-function, skip s-integration entirely
		double spt = m2*m2;
		double pstar = sqrt( ((Mres+m)*(Mres+m) - spt)*((Mres-m)*(Mres-m) - spt)/(2.0*Mres) );
		//cout << "pstar = " << pstar << endl;
		double gs = br/(4.*M_PI*pstar);
		double s_factor = gs;
		double Estar = sqrt(m*m + pstar*pstar);
		double psBmT = pstar / mT;
		double DeltaY = log(psBmT + sqrt(1.+psBmT*psBmT));
		double Yp = y + DeltaY;
		double Ym = y - DeltaY;
		double vsum = 0.0;
		for (int iv = 0; iv < nvpts; iv++)
		{
			double vpt = v_pts[iv];
			double Y = y + vpt*DeltaY;
			if (abs(Y) > PYmax) continue;
			PkY = Y;
			//cout << "Y = " << Y << endl;
			double mT_ch_Y_y = mT*cosh(vpt*DeltaY);
			double x2 = mT_ch_Y_y*mT_ch_Y_y - pT*pT;
			double v_factor = v_wts[iv]*DeltaY/sqrt(x2);
			double MTbar = Estar*Mres*mT_ch_Y_y/x2;
			double DeltaMT = Mres*pT*sqrt(Estar*Estar - x2)/x2;
			double MTp = MTbar + DeltaMT;
			double MTm = MTbar - DeltaMT;
			//cout << "MTm = " << MTm << endl;
			double zetasum = 0.0;
			//time (&rawtime);
			//timeinfo = localtime (&rawtime);
			//cerr << "Starting v-loop #" << iv << " at " << asctime(timeinfo);
			for (int izeta = 0; izeta < nzetapts; izeta++)
			{
				double zetapt = zeta_pts[izeta];
				double MT = MTbar + cos(zetapt)*DeltaMT;
				double zeta_factor = zeta_wts[izeta]*MT;
				double PT = sqrt(MT*MT - Mres*Mres);
				if (PT > PTmax) continue;
				PkT = PT;
				double temp_cos_PPhi_tilde = (mT*MT*cosh(Y-y) - Estar*Mres)/(pT*PT);
				//cout << "pT = " << pT << endl;
				//cout << "PT = " << PT << endl;
				double temp_sin_PPhi_tilde = sqrt(1. - temp_cos_PPhi_tilde*temp_cos_PPhi_tilde);
				//cout << "temp_cos_PPhi_tilde = " << temp_cos_PPhi_tilde << endl;
				//cout << "temp_sin_PPhi_tilde = " << temp_sin_PPhi_tilde << endl;
				double PPhi_tilde = place_in_range( atan2(temp_sin_PPhi_tilde, temp_cos_PPhi_tilde), PPhimin, PPhimax);
				double PPhi_tilde_shift = place_in_range( pphi + PPhi_tilde, PPhimin, PPhimax);
				PkPhi = PPhi_tilde_shift;
				double PPhi_tilde_shiftFLIP = place_in_range( pphi - PPhi_tilde, PPhimin, PPhimax);
				double * Pp = new double [4];
				double * Pm = new double [4];
				//probably not the most elegant set-up, but does the job for now...
				Pp[0] = MT * cosh(Y);
				Pp[1] = PT * cos(PPhi_tilde_shift);
				Pp[2] = PT * sin(PPhi_tilde_shift);
				Pp[3] = MT * sinh(Y);
				Pm[0] = Pp[0];
				Pm[1] = PT * cos(PPhi_tilde_shiftFLIP);
				Pm[2] = PT * sin(PPhi_tilde_shiftFLIP);
				Pm[3] = Pp[3];
				double ptausum = 0.0;
				for (int iptau = 0; iptau < nptaupts; iptau++)
				{
					double Csum = 0.0;
					double ptau_factor = Gamma * exp(-Gamma*ptau_pts[iptau]) * ptau_wts[iptau];
					for (int tempidx = 1; tempidx <= 2; tempidx++)
					{
						if (tempidx != 1)
							PkPhi = PPhi_tilde_shiftFLIP;		//also takes Pp --> Pm
						Csum += interpolate3D(PTpts, PPhipts, PYpts, resonance_spectra[reso_idx - 1],
									PkT, PkPhi, PkY, nPTpts, nPPhipts, nPYpts, 0, true, true);
						//cout << "Guessing this is where the problem is..." << endl
						//	<< "PkPhi = " << PkPhi << endl;
						//mytemp = interpolate3D(PTpts, PPhipts, PYpts, resonance_spectra[reso_idx - 1],
						//			PkT, PkPhi, PkY, nPTpts, nPPhipts, nPYpts, 0, true, true);
						//cerr << PkT << "   " << PkPhi << "   " << PkY << "   " << mytemp << endl;
						//Csum += mytemp;
					}
					ptausum += ptau_factor*Csum;
				}
				zetasum += zeta_factor*ptausum;
			}
			vsum += v_factor*zetasum;
		}
		ssum += Mres*s_factor*vsum;
	}
	else if (n_body == 3)
	{
		;
	}
	
	cout << "pion pT = " << pT << ": ssum = " << ssum << endl;
	return (ssum);
}

double get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-m)*(Mres-m);
	double sum = 0.;
	
	for (int is = 0; is < nspts; is++)
	{
		double sp = s_pts[is];
		double f1 = (Mres+m)*(Mres+m) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double g(double s)
{
	double g_res = br/(4.*M_PI);
	if (n_body == 3)
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		double g_res = pre_f * num / den;
	}

	return g_res;
}

//End of file
