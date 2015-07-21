#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>
#include<algorithm>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "../src/Stopwatch.h"
#include "../src/SV_parameters.h"
#include "../src/SV_readindata.h"
#include "../src/sourcevariances.h"
#include "../src/SV_generate_processing_record.h"
#include "../src/plumberglib.h"
#include "../src/sorter.h"
#include "compute_source_variances.h"

using namespace std;

int main(int argc, char *argv[])
{
   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   bool generatedcorrfuncs = false;
   string currentworkingdirectory = get_selfpath();
   int folderindex = get_folder_index(currentworkingdirectory);
   initialize_PRfile(currentworkingdirectory);

   ostringstream filename_stream;
   filename_stream << currentworkingdirectory << "/Processing_record.txt";
   ofstream output(filename_stream.str().c_str(), ios::app);

   output << "/**********Processing output**********/" << endl;
   output << "entering folder: " << currentworkingdirectory << endl;

   //load freeze out and particle information
   int FO_length = 0;
   int particle_idx=1;  //for pion+

   ostringstream decdatfile;
   output << "Loading the decoupling data...." << endl;
   decdatfile << currentworkingdirectory << "/decdat2.dat";
   output << decdatfile.str() << endl;
   FO_length=get_filelength(decdatfile.str().c_str());
   output << "Total number of freeze out fluid cell: " <<  FO_length << endl;

   //read the data arrays for the decoupling information
   FO_surf* FOsurf_ptr = new FO_surf[FO_length];
   read_decdat(FO_length, FOsurf_ptr, currentworkingdirectory, false);
   
   //read the positions of the freeze out surface
   read_surfdat(FO_length, FOsurf_ptr, currentworkingdirectory);
   
   //read the chemical potential on the freeze out surface
   int N_stableparticle;
   ifstream particletable("/home/plumberg.1/HBTPlumberg/EOS/EOS_particletable.dat");
   particletable >> N_stableparticle;
   double** particle_mu = new double* [N_stableparticle];
   for(int i=0; i<N_stableparticle; i++)
      particle_mu[i] = new double [FO_length];
   for(int i=0; i<N_stableparticle; i++)
      for(int j=0; j<FO_length; j++)
         particle_mu[i][j] = 0.0;
   if(N_stableparticle >0)
   {
      //if(hydropara_ptr->IEOS==7)       //for s95p_PCE
         read_decdat_mu(FO_length, N_stableparticle, particle_mu, currentworkingdirectory);
   }

   //read particle resonance decay table
   particle_info *particle = new particle_info [Maxparticle];
   int Nparticle=read_resonance(particle);
   output <<"read in total " << Nparticle << " particles!" << endl;
   output << "Calculating "<< particle[particle_idx].name << endl;
   if(N_stableparticle >0)
   {
      output << " EOS is partially chemical equilibrium " << endl;
      calculate_particle_mu(7, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
   }
   else
   {
	output << " EOS is chemical equilibrium. " << endl;
	for(int j=0; j<Nparticle; j++)
	{
		particle[j].mu = 0.0e0;
		for(int i=0; i<FO_length; i++)
			FOsurf_ptr[i].particle_mu[j] = 0.0e0;
	}
   }
	//calculate (semi-analytic approximation of) pure thermal spectra for all particle species
	calculate_thermal_particle_yield(Nparticle, particle, FOsurf_ptr[0].Tdec);
	//use this to estimate resonance-decay contributions from each particles species to final state particle, here, pion(+),
	//as well as set effective branching ratios
	compute_total_contribution_percentages(particle_idx, Nparticle, particle);
	//sort all particles by importance of their percentage contributions, then compute resonance SVs for only contributions up to some threshold
   sw.toc();
   output << "read in data finished!" << endl;
   output << "Used " << sw.takeTime() << " sec." << endl;

//	for(int i=0; i<FO_length; i++)
//	{
//		for(int j=0; j<Nparticle; j++)
//			cerr << FOsurf_ptr[i].particle_mu[j] << "    ";
//		cerr << endl;
//	}
//	if (1) exit(1);

   //HBT calculations begin ...
   double localy = 0.0e0;
   sw.tic();
     if(fabs(localy) > 1e-16)
     {
        output << "Case of y != 0 not yet supported.  Exiting..." << endl;
        return 0;
     }

//**********************************************************************************
//SELECT RESONANCES TO INCLUDE IN SOURCE VARIANCES CALCULATIONS
//sort resonances by importance and loop over all resonances needed to achieve a certain minimum percentage of total decay pions
	double threshold = 0.6;	//include only enough of the most important resonances to account for fixed fraction of total resonance-decay pion(+)s
				//threshold = 1.0 means include all resonance-decay pion(+)s,
				//threshold = 0.0 means include none of them
	vector<double> percent_contributions;
	for (int i = 0; i < Nparticle; i++)
		percent_contributions.push_back(particle[i].percent_contribution);
	vector<size_t> sorted_resonance_indices = ordered(percent_contributions);
	reverse(sorted_resonance_indices.begin(), sorted_resonance_indices.end());
	vector<int> chosen_resonance_indices;
	double running_total_percentage = 0.0;
	int count = 0;
	if (threshold < 1e-12)
	{
		output << "No resonances included." << endl;
	}
	else if (fabs(1. - threshold) < 1e-12)
	{
		for (int ii = 1; ii <= count; ii++)
			chosen_resonance_indices.push_back(sorted_resonance_indices[ii]);
		output << "All resonances included." << endl;
	}
	else
	{
		while (running_total_percentage <= threshold)
		{
			running_total_percentage += 0.01 * particle[sorted_resonance_indices[count]].percent_contribution;
			chosen_resonance_indices.push_back(sorted_resonance_indices[count]);
			count++;
		}
		if (chosen_resonance_indices.size() == 0)
		{
			output << "No resonances included!  Choose a higher threshold!" << endl;
			exit(1);
		}
		else
		{
			output << "Including the following " << count + 1 << " resonances (accounting for " << 100.*running_total_percentage
				<< "%, threshold = " << 100.*threshold << "%): " << endl;
			for (int ii = 1; ii <= count; ii++)
				output << "\t" << ii << ": " << particle[sorted_resonance_indices[ii - 1]].name << endl;
		}
	}
//END OF CODE SELECTING INCLUDED RESONANCES
//**********************************************************************************

//if (1) exit(1);

   SourceVariances Source_function(&particle[particle_idx], particle, Nparticle, FOsurf_ptr, chosen_resonance_indices);
   Source_function.Set_path(currentworkingdirectory);
   Source_function.Set_use_delta_f(true);

   Source_function.Update_sourcefunction(&particle[particle_idx], FO_length, particle_idx);

   Source_function.Set_ofstream(output);
   output << "Calculating HBT radii via source variances method..." << endl;
   Source_function.Analyze_sourcefunction(FOsurf_ptr);		//with previous function, this argument is redundant
   Source_function.Output_SVdN_dypTdpTdphi(folderindex);
   Source_function.Output_SVdN_dypTdpT(folderindex);
   //Source_function.Output_all_dN_dypTdpTdphi(folderindex);
   Source_function.Output_results(folderindex);
   output << "Finished calculating HBT radii via source variances method" << endl;



   sw.toc();
   output << "Finished in " << sw.takeTime() << " sec." << endl;
   sw_total.toc();
   output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;

   output << "/**********End of processing output**********/" << endl;

   output.close();

   checkforfiles_PRfile(currentworkingdirectory, folderindex, generatedcorrfuncs);

   finalize_PRfile(currentworkingdirectory);

   return 0;
}
