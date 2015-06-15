#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "../src/Stopwatch.h"
#include "../src/parameters.h"
#include "../src/readindata.h"
#include "../src/doHBT.h"
#include "../src/generate_processing_record.h"
#include "../src/plumberglib.h"
#include "resonance_estimates.h"

using namespace std;

int main(int argc, char *argv[])
{
   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   string currentworkingdirectory = get_selfpath();
   initialize_PRfile(currentworkingdirectory);

   ostringstream filename_stream;
   filename_stream << currentworkingdirectory << "/Processing_record.txt";
   ofstream output(filename_stream.str().c_str(), ios::app);

   output << "/**********Processing output**********/" << endl;
   output << "entering folder: " << currentworkingdirectory << endl;

   //load freeze out information
   int FO_length = 0;
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
   //print_particle_stability(particle, Nparticle);
   output <<"read in total " << Nparticle << " particles!" << endl;
   if(N_stableparticle >0)
   {
      output << " EOS is partially chemical equilibrium " << endl;
      calculate_particle_mu(7, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
   }
   else
   {
      output << " EOS is chemical equilibrium. " << endl;
      for(int i=0; i<FO_length; i++)
      for(int j=0; j<Nparticle; j++)
         FOsurf_ptr[i].particle_mu[j] = 0.0e0;
   }
   sw.toc();
   output << "read in data finished!" << endl;
   output << "Used " << sw.takeTime() << " sec." << endl;
   
   sw.tic();
	for (int ipi = 0; ipi < 100; ipi++)
		cout << ipi << "   " << particle[ipi].name << endl;
	//int particle_idx = 1;  //for pion+
	//int particle_idx = 4;  //for kaon+
	int particle_idx = 17;  //for proton
	double * all_particle_thermal = new double [Nparticle];
	double * percentages = new double [Nparticle];
	double * effective_widths = new double [Nparticle];
	estimate_resonance_thermal(Nparticle, particle, FOsurf_ptr[0].Tdec, all_particle_thermal);
	output << "Finished estimating thermal spectra of resonances" << endl;
	output << "Getting all relative contributions..." << endl;
	compute_total_contribution_percentages(particle_idx, Nparticle, particle, all_particle_thermal, percentages, effective_widths);
	output << "Finished computing total contribution percentages." << endl;
   ostringstream output_filename_stream;
   output_filename_stream << currentworkingdirectory << "/" << particle[particle_idx].name << "_output.dat";
   ofstream finaloutput(output_filename_stream.str().c_str(), ios::app);
	for (int i = 0; i < Nparticle; i++)
		finaloutput << particle[i].name << "   " << percentages[i] << "   " << all_particle_thermal[i]
			<< "   " << effective_widths[i] << "   " << effective_widths[i] * all_particle_thermal[i] << endl;
		//output << particle[i].name << ": " << percentages[i] << "%" << endl;
   sw.toc();
   output << "Finished in " << sw.takeTime() << " sec." << endl;
   sw_total.toc();
   output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;

   output << "/**********End of processing output**********/" << endl;

   output.close();

   finalize_PRfile(currentworkingdirectory);

   return 0;
}
