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
#include "../src/SV_parameters.h"
#include "../src/SV_readindata.h"
#include "../src/sourcevariances.h"
#include "../src/SV_generate_processing_record.h"
#include "../src/plumberglib.h"
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

   //HBT calculations begin ...
   int particle_idx=1;  //for pion+
   output << "Calculating "<< particle[particle_idx].name << endl;
   
   double localy = 0.0e0;
   sw.tic();
     if(fabs(localy) > 1e-16)
     {
        output << "not support y not equals 0 yet! Bye bye!" << endl;
        return 0;
     }

   SourceVariances Source_function(&particle[particle_idx]);
   Source_function.Set_path(currentworkingdirectory);
   Source_function.Set_use_delta_f(true);

   Source_function.Update_sourcefunction(&particle[particle_idx],FO_length, particle_idx);

   Source_function.Set_ofstream(output);
   output << "Calculating HBT radii via source variances method..." << endl;
   Source_function.Set_current_FOsurf_ptr(FOsurf_ptr);
   //Source_function.Analyze_sourcefunction_check(FOsurf_ptr);		//with previous function, this argument is redundant
   Source_function.Analyze_sourcefunction(FOsurf_ptr);		//with previous function, this argument is redundant
   Source_function.Output_SVdN_dypTdpTdphi(folderindex);
   Source_function.Output_SVdN_dypTdpT(folderindex);
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
