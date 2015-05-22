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
#include "process_event.h"

using namespace std;

int main(int argc, char *argv[])
{
   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   bool generatedcorrfuncs = false;
bool get_plane_angle_only = true;
   //instantiate doHBT class
   doHBT Source_function;

   //string currentworkingdirectory = mypath + patch::to_string(folderindex);
   //string currentworkingdirectory = ".";
   string currentworkingdirectory = get_selfpath();
   int folderindex = get_folder_index(currentworkingdirectory);
   initialize_PRfile(currentworkingdirectory);
   Source_function.Set_path(currentworkingdirectory);
   Source_function.Set_use_delta_f(false);

   ostringstream filename_stream;
   filename_stream << currentworkingdirectory << "/Processing_record.txt";
   ofstream output(filename_stream.str().c_str(), ios::app);

   output << "/**********Processing output**********/" << endl;
   //cout << "entering folder: " << currentworkingdirectory << endl;
   //output << "re-processing with delta_f = 0" << endl;
   output << "entering folder: " << currentworkingdirectory << endl;

   //load freeze out information
   int FO_length = 0;
   ostringstream decdatfile;
   //cout << "Loading the decoupling data...." << endl;
   output << "Loading the decoupling data...." << endl;
   decdatfile << currentworkingdirectory << "/decdat2.dat";
   //cout << decdatfile.str() << endl;
   output << decdatfile.str() << endl;
   FO_length=get_filelength(decdatfile.str().c_str());
   //cout << "Total number of freeze out fluid cell: " <<  FO_length << endl;
   output << "Total number of freeze out fluid cell: " <<  FO_length << endl;

//output << "am i even running the right code???" << endl;

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
   //cout <<"read in total " << Nparticle << " particles!" << endl;
   output <<"read in total " << Nparticle << " particles!" << endl;
   if(N_stableparticle >0)
   {
      //cout << " EOS is partially chemical equilibrium " << endl;
      output << " EOS is partially chemical equilibrium " << endl;
//cerr << "plumberg: (&FOsurf_ptr[0])->particle_mu[1] = " << (&FOsurf_ptr[0])->particle_mu[1] << endl;
      calculate_particle_mu(7, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
//cerr << "plumberg: (&FOsurf_ptr[0])->particle_mu[1] = " << (&FOsurf_ptr[0])->particle_mu[1] << endl;
   }
   else
   {
      //cout << " EOS is chemical equilibrium. " << endl;
      output << " EOS is chemical equilibrium. " << endl;
      for(int i=0; i<FO_length; i++)
      for(int j=0; j<Nparticle; j++)
         FOsurf_ptr[i].particle_mu[j] = 0.0e0;
   }
   sw.toc();
   //cout << "read in data finished!" << endl;
   output << "read in data finished!" << endl;
   //cout << "Used " << sw.takeTime() << " sec." << endl;
   output << "Used " << sw.takeTime() << " sec." << endl;

   //HBT calculations begin ...
   int particle_idx=1;  //for pion+
   //cout << "Calculating "<< particle[particle_idx].name << endl;
   output << "Calculating "<< particle[particle_idx].name << endl;
   
   double localy = 0.0e0;
   sw.tic();
     if(fabs(localy) > 1e-16)
     {
        //cout << "not support y not equals 0 yet! Bye bye!" << endl;
        output << "not support y not equals 0 yet! Bye bye!" << endl;
        return 0;
     }

//cerr << "plumberg: (&FOsurf_ptr[0])->particle_mu[particle_id] = " << (&FOsurf_ptr[0])->particle_mu[particle_idx] << endl;

   Source_function.Update_sourcefunction(&particle[particle_idx], FO_length, particle_idx);

   Source_function.Set_ofstream(output);
   //Source_function.Set_folderindex(folderindex);
   //Source_function.Set_global_currentworkingdirectory(currentworkingdirectory);

if (get_plane_angle_only)
{
	output << "Getting plane angle only to save time..." << endl;
	Source_function.Determine_plane_angle(FOsurf_ptr);
	Source_function.Output_ev_plane_psi(folderindex);
	Source_function.Output_ev_plane_psis(folderindex);
	Source_function.Output_dN_dypTdpTdphi(folderindex);
	Source_function.Output_dN_dypTdpT(folderindex);
	Source_function.Output_ev_mean_pT(folderindex);
	output << "Finished all." << endl;
	finalize_PRfile(currentworkingdirectory);
	output.close();
	
	return(0);
}

   output << "Calculating HBT radii via source variances method..." << endl;
   Source_function.Analyze_sourcefunction(FOsurf_ptr);
   Source_function.quick_Analyze_sourcefunction_vars();
   Source_function.Output_ev_plane_psis(folderindex);
   Source_function.Output_dN_dypTdpTdphi(folderindex);
   Source_function.Output_dN_dypTdpT(folderindex);
   Source_function.Output_ev_mean_pT(folderindex);
   Source_function.Output_ev_plane_psi(folderindex);
   Source_function.Output_results(folderindex);
   Source_function.Output_Svars_results(folderindex);
   Source_function.Output_ev_anisotropic_flows(folderindex);
   Source_function.Output_ev_anisotropic_flows_pTdiff(folderindex);
   output << "Finished calculating HBT radii via source variances method" << endl;



//cutting out GF calculation for now, just to process individual events via SV method
//   output << "Beginning calculation of HBT radii via Gaussian fit method" << endl;

//   Source_function.Get_GF_HBTradii(FOsurf_ptr, folderindex);

//   Source_function.Output_GF_results(folderindex);

/*   ostringstream OPHBTGFradii;
   OPHBTGFradii << currentworkingdirectory << "/HBTGFradii_ev" << folderindex << ".dat" ;
   ofstream OPfile(OPHBTGFradii.str().c_str());
   
   double localy = 0.0e0;

         OPfile << scientific << setprecision(7) << setw(15)
                << localp_T << "  " << localp_phi << "  "  
                << HBT_hadron.get_lambda_Correl() << "  " << HBT_hadron.get_lambda_Correl_err() << "  "
                << HBT_hadron.get_R2out_Correl() << "  " << HBT_hadron.get_R2out_Correl_err() << "  "
                << HBT_hadron.get_R2side_Correl() << "  " << HBT_hadron.get_R2side_Correl_err() << "  "
                << HBT_hadron.get_R2long_Correl() << "  " << HBT_hadron.get_R2long_Correl_err() << "  "
                << HBT_hadron.get_R2os_Correl() << "  " << HBT_hadron.get_R2os_Correl_err()
                << endl;
         sw.toc();
         //cout << "Finished in " << sw.takeTime() << " sec." << endl;
         cerr << "Finished in " << sw.takeTime() << " sec." << endl;

   OPfile.close();*/








   sw.toc();
   //cout << "Finished in " << sw.takeTime() << " sec." << endl;
   output << "Finished in " << sw.takeTime() << " sec." << endl;
   sw_total.toc();
   //cout << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;
   output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;

   output << "/**********End of processing output**********/" << endl;

   output.close();

   checkforfiles_PRfile(currentworkingdirectory, folderindex, generatedcorrfuncs);

   finalize_PRfile(currentworkingdirectory);

   return 0;
}
