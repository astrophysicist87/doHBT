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
#include "process_Cbar_and_Cavg_events.h"

using namespace std;

int main(int argc, char *argv[])
{
   if (argc==0)
   {
	cout << "Error: need to specify folder for processing!" << endl;
	exit(1);
   }

   //general pattern of typical results folder, without "-#" on end
   //usually something like "results"
   string folderstem = patch::to_string(argv[1]);
   //first event to include in ensemble
   int my_initial_event = atoi(argv[2]);
   //total number of event to include in ensemble
   int my_n_events = atoi(argv[3]);
   int my_final_event = my_initial_event + my_n_events - 1;

   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   string currentworkingdirectory = get_selfpath();
   string my_filename = "Processing_record_evs" + patch::to_string(my_initial_event) + "to" + patch::to_string(my_final_event) + ".txt";
   initialize_PRfile(currentworkingdirectory, my_filename);
   ostringstream filename_stream;
   filename_stream << currentworkingdirectory << "/" << my_filename;
   ofstream output(filename_stream.str().c_str(), ios::app);

   output << "/**********Processing output**********/" << endl;
   output << "Running in folder: " << currentworkingdirectory << endl;
   output << "Executable: " << argv[0] << endl;
   output << "Folderstem: " << argv[1] << endl;
   output << "Initial event: " << argv[2] << endl;
   output << "Number of events: " << argv[3] << endl;

   doHBT Source_function;

   /***Set requisite options for calculations***/
   //Since this code runs a level above the individual event results,
   //need to set a different variable than in process_event.cpp...
   Source_function.Set_runfolder(currentworkingdirectory);
   Source_function.Set_resultsfolder_stem(folderstem);
   Source_function.Set_use_delta_f(true);
   Source_function.Set_ofstream(output);
   Source_function.initial_event = my_initial_event;
   Source_function.n_events = my_n_events;
   //mass for pi^+ needs to be set explicitly here
   Source_function.Set_particle_mass(0.13957);

   //do the actual average HBT calculations...
   Source_function.Get_HBTradii_from_C_ev();
   output << "Finished." << endl;

   sw.toc();
   output << "Finished in " << sw.takeTime() << " sec." << endl;
   sw_total.toc();
   output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;
   output << "/**********End of processing output**********/" << endl;

   finalize_PRfile(currentworkingdirectory, my_filename);

   output.close();

   return 0;
}
