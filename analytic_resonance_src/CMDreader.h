#ifndef CMDREADER_H
#define CMDREADER_H

#include<string>
#include<sstream>
#include"corr_fcn_from_toyS_v2.h"
using namespace std;


int Read_in_parameters_from_CMD(int argc, char *argv[])
{
	eps_2_bar_cmd = eps_2_bar_default;
	v_2_bar_cmd = v_2_bar_default;
	psi_2_bar_cmd = psi_2_bar_default;
	eta_f_cmd = eta_f_default;
	Rad_cmd = Rad_default;
	T0_cmd = T0_default;

	//take arguments from command-line
	if (argc == 1) { //if no arguments included in command-line, run with defaults
		cout << "Running with defaults" << endl;
	}
	else if (is_multiple(argc , 2)) {
		//if incorrect number of arguments, end program and return (1)
		cerr << "Incorrect number of arguments: expected (function_call) (int param_key) (double param_val)" << endl
		     << "param_key: 1 - eps_2_bar, 2 - v_2_bar," << endl
		     << "3 - psi_2_bar, 4 - eta_f, 5 - Rad, 6 - T0" << endl;
	
		return (1);
	}
	else {//otherwise, update any parameters set from command-line and move on
	int arg_index = 0;
	do {
		switch (atoi(argv[arg_index+1]))
		{
			case 1:
				eps_2_bar_cmd = atof(argv[arg_index+2]);
				break;
			case 2:
				v_2_bar_cmd = atof(argv[arg_index+2]);
				break;
			case 3:
				psi_2_bar_cmd = atof(argv[arg_index+2])*PI;
				break;
			case 4:
				eta_f_cmd = atof(argv[arg_index+2]);
				break;
			case 5:
				Rad_cmd = atof(argv[arg_index+2]);
				break;
			case 6:
				T0_cmd = atof(argv[arg_index+2]);
				break;
		}
	
	arg_index += 2;  //advance to the next two input arguments
	  } while (arg_index < argc-1);
	}

	//finally, set actual parameters used in main code from finally updated command-line parameters
	eps_2_bar = eps_2_bar_cmd;
	v_2_bar = v_2_bar_cmd;
	psi_2_bar = psi_2_bar_cmd;
	eta_f = eta_f_cmd;
	Rad = Rad_cmd;
	T0 = T0_cmd;

	return (0);
}
	
#endif
