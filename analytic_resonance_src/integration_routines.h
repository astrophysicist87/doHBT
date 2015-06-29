#ifndef INTEGRATION_ROUTINES_H
#define INTEGRATION_ROUTINES_H

using namespace std;

double integrate(double(*function_original) (vector<double> *, void *), void * params_ptr,
		vector<double> * lower_limits_ptr, vector<double> * upper_limits_ptr, vector<int> * intervals);


double integrate(double(*function_original) (vector<double> *, void *), void * params_ptr,
		vector<double> * lower_limits_ptr, vector<double> * upper_limits_ptr, vector<int> * intervals)
{
	double result = 0.;
	const int dimint = (*lower_limits_ptr).size();		//the dimension of the integral we need to do
	vector< vector<double> * > all_local_xpts_ptrs (dimint);
	vector< vector<double> * > all_local_wts_ptrs (dimint);

	vector<double> halflengths (dimint), centers (dimint);	//vectors to hold the halflengths and center of each respective integration interval

	for (int i = 0; i < dimint; i++)
	{
		switch ((*intervals)[i])
		{
		case 0:	//finite interval (a,b)
			halflengths[i] = 0.5 * ((*upper_limits_ptr)[i] - (*lower_limits_ptr)[i]);
			centers[i] = 0.5 * ((*upper_limits_ptr)[i] + (*lower_limits_ptr)[i]);
			all_local_xpts_ptrs[i] = xi_ptr;
			all_local_wts_ptrs[i] = wi_ptr;
			break;
		case 1:	//half-infinite interval (0,inf)
			halflengths[i] = 1.;
			centers[i] = 0.;
			all_local_xpts_ptrs[i] = xi_0pinf_ptr;
			all_local_wts_ptrs[i] = wi_0pinf_ptr;
			break;
		case 2:	//full-infinite interval (-inf,inf)
			halflengths[i] = 1.;
			centers[i] = 0.;
			all_local_xpts_ptrs[i] = xi_minfpinf_ptr;
			all_local_wts_ptrs[i] = wi_minfpinf_ptr;
			break;
		}
	}
	for (int i = 0; i < pow(order, dimint); i++)
	{
		double wt_factor = 1.;
		vector<int> indices = convert_base(i, order, dimint);	//treats all points for multi-dimensional integration as one giant vector
									//indices vector now holds index value of loop variable i in order**dimint grid
									//total number of points is order**dimint
		vector<double> xlocs (dimint);
		for (int j = 0; j < dimint; j++)
		{
			xlocs[j] = centers[j] + halflengths[j] * (*(all_local_xpts_ptrs[j]))[indices[j]];
			wt_factor *= halflengths[j] * (*(all_local_wts_ptrs[j]))[indices[j]];
		}
		result += wt_factor * function_original(&xlocs, params_ptr);
	}

	return (result);
}

//End of file

#endif
