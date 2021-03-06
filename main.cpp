// DC9.cpp : 
/*	Code to estimate parameters in the presence of forecasts

Standard error. clear has to be added for DC7 and DC5
and some extra variables search for resize

INSTRUCTIONS:
1. The data must be arranged such that the socio-demographic variables appear first and then the distance from center followed by intensity
2. The square transform for discount factor is less robust than the exponential transformation due to their shapes. If Hessian is not invertible try fixing the discount factors
3. Note that other forecasts such as mandatory orders cannot be used since we'd then be assuming that we know the future.
4. The file names that end with 8 delete the first 8 time steps. This was found to give better predictions during LOOCV
5. The file names that ends with tod has the time of day variables. Identification issues are common because this data does not change over time. */

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <vector>
#include <chrono>
#include "Eigen//Dense"

using namespace Eigen;
using Eigen::MatrixXd;
using namespace std;

#include "households.h"
#include "functions.h"
#include "constants.h"

int main(int argc, char** argv)
{
	auto start = chrono::system_clock::now();
	//Writes console output to Output.txt
	/*FILE *stream;
	if ((stream = freopen("Output.txt", "w", stdout)) == NULL)
	exit(-1);*/

	//Define attributes of data
	//int num_households = 277;//This is hard coded to help with validation and needs to be modified for other datasets (iS it required?)
	vector<households> florence_household(CONST_NUM_HOUSEHOLDS_FLORENCE);
	//vector<households> michael_household(CONST_NUM_HOUSEHOLDS_MICHAEL);
	vector<households> household(CONST_NUM_HOUSEHOLDS_FLORENCE);
 	vector<string> covariate_name;//Name of variables
	set_household_data(florence_household, covariate_name, "florenceh");
	
	//set_household_data(michael_household, covariate_name, "michaelh");
	//set_intensity_forecasts(florence_household);
	//set_intensity_forecasts(michael_household);
	set_distance_forecasts(florence_household, "florence", "florence_locations");// set forecast values for Florence households
	//set_distance_forecasts(michael_household, "michael", "michael_locations");// set forecast values for Michael households
	
	
	//Assigning households based on hurricanes
	/*for (int i = 0; i < CONST_NUM_HOUSEHOLDS_FLORENCE + CONST_NUM_HOUSEHOLDS_MICHAEL; i++) {
		if (i <= CONST_NUM_HOUSEHOLDS_FLORENCE)
			household[i] = florence_household[i];
		    household = florence_household;
		if (i > CONST_NUM_HOUSEHOLDS_FLORENCE)
			household[i] = michael_household[i - CONST_NUM_HOUSEHOLDS_FLORENCE];
		    household = michael_household;
	}*/
	cout << "Displaying household data.." << endl;
	household = florence_household;
	household[1].display_data(covariate_name);

	cout << "order_perceivedness:" << endl;
	for (int i = 0; i < CONST_NUM_TIME_STEPS; i++) {
		cout << "time_step:" << i + 1 << endl;
		for (int k = 0; k < CONST_NUM_TIME_STEPS; k++) {
			cout << household[1].time_step[i].order_forecast[k] << " ";
		cout << endl;
		}
	}

	cout << "wind speed probabilities:" << endl;
	for (int i = 0; i < CONST_NUM_TIME_STEPS; i++) {
		cout << "time_step:" << i + 1 << endl;
		for (int k = 0; k < CONST_NUM_TIME_STEPS; k++) {
			cout << household[1].time_step[i].prob74_forecast[k] << " ";
			cout << endl;
		}
	}


	//Initialize the parameters to be estimated (we will then loop this within an optimization framework)
	vector<double> evac_estimate(CONST_NUM_COVARIATES, 1.0);//Evacuate
	vector<double> wait_estimate(CONST_NUM_COVARIATES, 1.0);//Defer decision to next time period
	vector<double> stay_estimate(CONST_NUM_COVARIATES, 1.0);//Stay indefinitely
	double discount_factor = 1.0;//Unconstrained discount variable. The actual discount_variable used is a transformed one

	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	household[1].display_value_functions();
	calculate_ev_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	household[1].display_derivatives();

	calculate_cv_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_probabilities(household);
	calculate_prob_partials(household);
	calculate_log_likelihood(household);
	calculate_likelihood_partials(household);

								 //Indicator variables used to indicate if parameters are to be included in the utility functions
	vector<int> evac_include;//Evacuate
	vector<int> wait_include;//Defer decision to next time period
	vector<int> stay_include;//Stay indefinitely

	set_include_indicators(evac_include, wait_include, stay_include, covariate_name);

	//Define vector of standard errors
	vector<double> standard_error;

	//Use the following function if the problem is being resolved to get the standard errors of the transformed variables
	//set_initial_values(evac_estimate, wait_estimate, stay_estimate, discount_factor);

	display_vector(evac_estimate, "Evacuate Parameters");
	display_vector(wait_estimate, "Wait Parameters");
	display_vector(stay_estimate, "Stay Parameters");
	cout << "Discount factor: " << discount_factor << endl;
	


	//Estimate parameters by maximizing log likelihood
	estimate_parameters(household, evac_estimate, wait_estimate, stay_estimate, evac_include, wait_include, stay_include, standard_error, discount_factor);

	//Print aggregate demand information
	vector<double> expected_demand(CONST_NUM_TIME_STEPS + 1, 0.0);
	double temp_demand;
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS + 1; t++)
		{
			temp_demand = 1;
			for (int s = 0; s < t; s++)
				temp_demand = temp_demand * household[i].time_step[s].prob_wait;

			if (t != CONST_NUM_TIME_STEPS)
				temp_demand = temp_demand * (1 - household[i].time_step[t].prob_wait);

			if (t == CONST_NUM_TIME_STEPS)
				cout << temp_demand << endl;

			expected_demand[t] += temp_demand;
		}
	}

	cout << endl << "Printing expected demand: " << endl;
	for (int i = 0; i < expected_demand.size(); i++)
		cout << expected_demand[i] << endl;

	//Display estimation results
	display_estimation_results(evac_estimate, wait_estimate, stay_estimate, evac_include, wait_include, stay_include, standard_error, covariate_name, discount_factor);

	//Use this function to write the estimate values to a text files that can be used to initialize values when resolving for estimating standard errors
	//output_parameter_values(evac_estimate, wait_estimate, stay_estimate, evac_include, wait_include, stay_include, discount_factor);

	//Calculate elasticities
	calculate_elasticities(household, evac_estimate, wait_estimate, stay_estimate, evac_include, wait_include, stay_include, discount_factor);

	/*double **avg_elasticity;//[t][j] averaged across all households
	avg_elasticity = new double*[CONST_NUM_TIME_STEPS];
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	avg_elasticity[t] = new double[CONST_NUM_COVARIATES];

	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
	avg_elasticity[t][j] = 0.0;
	for (int i = 0; i < num_households; i++)
	avg_elasticity[t][j] += household[i].time_step[t].covar_der_prob[j];

	avg_elasticity[t][j] = avg_elasticity[t][j] / num_households;
	}
	}

	//Print average elasticities
	cout << "Printing average elasticities: " <<  endl;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	cout << setw(CONST_EXTRA_WIDE_SPACING) << covariate_name[j];
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
	cout << endl << t;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	cout << setw(CONST_WIDE_SPACING) << avg_elasticity[t][j];
	}
	cout << endl;*/

	//Write the elasticities of number of years of residency
	/*household[0].output_elasticities(covariate_name, 7);
	household[182].output_elasticities(covariate_name, 7);
	household[110].output_elasticities(covariate_name, 7);

	//Print elasticitices and write the elasticities of distance from the center to a file
	household[0].output_elasticities(covariate_name, 8);
	household[182].output_elasticities(covariate_name, 8);
	household[110].output_elasticities(covariate_name, 8);

	//Finite difference elasticities wrt to Number of Vehicles (index = 1)
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 1, 0);
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 1, 182);
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 1, 110);

	//Finite difference elasticities wrt to Household Size (index = 2)
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 2, 0);
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 2, 182);
	output_finite_difference(household, evac_estimate, wait_estimate, stay_estimate, discount_factor, covariate_name, 2, 110);

	//Explaining elasticities
	display_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);

	//Vaidate results using leave-one-out cross validation
	validate_model(evac_include, wait_include, stay_include, covariate_name);*/

	//Compute running time 
	auto end = chrono::system_clock::now();
	chrono::duration<double> diff = end - start;
	cout << endl << "Total running time of the algorithm is: " << diff.count() << " seconds" << endl;

	system("pause");
	return 0;
}