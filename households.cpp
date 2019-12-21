#include "stdafx.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

#include "households.h"
#include "constants.h"

households::households()
{

}

households::~households()
{

}

void households::display_data(vector<string> &covariate_name) const
{
	cout << setprecision(7);
	cout << "Household: " << id << endl;

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		cout << covariate_name[j] << ": ";
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			cout << time_step[t].covariate[j] << " ";
		cout << endl;
	}
	cout << "Action: ";
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		cout << time_step[t].action << " ";
	cout << endl << endl;


	/*cout << "Intensity forecasts: " << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << "Period " << t + 1 << endl;
		for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
		{
			cout << u + 1 << " ";
			for (int i = 0; i < 6; i++)
				cout << time_step[t].intensity_forecast[u][i] << " ";
			cout << endl;
		}
	}*/
	cout << endl;
}

void households::display_value_functions() const
{
	cout << setprecision(7);
	cout << endl << "Household: " << id << endl;

	cout << setw(CONST_WIDE_SPACING) << "time" << setw(CONST_WIDE_SPACING) << "cv_evac" << setw(CONST_WIDE_SPACING) << "cv_wait" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << setw(CONST_WIDE_SPACING) << t;
		cout << setw(CONST_WIDE_SPACING) << time_step[t].cv_evac;
		cout << setw(CONST_WIDE_SPACING) << time_step[t].cv_wait << endl;
	}

	cout << endl << "Printing exante value functions" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << "Time step: " << t << endl;
		for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
		{

			cout << setw(CONST_WIDE_SPACING) << time_step[t].exante_value[u];
			cout << endl;
		}
	}
}

void households::display_derivatives() const
{
	cout << "Displaying derivatives: ";
	cout << setprecision(7);
	cout << endl << "Household: " << id << endl;

	cout << setw(CONST_WIDE_SPACING) << "cv_evac_der" << setw(CONST_WIDE_SPACING) << "cv_wait_der" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		for (int j = 0; j < time_step[t].wait_der_cv_evac.size(); j++)
		{
			cout << " " << time_step[t].wait_der_cv_evac[j];
		}

		cout << "||";
		for (int j = 0; j < time_step[t].wait_der_cv_wait.size(); j++)
		{
			cout << " " << time_step[t].wait_der_cv_wait[j];
		}
		cout << endl;
	}

	cout << endl << "Printing derivatives of exante value functions wrt constants" << endl;
	for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
	{
	cout << "Future time step: " << u << endl;
	
	for (int j = 0; j < CONST_NUM_COVARIATES;j++)
		if(j<=1 || j==14)
			cout << setw(CONST_WIDE_SPACING) << time_step[0].evac_der_exante_value[u][j];
	    cout << endl;
		
	    
	
	}
	for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
	{
		cout << "Future time step " << u << ": ";
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
		if (j <= 1 || j == 14)
			cout << setw(CONST_WIDE_SPACING) << time_step[0].wait_der_exante_value[u][j];
	cout << endl;
	}
}

void households::display_probabilities() const
{
	cout << "Displaying probabilities: " << endl;
	cout << setprecision(7);
	cout << "Household: " << id << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		cout << setw(CONST_WIDE_SPACING) << t << ": " << time_step[t].prob_wait << endl;
	cout << endl;
}

void households::display_prob_partials() const
{
	cout << setprecision(7);
	cout << setw(CONST_WIDE_SPACING) << "Household: " << id << endl;

	cout << "-------------------" << endl;
	cout << "Prob_wait partials:" << endl;
	cout << "-------------------" << endl;

	cout << "d prob_wait / d evac coefficients:" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << setw(CONST_WIDE_SPACING) << t << ": ";
		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			cout << setw(CONST_WIDE_SPACING) << time_step[t].evac_der_prob_wait[j];
		cout << endl;
	}

	cout << "d prob_wait / d wait coefficients:" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << setw(CONST_WIDE_SPACING) << t << ": ";
		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			cout << setw(CONST_WIDE_SPACING) << time_step[t].wait_der_prob_wait[j];
		cout << endl;
	}

	cout << "d prob_wait / d stay coefficients:" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << setw(CONST_WIDE_SPACING) << t << ": ";
		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			cout << setw(CONST_WIDE_SPACING) << time_step[t].stay_der_prob_wait[j];
		cout << endl;
	}

	cout << "d prob_wait / d discount factor: " << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << setw(CONST_WIDE_SPACING) << t << ": ";
		cout << setw(CONST_WIDE_SPACING) << time_step[t].dfact_der_prob_wait;
		cout << endl;
	}
}

void households::display_elasticities(vector<string> &covariate_name) const
{
	cout << "Displaying elasticities: " << endl;
	cout << setprecision(7);
	cout << "Household: " << id << endl;
	cout << "Time";
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
		cout << setw(CONST_EXTRA_WIDE_SPACING) << covariate_name[j];
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		cout << endl << t;
		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			cout << setw(CONST_WIDE_SPACING) << time_step[t].covar_der_prob[j];
	}
	cout << endl;
}

void households::output_elasticities(vector<string> &covariate_name, int covariate_index) const
{
	//Write elasticities of a certain variable to an external file
	ofstream outfile("./Results/" + covariate_name[covariate_index] + "_" + to_string(id - 1) + ".dat", std::ios::trunc);
	outfile.precision(7);
	if (outfile.is_open())
	{
		outfile << "Time \t Prob" << endl;
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			outfile << t + 1 << "\t";
			outfile << -time_step[t].covar_der_prob[covariate_index];
			outfile << endl;
		}
	}
	else
		cout << "ERROR: Could not create output dat file " << endl;
}