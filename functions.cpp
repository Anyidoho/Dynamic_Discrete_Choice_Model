#include "stdafx.h"

#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <queue>
#include <deque>
#include <cmath>
#include <numeric>
#include <string>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

#include "functions.h"
#include "households.h"
#include "constants.h"

//Function that converts text to columns with two possible delimiters (can be exteded to include more delimiters)
vector<string> split_string(const char *line, char c1, char c2)
{
	vector<string> result;

	do
	{
		const char *begin = line;

		while (*line != c1 && *line != c2 && *line)
			line++;

		result.push_back(string(begin, line));
	} while (0 != *line++);

	return result;
}

//Function to compute the distance between two lat longs using the Haversine formula
double compute_distance(double lat1, double long1, double lat2, double  long2)
{
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = (lat1 * 3.14159265358979323846 / 180);
	lon1r = (long1 * 3.14159265358979323846 / 180);
	lat2r = (lat2 * 3.14159265358979323846 / 180);
	lon2r = (long2 * 3.14159265358979323846 / 180);

	u = sin((lat2r - lat1r) / 2);
	v = sin((lon2r - lon1r) / 2);

	return (2.0 * 6371.0 * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v)))*0.01;
}

double product_rule(vector<double> &prob, vector<double> &der_prob)
{
	double sum = 0;
	double product;
	for (int i = 0; i < der_prob.size(); i++)
	{
		product = 1;
		for (int j = 0; j < prob.size(); j++)
		{
			if (i != j)
				product = product * prob[j];
		}
		sum += der_prob[i] * product;
	}

	return sum;
}

void set_household_data(vector<households> &household, vector<string> &covariate_name, string file_name)
{
	cout << "Reading household data.." << endl;

	//Initialize sizes
	for (int i = 0; i < household.size(); i++)
	{
		household[i].id = i + 1;
		household[i].time_step.resize(CONST_NUM_TIME_STEPS);

		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].covariate.resize(CONST_NUM_COVARIATES);
		//order forecast
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].order_forecast.resize(CONST_NUM_TIME_STEPS);

		//mandatory order
		
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].manorder_forecast.resize(CONST_NUM_TIME_STEPS);

		//voluntary order
		
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].volorder_forecast.resize(CONST_NUM_TIME_STEPS);

		//distance to landfall
		
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].distlandfall_forecast.resize(CONST_NUM_TIME_STEPS);

		//intensity forecast
	
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].intensity_forecast.resize(CONST_NUM_TIME_STEPS);

		//flood level forecast
		
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].floodlevel_forecast.resize(CONST_NUM_TIME_STEPS);

		//wind speed probability forecast
		
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].prob74_forecast.resize(CONST_NUM_TIME_STEPS);

	}
	

	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].exante_value = new double[CONST_NUM_TIME_STEPS];

			household[i].time_step[t].evac_der_exante_value = new double*[CONST_NUM_TIME_STEPS];
			household[i].time_step[t].wait_der_exante_value = new double*[CONST_NUM_TIME_STEPS];
			household[i].time_step[t].stay_der_exante_value = new double*[CONST_NUM_TIME_STEPS];

			household[i].time_step[t].dfact_der_exante_value = new double[CONST_NUM_TIME_STEPS];

			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
			{
				household[i].time_step[t].evac_der_exante_value[u] = new double[CONST_NUM_COVARIATES];
				household[i].time_step[t].wait_der_exante_value[u] = new double[CONST_NUM_COVARIATES];
				household[i].time_step[t].stay_der_exante_value[u] = new double[CONST_NUM_COVARIATES];			
				
			}
		}
	}

	int line_number = 0;
	int start_line = 2;
	int index, ts_index;
	string line, temp_number;
	vector<string> token;
	ifstream file(file_name+".txt");//This is the text file in which income indicators are replaced by mean values

	covariate_name.clear();
	if (file.is_open())
	{
		do
		{
			line_number++;
			getline(file, line);
			//Store the names of the covariates
			if (line_number == 1)
			{
				//cout << line << endl;
				token = split_string(line.c_str(), '\t', ' ');
				covariate_name.push_back("Constant");
				for (int j = 2; j < CONST_NUM_COVARIATES + 1; j++)
					covariate_name.push_back(token[j].c_str());
			}
			//Store the covariates
			else if (line_number >= start_line)
			{
				//cout << line << endl;
				token = split_string(line.c_str(), '\t', ' ');
				index = atoi(token[0].c_str()) - 1;
				ts_index = atoi(token[1].c_str()) - 1;

				//Store constant
				household[index].time_step[ts_index].covariate[0] = 1;

				//Store other covariates
				for (int j = 1; j < CONST_NUM_COVARIATES; j++)
					household[index].time_step[ts_index].covariate[j] = atof(token[j + 1].c_str());
				//orderforecast
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if(k< ts_index)
						household[index].time_step[ts_index].order_forecast[k] = -1;
					else
						household[index].time_step[ts_index].order_forecast[k]= atof(token[12].c_str()); 

				//mandatory forecast
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].manorder_forecast[k] = -1;
					else
					household[index].time_step[ts_index].manorder_forecast[k] = atof(token[13].c_str());

				//voluntary order forecast
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].volorder_forecast[k] = -1;
					else
					household[index].time_step[ts_index].volorder_forecast[k] = atof(token[14].c_str());
				
				//distance_landfall
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].distlandfall_forecast[k] = -1;
					else
					household[index].time_step[ts_index].distlandfall_forecast[k] = atof(token[15].c_str());
				
				//intensity forecast
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].intensity_forecast[k] = -1;
					else
					household[index].time_step[ts_index].intensity_forecast[k] = atof(token[16].c_str());
				
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].floodlevel_forecast[k] = -1;
					else
					household[index].time_step[ts_index].floodlevel_forecast[k] = atof(token[17].c_str());
				
				for (int k = 0; k < CONST_NUM_TIME_STEPS; k++)
					if (k < ts_index)
						household[index].time_step[ts_index].prob74_forecast[k] = -1;
					else
					household[index].time_step[ts_index].prob74_forecast[k] = atof(token[19].c_str());
				//Store action;
				household[index].time_step[ts_index].action = atoi(token[CONST_NUM_COVARIATES + 1].c_str());

			}
		} while (file.good());
	}
	else
		cout << endl << "Unable to open input file";
	file.close();
}
/*
void set_intensity_forecasts(vector<households> &household)
{
	cout << "Reading intensity forecasts.." << endl;

	//Initialize sizes of intensity forecasts
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].intensity_forecast = new double*[CONST_NUM_TIME_STEPS];
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
				household[i].time_step[t].intensity_forecast[u] = new double[6];
		}
	}

	//Initialize sizes of intensity forecasts
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
				for (int s = 0; s < CONST_NUM_STATES; s++)
					household[i].time_step[t].intensity_forecast[u][s] = 0;
		}
	}

	//Read data from text files
	/*int pos;//Position of current time period in the text file
	int t, u;
	string line;
	vector<string> token;
	ifstream file("./intensity_forecast.txt");

	while (getline(file, line))
	{
		//Record the current time period
		if (line.find("Period") != string::npos)
		{
			pos = line.find("Period");
			stringstream(line.substr(pos + 6)) >> t;
			t = t - 1;
		}
		else
		{
			token = split_string(line.c_str(), '\t', ' ');
			u = atoi(token[0].c_str()) - 1;

			for (int h = 0; h < household.size(); h++)
				for (int i = 0; i < CONST_NUM_STATES; i++)
					household[h].time_step[t].intensity_forecast[u][i] = atof(token[i + 1].c_str());
		}
	}
	file.close();*
}
*/
void set_distance_forecasts(vector<households> &household, string file_name1, string file_name2)// added new argument as file name
{
	cout << "Reading distance forecasts.." << endl;

	//Initialize sizes of distance forecasts
	for (int i = 0; i < household.size(); i++)
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].distance_forecast.resize(CONST_NUM_TIME_STEPS);

	//Store the predicted location of the hurricane
	double*** pred_location;
	pred_location = new double**[CONST_NUM_TIME_STEPS];
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		pred_location[t] = new double*[CONST_NUM_TIME_STEPS];
		for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
			pred_location[t][u] = new double[2];
	}

	//Read future distances of the center of the hurricane in a separate variable
	int pos;//Position of current time period in the text file
	int t, u;
	string line;
	vector<string> token;
	ifstream file( file_name1 + ".txt"); //pass in file name as argument

	while (getline(file, line))
	{
		//Record the current time period
		if (line.find("Period") != string::npos)
		{
			pos = line.find("Period");
			stringstream(line.substr(pos + 6)) >> t;
			t = t - 1;
		}
		else
		{
			token = split_string(line.c_str(), '\t', ' ');
			u = atoi(token[0].c_str()) - 1;

			pred_location[t][u][0] = atof(token[1].c_str());
			pred_location[t][u][1] = -atof(token[2].c_str());
		}
	}
	file.close();
	//Print read values 
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
	cout << "Period: " << t << endl;
	for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
	cout << u << "\t" << pred_location[t][u][0] << "\t" << pred_location[t][u][1] << endl;
	}

	//Store household lat longs 
	int h;
	ifstream file2(file_name2 + ".txt");

	while (getline(file2, line))
	{
		token = split_string(line.c_str(), '\t', ' ');
		h = atoi(token[0].c_str()) - 1;

		household[h].latitude = atof(token[1].c_str());
		household[h].longitude = atof(token[2].c_str());
	}
	file2.close();

	//Print household lat longs
	/*cout << "Printing lat longs of households: " << endl;
	for (int i = 0; i < household.size(); i++)
	{
	cout << i + 1 << "\t" << household[i].latitude << "\t" << household[i].longitude << endl;
	}*/

	//Compute distances from each household and store them
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
			{
				if (u >= t)
					household[i].time_step[t].distance_forecast[u] = compute_distance(household[i].latitude, household[i].longitude, pred_location[t][u][0], pred_location[t][u][1]);
				else
					household[i].time_step[t].distance_forecast[u] = -1;
			}
		}
	}

	//Print distances for one household
	cout << "Printing distances in km from household 1" << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
	cout << "Period " << t << endl;
	for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
	cout << household[1].time_step[t].distance_forecast[u] << endl;
	}
}



void set_include_indicators(vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, vector<string> &covariate_name)
{
	cout << "Reading include indicators.." << endl;

	int line_number = 0;
	int start_line = 2;
	int index, ts_index;
	string line, temp_number;
	vector<string> token;
	ifstream file("./include.txt");
	if (file.is_open())
	{

		do
		{
			line_number++;
			getline(file, line);
			//Store the covariates
			if (line_number >= start_line)
			{
				cout << line << endl;
				token = split_string(line.c_str(), '\t', ' ');
				evac_include.push_back(atoi(token[1].c_str()));
				wait_include.push_back(atoi(token[2].c_str()));
				stay_include.push_back(atoi(token[3].c_str()));
			}
		} while (file.good());
	}
	else
		cout << endl << "Unable to open input file";
	file.close();
	cout << "Printing include indicators.." << endl;
	cout << setw(CONST_EXTRA_WIDE_SPACING) << " ";
	cout << setw(CONST_WIDE_SPACING) << "Evacuate";
	cout << setw(CONST_WIDE_SPACING) << "Wait";
	cout << setw(CONST_WIDE_SPACING) << "Stay" << endl;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		cout << setw(CONST_EXTRA_WIDE_SPACING) << covariate_name[j];
		cout << setw(CONST_WIDE_SPACING) << evac_include[j];
		cout << setw(CONST_WIDE_SPACING) << wait_include[j];
		cout << setw(CONST_WIDE_SPACING) << stay_include[j] << endl;
	}
	cout << endl;
}

void calculate_value_functions(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	//Step 0: Initialize sizes and values
	for (int i = 0; i < household.size(); i++)
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
					household[i].time_step[t].exante_value[u] = 0.0;

	//Step 1: Set discount factors
	double trans_dfact;
#ifdef CONST_NO_TRANSFORM
	trans_dfact = discount_factor;
#endif

#ifdef CONST_LOGIT_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + exp(-discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + pow(discount_factor, -2)));
#endif

	//Step 2: Initialize values to zeros
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].cv_evac = 0.0;
			household[i].time_step[t].cv_wait = 0.0;
		}
	}

	//Step 3: For every household, update the ex ante value functions (Backward Induction)
	double temp, temp_e, temp_w;
	for (int i = 0; i < household.size(); i++) //For every household
	{
		//int t = 0;
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++) //For every time step
		{
			for (int u = CONST_NUM_TIME_STEPS - 1; u > t; u--) //Scan previous time steps to compute ex ante value function
			{
				
					if (u == CONST_NUM_TIME_STEPS - 1)
					{
						temp = 0.0;
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += stay_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						temp = log(temp_e + temp_w);
						temp += CONST_EULER;
						household[i].time_step[t].exante_value[u] = temp;
					}
					else
					{
						temp = 0.0;
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

					
						
						temp_w +=household[i].time_step[t].exante_value[u + 1];
						temp_w = trans_dfact * temp_w;
						
						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += wait_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						temp = log(temp_e + temp_w);
						temp += CONST_EULER;
						household[i].time_step[t].exante_value[u] = temp;
					}
				}
			}
     	}
	

	//Step 4:  Update the conditional value functions 
	//Conditional value functions of the evac option
	for (int i = 0; i < household.size(); i++)
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				household[i].time_step[t].cv_evac += evac_estimate[j] * household[i].time_step[t].covariate[j];

	//Conditional value functions of the wait/stay option
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			if (t == CONST_NUM_TIME_STEPS - 1)//If it is the last time period
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
					household[i].time_step[t].cv_wait += stay_estimate[j] * household[i].time_step[t].covariate[j];
			}
			else
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
					household[i].time_step[t].cv_wait += wait_estimate[j] * household[i].time_step[t].covariate[j];
				    household[i].time_step[t].cv_wait += trans_dfact * household[i].time_step[t].exante_value[t + 1];
			}
		}
	}
}

//Derivatives of ex ante value functions
void calculate_ev_partials(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	//Step 1: Set discount factors
	double trans_dfact, der_trans_dfact;
#ifdef CONST_NO_TRANSFORM
	trans_dfact = discount_factor;
	der_trans_dfact = 1;
#endif

#ifdef CONST_LOGIT_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + exp(-discount_factor)));
	der_trans_dfact = ((double) 1.0 / (2 + exp(discount_factor) + exp(-discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + pow(discount_factor, -2)));
	der_trans_dfact = ((double) 2.0 / (pow(discount_factor, -1) + 2 * discount_factor + pow(discount_factor, 3)));
#endif

	//Step 0: Initialize sizes and values
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
			{
					for (int j = 0; j < CONST_NUM_COVARIATES; j++)
					{
						household[i].time_step[t].evac_der_exante_value[u][j] = 0.0;
						household[i].time_step[t].wait_der_exante_value[u][j] = 0.0;
						household[i].time_step[t].stay_der_exante_value[u][j] = 0.0;
					}

					household[i].time_step[t].dfact_der_exante_value[u] = 0.0;
				}
			}
		
	}

	//Step 3: For every household, update the ex ante value derivatives (Backward Induction)
	double temp_e, temp_w, constant, temp_der, temp_der1, temp_der2;
	for (int i = 0; i < household.size(); i++) //For every household
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++) //For every time step
		{
			for (int u = CONST_NUM_TIME_STEPS - 1; u > t; u--) //Scan previous time steps to compute ex ante value function derivatives
			{
				
					if (u == CONST_NUM_TIME_STEPS - 1)
					{
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
					    temp_e += evac_estimate[CONST_NUM_SD_COVARIATES] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += stay_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						constant = (double) 1.0 / (temp_e + temp_w);

						//if (i == 0 && s == 0)
						//	cout << "DEBUG: " << constant << endl;

						//Derivatives wrt to evac coefficients 
						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							household[i].time_step[t].evac_der_exante_value[u][j] = constant * (household[i].time_step[t].covariate[j] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES] = constant * (household[i].time_step[t].order_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 1] = constant * (household[i].time_step[t].manorder_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 2] = constant * (household[i].time_step[t].volorder_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 3] = constant * (household[i].time_step[t].distlandfall_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 4] = constant * (household[i].time_step[t].intensity_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 5] = constant * (household[i].time_step[t].floodlevel_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 6] = constant * (household[i].time_step[t].distance_forecast[u] * temp_e);
						household[i].time_step[t].evac_der_exante_value[u][CONST_NUM_SD_COVARIATES + 7] = constant * (household[i].time_step[t].prob74_forecast[u] * temp_e);
					
						//if (i == 0 && s == 0 && t == 0)
						//	cout << "BINGO: " << u << " " << household[i].time_step[t].evac_der_exante_value[u][s][0];

						//Derivatives wrt to stay coefficients 
						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							household[i].time_step[t].stay_der_exante_value[u][j] = constant * (household[i].time_step[t].covariate[j] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES] = constant * (household[i].time_step[t].order_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 1] = constant * (household[i].time_step[t].manorder_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 2] = constant * (household[i].time_step[t].volorder_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 3] = constant * (household[i].time_step[t].distlandfall_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 4] = constant * (household[i].time_step[t].intensity_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 5] = constant * (household[i].time_step[t].floodlevel_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 6] = constant * (household[i].time_step[t].distance_forecast[u] * temp_w);
						household[i].time_step[t].stay_der_exante_value[u][CONST_NUM_SD_COVARIATES + 7] = constant * (household[i].time_step[t].prob74_forecast[u] * temp_w);

						
					}
					else
					{
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

						
						temp_w +=household[i].time_step[t].exante_value[u + 1];
						temp_w = trans_dfact * temp_w;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += wait_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES] * household[i].time_step[t].order_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						constant = (double) 1.0 / (temp_e + temp_w);

						//Derivatives wrt to evac coefficients
						temp_der = 0.0;
						for (int j = 0; j < CONST_NUM_COVARIATES; j++)
						{
							temp_der += household[i].time_step[t].evac_der_exante_value[u + 1][j];
							temp_der = trans_dfact * temp_der;
							
							if (j < CONST_NUM_SD_COVARIATES)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].covariate[j]+ temp_w * temp_der );

							if (j == CONST_NUM_SD_COVARIATES)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].order_forecast[u] + temp_w * temp_der);
							
							if (j == CONST_NUM_SD_COVARIATES + 1)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].manorder_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 2)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].volorder_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 3)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].distlandfall_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 4)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].intensity_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 5)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].floodlevel_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 6)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].distance_forecast[u] + temp_w * temp_der);

							if (j == CONST_NUM_SD_COVARIATES + 7)
								household[i].time_step[t].evac_der_exante_value[u][j] = constant * (temp_e * household[i].time_step[t].prob74_forecast[u] + temp_w * temp_der);

				
						}

						//Derivatives wrt to wait coefficients
						temp_der = 0.0;
						for (int j = 0; j < CONST_NUM_COVARIATES; j++)
						{
							temp_der += household[i].time_step[t].wait_der_exante_value[u + 1][j];
							temp_der = trans_dfact * temp_der;
							if (j < CONST_NUM_SD_COVARIATES)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].covariate[j] + temp_der));

							if (j == CONST_NUM_SD_COVARIATES)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].order_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 1)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].manorder_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 2)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].volorder_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 3)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].distlandfall_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 4)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].intensity_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 5)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].floodlevel_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 6)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].distance_forecast[u] + temp_der));
							
							if (j == CONST_NUM_SD_COVARIATES + 7)
								household[i].time_step[t].wait_der_exante_value[u][j] = constant * (temp_w * (household[i].time_step[t].prob74_forecast[u] + temp_der ));

	
						}

						//Derivatives wrt to stay coefficients
						temp_der = 0.0;
						for (int j = 0; j < CONST_NUM_COVARIATES; j++)
						{
							temp_der += household[i].time_step[t].stay_der_exante_value[u + 1][j];
						temp_der = trans_dfact * temp_der;	
				        household[i].time_step[t].stay_der_exante_value[u][j] = constant * (temp_w * temp_der);
						}

						//Derivatives wrt to discount factor 
						
						temp_der = 0.0;
						temp_der1 = 0.0;
						temp_der2 = 0.0;

						
						temp_der1 += household[i].time_step[t].dfact_der_exante_value[u + 1];
						temp_der1 = trans_dfact * temp_der1;

						
						temp_der2 +=household[i].time_step[t].exante_value[u + 1];
						temp_der2 = der_trans_dfact * temp_der2;
						household[i].time_step[t].dfact_der_exante_value[u] = constant * temp_w * (temp_der1 + temp_der2);
					}
				}
			}
		}
	
}

//Derivatives of conditional value functions
void calculate_cv_partials(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	double constant;
	double trans_dfact, der_trans_dfact;//Transformed discount factor and derivative of transformed discount factor
	double temp_der1, temp_der2;

#ifdef CONST_NO_TRANSFORM
	trans_dfact = discount_factor;
	der_trans_dfact = 1;
#endif

#ifdef CONST_LOGIT_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + exp(-discount_factor)));
	der_trans_dfact = ((double) 1.0 / (2 + exp(discount_factor) + exp(-discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + pow(discount_factor, -2)));
	der_trans_dfact = ((double) 2.0 / (pow(discount_factor, -1) + 2 * discount_factor + pow(discount_factor, 3)));
#endif

	//Compute partial derivatives of cv_evac 
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = CONST_NUM_TIME_STEPS - 1; t >= 0; t--)
		{
			//Initialize sizes
			household[i].time_step[t].evac_der_cv_evac.clear();
			household[i].time_step[t].wait_der_cv_evac.clear();
			household[i].time_step[t].stay_der_cv_evac.clear();

			household[i].time_step[t].evac_der_cv_evac.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].wait_der_cv_evac.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].stay_der_cv_evac.resize(CONST_NUM_COVARIATES, 0.0);

			household[i].time_step[t].dfact_der_cv_evac = 0.0;

			//Update derivatives
			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				household[i].time_step[t].evac_der_cv_evac[j] = household[i].time_step[t].covariate[j];
		}
	}

	//Compute partial derivatives of cv_wait - For every household, update the derivatives using backward induction
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = CONST_NUM_TIME_STEPS - 1; t >= 0; t--)
		{
			//Initialize sizes
			household[i].time_step[t].evac_der_cv_wait.clear();
			household[i].time_step[t].wait_der_cv_wait.clear();
			household[i].time_step[t].stay_der_cv_wait.clear();

			household[i].time_step[t].evac_der_cv_wait.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].wait_der_cv_wait.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].stay_der_cv_wait.resize(CONST_NUM_COVARIATES, 0.0);

			household[i].time_step[t].dfact_der_cv_wait = 0.0;

			//Update derivatives
			if (t == CONST_NUM_TIME_STEPS - 1)//If it is the last time period
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
					household[i].time_step[t].stay_der_cv_wait[j] = household[i].time_step[t].covariate[j];
			}
			else
			{
				//Partial derivatives wrt to evac coefficients
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				{
					
						/*if (i == 0 && t == 14 && j == 0)
						{
						cout << " CHECK INTENSITY: " << household[i].time_step[t].intensity_forecast[t + 1][s];
						cout << " CHECK DER: " << household[i].time_step[t].evac_der_exante_value[t + 1][s][j] << endl;
						}*/
						household[i].time_step[t].evac_der_cv_wait[j] = trans_dfact * ( household[i].time_step[t].evac_der_exante_value[t + 1][j]);
					}
				}

				//Partial derivatives wrt to wait coefficients
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				{
					household[i].time_step[t].wait_der_cv_wait[j] += trans_dfact * (household[i].time_step[t].wait_der_exante_value[t + 1][j]);
					household[i].time_step[t].wait_der_cv_wait[j] += household[i].time_step[t].covariate[j];
				}

				//Partial derivatives wrt to stay coefficients
				for (int j = 0; j < CONST_NUM_COVARIATES; j++) {
					household[i].time_step[t].stay_der_cv_wait[j] += trans_dfact * (household[i].time_step[t].stay_der_exante_value[t + 1][j]);
				}
				//Partial derivatives wrt to discount factor
				temp_der1 = 0.0;
				temp_der2 = 0.0;
				
				temp_der1 += household[i].time_step[t].dfact_der_exante_value[t + 1];
				temp_der1 = trans_dfact * temp_der1;

				
				temp_der2 +=  household[i].time_step[t].exante_value[t + 1];
				temp_der2 = der_trans_dfact * temp_der2;

				household[i].time_step[t].dfact_der_cv_wait = temp_der1 + temp_der2;
			}
		}
	
}

void calculate_probabilities(vector<households> &household)
{
	//Estimate the probabilities of waiting in each time period
	for (int i = 0; i < household.size(); i++)
	{
		//Re-initialize all values to zeros
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[i].time_step[t].prob_wait = 0.0;

		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].prob_wait = (double)exp(household[i].time_step[t].cv_wait) / (exp(household[i].time_step[t].cv_evac) + exp(household[i].time_step[t].cv_wait));
			//if (household[i].time_step[t].action == 1) //We can terminate early but we are commenting this for the macroscopic validation
			//	break;
		}
	}
}

void calculate_elasticity_probabilities(vector<households> &household)
{
	//Estimate the probabilities of waiting in each time period
	for (int i = 0; i < household.size(); i++)
	{
		//Re-initialize all values to zeros
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].prob_wait = 0.0;
			household[i].time_step[t].prob_wait = (double)exp(household[i].time_step[t].cv_wait) / (exp(household[i].time_step[t].cv_evac) + exp(household[i].time_step[t].cv_wait));
		}
	}
}

void calculate_prob_partials(vector<households> &household)
{
	double constant;

	//For every household, update the derivatives of probabilities
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			//Initialize sizes
			household[i].time_step[t].evac_der_prob_wait.clear();
			household[i].time_step[t].wait_der_prob_wait.clear();
			household[i].time_step[t].stay_der_prob_wait.clear();

			household[i].time_step[t].evac_der_prob_wait.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].wait_der_prob_wait.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].stay_der_prob_wait.resize(CONST_NUM_COVARIATES, 0.0);

			//household[i].time_step[t].dfact_der_prob_wait = 0.0;

			//Calculate probabilities
			constant = pow((1 + exp(household[i].time_step[t].cv_evac - household[i].time_step[t].cv_wait)), 2);
			constant = (double)-1.0 / constant;
			constant = (exp(household[i].time_step[t].cv_evac - household[i].time_step[t].cv_wait)) * constant;

			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			{
				household[i].time_step[t].evac_der_prob_wait[j] = constant * (household[i].time_step[t].evac_der_cv_evac[j] - household[i].time_step[t].evac_der_cv_wait[j]);
				household[i].time_step[t].wait_der_prob_wait[j] = constant * (household[i].time_step[t].wait_der_cv_evac[j] - household[i].time_step[t].wait_der_cv_wait[j]);
				household[i].time_step[t].stay_der_prob_wait[j] = constant * (household[i].time_step[t].stay_der_cv_evac[j] - household[i].time_step[t].stay_der_cv_wait[j]);
			}

			household[i].time_step[t].dfact_der_prob_wait = constant * (household[i].time_step[t].dfact_der_cv_evac - household[i].time_step[t].dfact_der_cv_wait);
		}
	}
}

void calculate_log_likelihood(vector<households> &household)
{
	//Estimate the probabilities of the policy and log-likelihood
	for (int i = 0; i < household.size(); i++)
	{
		household[i].prob_policy = 1;

		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			if (household[i].time_step[t].action == 1)
			{
				household[i].prob_policy = household[i].prob_policy * (1.0 - household[i].time_step[t].prob_wait);
				break;
			}
			else
				household[i].prob_policy = household[i].prob_policy * household[i].time_step[t].prob_wait;
		}
		household[i].log_likelihood = log(household[i].prob_policy);
	}
}

void calculate_likelihood_partials(vector<households> &household)
{
	//Re-initialize sizes of probabilities
	for (int i = 0; i < household.size(); i++)
	{
		household[i].evac_der_lnprob.clear();
		household[i].wait_der_lnprob.clear();
		household[i].stay_der_lnprob.clear();

		household[i].evac_der_lnprob.resize(CONST_NUM_COVARIATES, 0.0);
		household[i].wait_der_lnprob.resize(CONST_NUM_COVARIATES, 0.0);
		household[i].stay_der_lnprob.resize(CONST_NUM_COVARIATES, 0.0);

		household[i].dfact_der_lnprob = 0.0;
	}

	//Compute derivatives wrt coefficients
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			if (household[i].time_step[t].action == 1)
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				{
					household[i].evac_der_lnprob[j] += ((double) 1.0 / (1 - household[i].time_step[t].prob_wait)) * (-household[i].time_step[t].evac_der_prob_wait[j]);
					household[i].wait_der_lnprob[j] += ((double) 1.0 / (1 - household[i].time_step[t].prob_wait)) * (-household[i].time_step[t].wait_der_prob_wait[j]);
					household[i].stay_der_lnprob[j] += ((double) 1.0 / (1 - household[i].time_step[t].prob_wait)) * (-household[i].time_step[t].stay_der_prob_wait[j]);
				}

				household[i].dfact_der_lnprob += ((double) 1.0 / (1 - household[i].time_step[t].prob_wait)) * (-household[i].time_step[t].dfact_der_prob_wait);
				break;
			}
			else
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				{
					household[i].evac_der_lnprob[j] += ((double) 1.0 / household[i].time_step[t].prob_wait) * household[i].time_step[t].evac_der_prob_wait[j];
					household[i].wait_der_lnprob[j] += ((double) 1.0 / household[i].time_step[t].prob_wait) * household[i].time_step[t].wait_der_prob_wait[j];
					household[i].stay_der_lnprob[j] += ((double) 1.0 / household[i].time_step[t].prob_wait) * household[i].time_step[t].stay_der_prob_wait[j];
				}

				household[i].dfact_der_lnprob += ((double) 1.0 / household[i].time_step[t].prob_wait) * household[i].time_step[t].dfact_der_prob_wait;
			}
		}
	}
}

void calculate_score(vector<households> &household, vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include)
{
	for (int i = 0; i < household.size(); i++)
	{
		household[i].score.clear();

		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			if (evac_include[j] == 1)
				household[i].score.push_back(household[i].evac_der_lnprob[j]);

		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			if (wait_include[j] == 1)
				household[i].score.push_back(household[i].wait_der_lnprob[j]);

		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			if (stay_include[j] == 1)
				household[i].score.push_back(household[i].stay_der_lnprob[j]);

#ifdef CONST_OPTIMIZE_DFACT
		household[i].score.push_back(household[i].dfact_der_lnprob);
#endif
	}
}

void calculate_outer_prod(vector<households> &household, int num_variables)
{
	for (int i = 0; i < household.size(); i++)
	{
		//Define outer product for the household
		household[i].outer_product = new double*[num_variables];
		for (int n = 0; n < num_variables; n++)
			household[i].outer_product[n] = new double[num_variables];

		for (int n = 0; n < num_variables; n++)
			for (int m = 0; m < num_variables; m++)
				household[i].outer_product[n][m] = household[i].score[n] * household[i].score[m];

		//Print outer_product matrices
		/*
		cout << setprecision(12);
		for (int i = 111; i < 112; i++)
		{
		cout << "Printing information for household " << i + 1 << endl;
		display_vector(household[i].score, "Gradient");
		display_matrix(household[i].outer_product, num_variables, "Outer Product");
		}*/
	}
}

void calculate_direction(vector<households> &household, vector<double> &temp_direction, vector<double> &standard_error, double &gap, int num_variables)
{
	MatrixXd outer_prod(num_variables, num_variables);//avareage of outer product (B matrix in Train's notation)
	MatrixXd identity(num_variables, num_variables);//define identity matrix for more accurate matrix inversion
	MatrixXd outer_prod_inverse(num_variables, num_variables);
	VectorXd gradient(num_variables);//Gradient of the objective wrt the estimates
	VectorXd direction(num_variables);

	temp_direction.clear();

	//Construct the gradient
	for (int n = 0; n < num_variables; n++)
		gradient(n) = 0.0;

	for (int i = 0; i < household.size(); i++)
		for (int n = 0; n < num_variables; n++)
			gradient(n) += household[i].score[n];

	for (int n = 0; n < num_variables; n++)
		gradient(n) = (double)gradient(n) / household.size();

	//Compute the average of the outer products (which gives the apporximate hessian)
	for (int n = 0; n < num_variables; n++)
		for (int m = 0; m < num_variables; m++)
			outer_prod(n, m) = 0.0;

	//Take the average of household outer products
	for (int i = 0; i < household.size(); i++)
		for (int n = 0; n < num_variables; n++)
			for (int m = 0; m < num_variables; m++)
				outer_prod(n, m) += household[i].outer_product[n][m];

	for (int n = 0; n < num_variables; n++)
		for (int m = 0; m < num_variables; m++)
			outer_prod(n, m) = (double)outer_prod(n, m) / household.size();

	//cout << endl << "Outer product: " << endl;
	//cout << outer_prod;

	//Define identity matrix
	for (int n = 0; n < num_variables; n++)
		for (int m = 0; m < num_variables; m++)
			identity(n, m) = 0;

	for (int n = 0; n < num_variables; n++)
		identity(n, n) = 1;

	//outer_prod_inverse = outer_prod.inverse();

	FullPivLU<MatrixXd> operation(outer_prod);
	outer_prod_inverse = operation.inverse();

	//outer_prod_inverse = outer_prod.inverse();

	cout << setprecision(12);

#if DEBUG
	//cout << endl << endl << "The approximation of the Hessian is: " << endl << outer_prod << endl;
	//cout << endl << "Approximate Hessian Inverse: " << endl;
	//cout << outer_prod_inverse;
	//cout << endl << endl << "Checking if inverse is correct: " << endl << outer_prod * outer_prod_inverse << endl;
	cout << " Hessian invertible? " << operation.isInvertible();
#endif

	for (int n = 0; n < num_variables; n++)
		standard_error[n] = sqrt((double)outer_prod_inverse(n, n) / household.size());

	direction = outer_prod.inverse() * gradient;

	gap = gradient.dot(direction);
	cout << " Gap = " << gap << endl;

	//cout << endl << "Gradient: " << endl;
	//cout << gradient << endl;

	//cout << endl << "Direction: " << endl;
	//cout << direction << endl;

	//Copy the direction to the temp_direction vector
	for (int i = 0; i < num_variables; i++)
		temp_direction.push_back(direction(i));
}

void estimate_parameters(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, vector<double> &standard_error, double &discount_factor)
{
	cout << "Begining estimation.." << endl;
	double objective, old_objective, step_size;
	vector<double> direction;

	//--------Step 0: Determine the number of variables----------
	int num_variables = 0;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (evac_include[j] == 1)
			num_variables++;
		if (wait_include[j] == 1)
			num_variables++;
		if (stay_include[j] == 1)
			num_variables++;
	}
#ifdef CONST_OPTIMIZE_DFACT
	num_variables++;//Accounts for discount factor
#endif

					//------Step 1: Initialize the variables to be optimized--------
	VectorXd lambda(num_variables); //Values of the estimates
	VectorXd old_lambda(num_variables); //Temporary values of the estimates

	standard_error.clear();
	standard_error.resize(num_variables, 0.0);//Standard errors of estimates

											  //Intitalize lambda
	int index = 0;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (evac_include[j] == 1)
		{
			lambda(index) = evac_estimate[j];
			index++;
		}
	}

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (wait_include[j] == 1)
		{
			lambda(index) = wait_estimate[j];
			index++;
		}
	}

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (stay_include[j] == 1)
		{
			lambda(index) = stay_estimate[j];
			index++;
		}
	}
#ifdef CONST_OPTIMIZE_DFACT
	lambda(index) = discount_factor;
#endif

	//------------Step 3: Perform gradient descent iterations--------------
	double gap = CONST_INFTY;
	for (int iter = 0; iter < CONST_MAX_ITERATIONS; iter++)
	{
		cout << "Iteration: " << iter;
		cout << " Discount factor: " << discount_factor;
		//Solve the dynamic program, compute probabilities and associated derivatives
		calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
		calculate_ev_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
		calculate_cv_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
		calculate_probabilities(household);
		calculate_prob_partials(household);
		calculate_log_likelihood(household);
		calculate_likelihood_partials(household);

		/*household[9].display_value_functions();
		household[9].display_probabilities();
		household[9].display_derivatives();
		household[9].display_prob_partials();*/

		objective = 0.0;
		for (int i = 0; i < household.size(); i++)
			objective += household[i].log_likelihood;
		cout << " Objective = " << objective;

		/*cout << "Printing lambda (estimates): " << endl;
		cout << lambda << endl;*/

		calculate_score(household, evac_include, wait_include, stay_include);
		calculate_outer_prod(household, num_variables);
		calculate_direction(household, direction, standard_error, gap, num_variables);

		step_size = 2;
		old_objective = objective;
		old_lambda = lambda;

		//Half step size in each iteration until the objective increases
		//cout << endl << "Computing step size.." << "Current objective: " << objective << endl;
		while (objective < old_objective + CONST_LARGE_EPSILON && step_size > CONST_EPSILON)
		{
			//Update lambda's
			step_size = (double)step_size / 2.0;
			//cout << "Step size: " << step_size;

			for (int n = 0; n < num_variables; n++)
				lambda(n) = old_lambda(n) + step_size * direction[n];

			//Update coefficients and discount factor
			index = 0;
			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			{
				if (evac_include[j] == 1)
				{
					evac_estimate[j] = lambda(index);
					index++;
				}
			}

			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			{
				if (wait_include[j] == 1)
				{
					wait_estimate[j] = lambda(index);
					index++;
				}
			}

			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
			{
				if (stay_include[j] == 1)
				{
					stay_estimate[j] = lambda(index);
					index++;
				}
			}
#ifdef CONST_OPTIMIZE_DFACT
			discount_factor = lambda(index);
#endif

			//Compute the new objective after moving to the new lambda
			calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
			calculate_probabilities(household);
			calculate_log_likelihood(household);

			objective = 0.0;
			for (int i = 0; i < household.size(); i++)
				objective += household[i].log_likelihood;
			//cout << " Objective: " << objective << endl;
		}

		if (gap <= CONST_CONVERGENCE)
			break;
	}

	//Compute standard errors at the final values of the estimates
	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_ev_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_cv_partials(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_probabilities(household);
	calculate_prob_partials(household);
	calculate_log_likelihood(household);
	calculate_likelihood_partials(household);

	objective = 0.0;
	for (int i = 0; i < household.size(); i++)
		objective += household[i].log_likelihood;
	cout << endl << "Final Objective = " << objective;

	calculate_score(household, evac_include, wait_include, stay_include);
	calculate_outer_prod(household, num_variables);
	calculate_direction(household, direction, standard_error, gap, num_variables);

	cout << endl << "Displaying standard errors.." << endl;
	for (int n = 0; n < num_variables; n++)
		cout << standard_error[n] << endl;

	//Display Estimates
	/*cout << endl << "Displaying final estimates: " << endl;
	cout << lambda << endl;*/
}

void display_estimation_results(vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, vector<double> &standard_error, vector<string> &covariate_name, double &discount_factor)
{
	vector<double> evac_zvalue(CONST_NUM_COVARIATES, 0);
	vector<double> wait_zvalue(CONST_NUM_COVARIATES, 0);
	vector<double> stay_zvalue(CONST_NUM_COVARIATES, 0);

	//Populate z values
	int index = 0;
	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (evac_include[j] == 1)
		{
			evac_zvalue[j] = evac_estimate[j] / standard_error[index];
			index++;
		}
	}

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (wait_include[j] == 1)
		{
			wait_zvalue[j] = wait_estimate[j] / standard_error[index];
			index++;
		}
	}

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		if (stay_include[j] == 1)
		{
			stay_zvalue[j] = stay_estimate[j] / standard_error[index];
			index++;
		}
	}

	cout << endl << "Displaying estimation results.." << endl;
	cout << setw(CONST_EXTRA_WIDE_SPACING) << " ";
	cout << setw(CONST_WIDE_SPACING) << "Evacuate";
	cout << setw(CONST_WIDE_SPACING) << "z-value";
	cout << setw(CONST_WIDE_SPACING) << "Wait";
	cout << setw(CONST_WIDE_SPACING) << "z-value";
	cout << setw(CONST_WIDE_SPACING) << "Stay";
	cout << setw(CONST_WIDE_SPACING) << "z-value" << endl;

	cout << setprecision(7);

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	{
		cout << setw(CONST_EXTRA_WIDE_SPACING) << covariate_name[j];

		if (evac_include[j] == 1)
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << evac_estimate[j];
			cout << setw(CONST_EXTRA_WIDE_SPACING) << evac_zvalue[j];
		}
		else
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
		}
		if (wait_include[j] == 1)
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << wait_estimate[j];
			cout << setw(CONST_EXTRA_WIDE_SPACING) << wait_zvalue[j];
		}
		else
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
		}

		if (stay_include[j] == 1)
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << stay_estimate[j];
			cout << setw(CONST_EXTRA_WIDE_SPACING) << stay_zvalue[j];
		}
		else
		{
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
			cout << setw(CONST_EXTRA_WIDE_SPACING) << "N/A";
		}

		cout << endl;
	}
	cout << "Estimate of unconstrained discount factor: " << discount_factor << endl;

#ifdef CONST_OPTIMIZE_DFACT
#ifdef CONST_NO_TRANSFORM
	cout << "z-value of discount factor: " << (double)discount_factor / standard_error[index] << endl;
#endif 

#ifdef CONST_LOGIT_TRANSFORM
	cout << "Estimate of actual discount factor: " << ((double)exp(discount_factor) / (1 + exp(discount_factor))) << endl;
	cout << "The z-value needs to be estimated using the delta method (or re-run the model by fixing it at the current value)" << endl;
#endif

#ifdef CONST_SQUARE_TRANSFORM
	cout << "Estimate of actual discount factor: " << ((double)pow(discount_factor, 2) / (1 + pow(discount_factor, 2))) << endl;
	cout << "The z-value needs to be estimated using the delta method (or re-run the model by fixing it at the current value)" << endl;
#endif
#endif

	//Print results for fixing variables and computing the standard error of the discount factor
	/*for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	cout << "evac_estimate[" << j << "] = " << evac_estimate[j] << ";" << endl;

	for (int j = 0; j < CONST_NUM_COVARIATES; j++)
	cout << "wait_estimate[" << j << "] = " << wait_estimate[j] << ";" << endl;*/
}

void calculate_elasticities(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, double &discount_factor)
{
	//Step 1: Set discount factors
	double trans_dfact;
#ifdef CONST_NO_TRANSFORM
	trans_dfact = discount_factor;
#endif

#ifdef CONST_LOGIT_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + exp(-discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + pow(discount_factor, -2)));
#endif

	//Initialize size  
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			household[i].time_step[t].covar_der_exante_value = new double*[CONST_NUM_TIME_STEPS];
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
			{
			     household[i].time_step[t].covar_der_exante_value[u] = new double[CONST_NUM_COVARIATES];
			}
		}
	}

	//Initialize values 
	for (int i = 0; i < household.size(); i++)
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			for (int u = 0; u < CONST_NUM_TIME_STEPS; u++)
					for (int j = 0; j < CONST_NUM_COVARIATES; j++)
						household[i].time_step[t].covar_der_exante_value[u][j] = 0.0;

	//Find derivatives of exante value functions using backward induction
	double constant, temp_e, temp_w, temp_der;
	for (int i = 0; i < household.size(); i++) //For every household
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++) //For every time step
		{
			for (int u = CONST_NUM_TIME_STEPS - 1; u > t; u--) //Scan previous time steps to compute ex ante value function derivatives
			{
				
					if (u == CONST_NUM_TIME_STEPS - 1)
					{
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += stay_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += stay_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						constant = (double) 1.0 / (temp_e + temp_w);

						//Derivatives wrt to sd covariates (others are anyway initialized to zeros)
						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							household[i].time_step[t].covar_der_exante_value[u][j] = constant * (temp_e * evac_estimate[j] + temp_w * stay_estimate[j]);
					}
					else
					{
						temp_e = 0.0;
						temp_w = 0.0;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_e += evac_estimate[j] * household[i].time_step[u].covariate[j];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_e += evac_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_e = exp(temp_e);

						
						temp_w += household[i].time_step[t].prob74_forecast[u + 1] * household[i].time_step[t].exante_value[u + 1];
						temp_w = trans_dfact * temp_w;

						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
							temp_w += wait_estimate[j] * household[i].time_step[u].covariate[j];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES ] * household[i].time_step[t].order_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 1] * household[i].time_step[t].manorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 2] * household[i].time_step[t].volorder_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 3] * household[i].time_step[t].distlandfall_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 4] * household[i].time_step[t].intensity_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 5] * household[i].time_step[t].floodlevel_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 6] * household[i].time_step[t].distance_forecast[u];
						temp_w += wait_estimate[CONST_NUM_SD_COVARIATES + 7] * household[i].time_step[t].prob74_forecast[u];
						temp_w = exp(temp_w);

						constant = (double) 1.0 / (temp_e + temp_w);

						//Derivatives wrt to sd covariates (others are anyway initialized to zeros)
						for (int j = 0; j < CONST_NUM_SD_COVARIATES; j++)
						{
							temp_der = 0.0;
					
							temp_der += household[i].time_step[t].prob74_forecast[u + 1] * household[i].time_step[t].covar_der_exante_value[u + 1][j];
							temp_der = trans_dfact * temp_der;

							household[i].time_step[t].covar_der_exante_value[u][j] = constant * (temp_e * evac_estimate[j] + temp_w * (wait_estimate[j] + temp_der));
						}
					}
				}
		}
	}

	//Derivatives of conditional value functions
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = CONST_NUM_TIME_STEPS - 1; t >= 0; t--)
		{
			household[i].time_step[t].covar_der_cv_evac.resize(CONST_NUM_COVARIATES, 0.0);
			household[i].time_step[t].covar_der_cv_wait.resize(CONST_NUM_COVARIATES, 0.0);

			//Update derivatives of cv_evac
			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				household[i].time_step[t].covar_der_cv_evac[j] = evac_estimate[j];

			//Update derivatives of cv_wait
			if (t == CONST_NUM_TIME_STEPS - 1)
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
					household[i].time_step[t].covar_der_cv_wait[j] = stay_estimate[j];
			}
			else
			{
				for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				{
					
					household[i].time_step[t].covar_der_cv_wait[j] += trans_dfact * (household[i].time_step[t].prob74_forecast[t + 1] * household[i].time_step[t].covar_der_exante_value[t + 1][j]);
					household[i].time_step[t].covar_der_cv_wait[j] += wait_estimate[j];
				}
			}
		}
	}

	//Derivatives of probabilities
	for (int i = 0; i < household.size(); i++)
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			//Initialize sizes
			household[i].time_step[t].covar_der_prob.resize(CONST_NUM_COVARIATES, 0.0);

			//Calculate derivatives
			constant = pow((1 + exp(household[i].time_step[t].cv_evac - household[i].time_step[t].cv_wait)), 2);
			constant = (double)-1.0 / constant;
			constant = (exp(household[i].time_step[t].cv_evac - household[i].time_step[t].cv_wait)) * constant;

			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				household[i].time_step[t].covar_der_prob[j] = constant * (household[i].time_step[t].covar_der_cv_evac[j] - household[i].time_step[t].covar_der_cv_wait[j]);
		}
	}
}

void display_value_functions(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	//Step 1: Set discount factors
	double trans_dfact;
#ifdef CONST_NO_TRANSFORM
	trans_dfact = discount_factor;
#endif

#ifdef CONST_LOGIT_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + exp(-discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
	trans_dfact = ((double) 1.0 / (1 + pow(discount_factor, -2)));
#endif

	cout << "Printing value functions: " << endl;
	double temp = 0.0;
	for (int t = CONST_NUM_TIME_STEPS - 1; t >= 0; t--)
	{
		cout << t;
		if (t == CONST_NUM_TIME_STEPS - 1)
		{
			cout << "\t" << household[0].time_step[t].cv_evac;
			cout << "\t" << household[0].time_step[t].cv_wait;
			cout << "\t" << 0 << endl;
		}
		else
		{
			temp = 0.0;
			cout << "\t" << household[0].time_step[t].cv_evac;
			for (int j = 0; j < CONST_NUM_COVARIATES; j++)
				temp += wait_estimate[j] * household[0].time_step[t].covariate[j];
			cout << "\t" << household[0].time_step[t].cv_wait;
			cout << "\t" << temp;
			temp = 0.0;
			
			temp += trans_dfact * household[0].time_step[t].prob74_forecast[t + 1] * household[0].time_step[t].exante_value[t + 1];
			cout << "\t" << temp << endl;
		}
	}

	cout << "Printing derivatives of conditinal value functions wrt to covariates: " << endl;
	for (int t = CONST_NUM_TIME_STEPS - 1; t >= 0; t--)
	{
		cout << t;
		cout << "\t" << household[0].time_step[t].covar_der_cv_evac[2];
		if (t == CONST_NUM_TIME_STEPS - 1)
		{
			cout << "\t" << household[0].time_step[t].covar_der_cv_wait[2] << endl;
		}
		else
		{
			temp = 0.0;
			
			temp += trans_dfact * (household[0].time_step[t].prob74_forecast[t + 1] * household[0].time_step[t].covar_der_exante_value[t + 1][2]);
			cout << "\t" << household[0].time_step[t].covar_der_cv_wait[2] << "\t" << wait_estimate[2] << "\t" << temp << endl;
		}
	}
}

void output_finite_difference(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor, vector<string> &covariate_name, int covariate_index, int household_index)
{
	vector<double> prob_diff(CONST_NUM_TIME_STEPS, 0.0);//Difference in probabilities after and before for Prob(evac) not Prob(wait)

	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_elasticity_probabilities(household);
	//cout << "Printing data for estimating finite differences for " << covariate_name[covariate_index] << ": ";
	//household[household_index].display_probabilities();

	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		prob_diff[t] = household[household_index].time_step[t].prob_wait;

	//cout << "Current " << covariate_name[covariate_index] << " value: " << household[household_index].time_step[0].covariate[covariate_index] << endl;
	//cout << "Increase the covariate by 1: " << endl;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		household[household_index].time_step[t].covariate[covariate_index] += 1;
	//cout << "New " << covariate_name[covariate_index] << " value: " << household[household_index].time_step[0].covariate[covariate_index] << endl;

	//Recompute probabiltiies
	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_elasticity_probabilities(household);
	//household[household_index].display_probabilities();

	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		prob_diff[t] -= household[household_index].time_step[t].prob_wait;

	//Write different to dat file (We will print the derivatives wrt to probability of evacuation as reported in the journal paper)
	ofstream outfile("./Results/" + covariate_name[covariate_index] + "_" + to_string(household_index) + ".dat", std::ios::trunc);
	outfile.precision(5);
	if (outfile.is_open())
	{
		outfile << "Time \t Prob" << endl;
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		{
			outfile << t + 1 << "\t";
			outfile << prob_diff[t];
			outfile << endl;
		}
	}
	else
		cout << "ERROR: Could not create output dat file " << endl;

	//Reset things back to how it was
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		household[household_index].time_step[t].covariate[covariate_index] -= 1;

	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);
	calculate_elasticity_probabilities(household);
}

void output_parameter_values(vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, double &discount_factor)
{
	ofstream outfile("./initial_parameter_values.txt");
	outfile.precision(12);
	if (outfile.is_open())
	{
		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
		{
			if (evac_include[j] == 1)
				outfile << "1" << "\t" << j << "\t" << evac_estimate[j] << endl;
		}

		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
		{
			if (wait_include[j] == 1)
				outfile << "2" << "\t" << j << "\t" << wait_estimate[j] << endl;
		}

		for (int j = 0; j < CONST_NUM_COVARIATES; j++)
		{
			if (stay_include[j] == 1)
				outfile << "3" << "\t" << j << "\t" << stay_estimate[j] << endl;
		}

#ifdef CONST_NO_TRANSFORM
		outfile << "4" << "\t" << discount_factor;
#endif 

#ifdef CONST_LOGIT_TRANSFORM
		outfile << "4" << "\t" << ((double)exp(discount_factor) / (1 + exp(discount_factor)));
#endif

#ifdef CONST_SQUARE_TRANSFORM
		outfile << "4" << "\t" << ((double)pow(discount_factor, 2) / (1 + pow(discount_factor, 2)));
#endif
	}
	else
		cout << "ERROR: Could not create output parameter values file " << endl;

	outfile.close();
}

void set_initial_values(vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	string line;
	vector<string> token;
	ifstream file("./initial_parameter_values.txt");

	while (getline(file, line))
	{
		token = split_string(line.c_str(), '\t', ' ');

		//Store evac parameters
		if (atoi(token[0].c_str()) == 1)
			evac_estimate[atoi(token[1].c_str())] = atof(token[2].c_str());

		//Store wait parameters
		if (atoi(token[0].c_str()) == 2)
			wait_estimate[atoi(token[1].c_str())] = atof(token[2].c_str());

		//Store stay parameters
		if (atoi(token[0].c_str()) == 3)
			stay_estimate[atoi(token[1].c_str())] = atof(token[2].c_str());

		//Store discount factor 
		if (atoi(token[0].c_str()) == 4)
			discount_factor = atof(token[1].c_str());
	}
	file.close();
}

void predict_choice_probabilities(vector<households> &household, vector<double> &evac_estimate, vector<double> &wait_estimate, vector<double> &stay_estimate, double &discount_factor)
{
	//Should we require k-fold validation, minor modifications can be made to this function
	//Can we vizualize the prediction probability on the cumulative demand curve? 
	calculate_value_functions(household, evac_estimate, wait_estimate, stay_estimate, discount_factor);

	//Estimate choice probabilities at each time step
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		household[0].time_step[t].prob_wait = 0.0;

	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		household[0].time_step[t].prob_wait = (double)exp(household[0].time_step[t].cv_wait) / (exp(household[0].time_step[t].cv_evac) + exp(household[0].time_step[t].cv_wait));
		//if (household[0].time_step[t].action == 1)
		//	break;
	}

	//Estimate the probabilities of the observed policy
	household[0].prob_policy = 1;

	//Old function (exact probability)
	/*
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
	if (household[0].time_step[t].action == 1)
	{
	household[0].prob_policy = household[0].prob_policy * (1.0 - household[0].time_step[t].prob_wait);
	break;
	}
	else
	household[0].prob_policy = household[0].prob_policy * household[0].time_step[t].prob_wait;
	}*/

	//Print probability of wait
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		cout << t << " " << household[0].time_step[t].prob_wait << endl;

	//Calculating the probability that the household will evacuate
	double temp_stay_prob = 1.0;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
		temp_stay_prob = temp_stay_prob * household[0].time_step[t].prob_wait;

	//Add to the above probability, the probability of evacuating in adjacent time steps
	int evac_time_index = -1;
	for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
	{
		if (household[0].time_step[t].action == 1)
		{
			evac_time_index = t;
			break;
		}
	}

	ofstream outfile_prob("evac_prob_comparison.txt", std::ios::app);
	outfile_prob.precision(5);
	if (outfile_prob.is_open())
		outfile_prob << evac_time_index << " " << 1 - temp_stay_prob << endl;
	else
		cout << "ERROR: Could not create output dat file " << endl;
	outfile_prob.close();

	int index;
	double temp_prob = 1.0;
	double left, middle, right;
	if (evac_time_index != -1) //If the household evacuates
	{
		if (evac_time_index == 0)
			index = evac_time_index + 1;
		else if (evac_time_index == CONST_NUM_TIME_STEPS - 1)
			index = evac_time_index - 1;
		else
			index = evac_time_index;

		for (int t = 0; t < index - 1; t++)
			temp_prob = temp_prob * household[0].time_step[t].prob_wait;
		left = 1 - household[0].time_step[index - 1].prob_wait;
		middle = household[0].time_step[index - 1].prob_wait * (1 - household[0].time_step[index].prob_wait);
		right = household[0].time_step[index - 1].prob_wait * household[0].time_step[index].prob_wait * (1 - household[0].time_step[index + 1].prob_wait);
		household[0].prob_policy = temp_prob * (left + middle + right);
	}
	else //If the household stays
	{
		for (int t = 0; t < CONST_NUM_TIME_STEPS; t++)
			household[0].prob_policy = household[0].prob_policy * household[0].time_step[t].prob_wait;
	}

	household[0].log_likelihood = log(household[0].prob_policy);

	cout << "Predicted probability of the LOOCV household's sequence of actions: " << household[0].prob_policy << endl;
	cout << "The corresponding log likelihood value is: " << household[0].log_likelihood << endl;

	//Write the probabilities to a file
	ofstream outfile("predicted_probabilities.txt", std::ios::app);
	outfile.precision(5);
	if (outfile.is_open())
		outfile << household[0].prob_policy << endl;
	else
		cout << "ERROR: Could not create output dat file " << endl;
	outfile.close();
}

void validate_model(vector<int> &evac_include, vector<int> &wait_include, vector<int> &stay_include, vector<string> &covariate_name)
{
	/*vector<households> valid_household(CONST_NUM_HOUSEHOLDS_FLORENCE+ CONST_NUM_HOUSEHOLDS_MICHAEL);
	vector<double> standard_error;
	vector<households> loocv_household;

	vector<double> valid_evac_estimate(CONST_NUM_COVARIATES, 0.0);//Evacuate
	vector<double> valid_wait_estimate(CONST_NUM_COVARIATES, 0.0);//Defer decision to next time period
	vector<double> valid_stay_estimate(CONST_NUM_COVARIATES, 0.0);//Stay indefinitely
	double valid_discount_factor = 1;//Unconstrained discount variable. The actual discount_variable used is a transformed one

	vector<households> household(CONST_NUM_HOUSEHOLDS_FLORENCE + CONST_NUM_HOUSEHOLDS_MICHAEL);
	set_household_data(household, covariate_name);
	set_intensity_forecasts(household);
	set_distance_forecasts(household);

	//int a;
	double temp_demand;
	vector<double> expected_demand(CONST_NUM_TIME_STEPS + 1, 0.0);
	for (int i = 0; i < household.size(); i++)
	{
		valid_household.clear();
		valid_household.resize(CONST_NUM_HOUSEHOLDS_FLORENCE + CONST_NUM_HOUSEHOLDS_MICHAEL);

		valid_household = household;
		valid_household.erase(valid_household.begin() + i);

		cout << "Printing size of training data: " << valid_household.size() << " after removing household " << i + 1 << endl;
		//cin >> a;

		loocv_household.clear();
		loocv_household.push_back(household[i]);

		loocv_household[0].display_data(covariate_name); //Change loocv's id

		cout << "Information on other households: " << endl;
		for (int j = 0; j < valid_household.size(); j++)
		{
			if (j >= i)
				valid_household[j].id--;
			cout << valid_household[j].id << "\t";
			for (int k = 0; k < CONST_NUM_SD_COVARIATES; k++)
				cout << valid_household[j].time_step[0].covariate[k] << "\t";

			cout << endl;
		}

		valid_evac_estimate.clear();
		valid_wait_estimate.clear();
		valid_stay_estimate.clear();

		valid_evac_estimate.resize(CONST_NUM_COVARIATES, 0.0);
		valid_wait_estimate.resize(CONST_NUM_COVARIATES, 0.0);
		valid_stay_estimate.resize(CONST_NUM_COVARIATES, 0.0);
		valid_discount_factor = 1;

		display_vector(valid_evac_estimate, "Evacuate Parameters");
		display_vector(valid_wait_estimate, "Wait Parameters");
		display_vector(valid_stay_estimate, "Stay Parameters");

		//Estimate parameters by maximizing log likelihood
		estimate_parameters(valid_household, valid_evac_estimate, valid_wait_estimate, valid_stay_estimate, evac_include, wait_include, stay_include, standard_error, valid_discount_factor);

		//Display estimation results
		display_estimation_results(valid_evac_estimate, valid_wait_estimate, valid_stay_estimate, evac_include, wait_include, stay_include, standard_error, covariate_name, valid_discount_factor);

		//Estimate the preditive abilities of the model
		predict_choice_probabilities(loocv_household, valid_evac_estimate, valid_wait_estimate, valid_stay_estimate, valid_discount_factor);

		//Print aggregate demand information
		for (int t = 0; t < CONST_NUM_TIME_STEPS + 1; t++)
		{
			temp_demand = 1;
			for (int s = 0; s < t; s++)
				temp_demand = temp_demand * loocv_household[0].time_step[s].prob_wait;

			if (t != CONST_NUM_TIME_STEPS)
				temp_demand = temp_demand * (1 - loocv_household[0].time_step[t].prob_wait);

			expected_demand[t] += temp_demand;
		}
	}

	cout << endl << "Printing expected demand using LOOCV households: " << endl;
	for (int i = 0; i < expected_demand.size(); i++)
		cout << expected_demand[i] << endl;*/
}