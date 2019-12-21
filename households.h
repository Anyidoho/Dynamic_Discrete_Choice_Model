#ifndef households_h
#define households_h

//Forward declaration of classes
/*class states //Is initialized in the calculate_value_functions code
{
public:
int intensity;//varies from 0 to 5

//---------Value function variables (are populated in the function calculate_value_functions())------
double cv_evac;
double cv_wait;
vector<double> exante_value;//vector of exante values for all future time periods

//---Derivatives of the conditional values (are populated in the function calculate_cv_parials())----
vector<double> evac_der_cv_evac;//Derivative of the conditional value of evac wrt 'evac' coefficients
vector<double> wait_der_cv_evac;//Derivative of the conditional value of evac wrt 'wait' coefficients
vector<double> stay_der_cv_evac;//Derivative of the conditional value of evac wrt 'stay' coefficients
double iparam_der_cv_evac;//Derivative of the conditional value of evac wrt to the intensity paramter

vector<double> evac_der_cv_wait;//Derivative of the conditional value of wait wrt 'evac' coefficients
vector<double> wait_der_cv_wait;//Derivative of the conditional value of wait wrt 'wait' coefficients
vector<double> stay_der_cv_wait;//Derivative of the conditional value of wait wrt 'stay' coefficients
double iparam_der_cv_wait;//Derivative of the conditional value of wait wrt to the intensity paramter
};*/

class time_steps
{
public:
	//Variable data
	int action;//action taken by the household. Can be 0 or 1
	vector<double> covariate;//Independent variable data

							 //vector<states> state;//Length of the vector is CONST_NUM_STATES
							 //----------Value function variables (are populated in the function calculate_value_functions())------
	double cv_evac;
	double cv_wait;

	double* exante_value;//[u] - u is the future time step 
	double** evac_der_exante_value;//[u][s][j] - u is future time step, j is covariate index
	double** wait_der_exante_value;
	double** stay_der_exante_value;
	double* dfact_der_exante_value;

	//Derivatives of the conditional values (are populated in the function calculate_cv_parials())
	vector<double> evac_der_cv_evac;//Derivative of the conditional value of evac wrt 'evac' coefficients
	vector<double> wait_der_cv_evac;//Derivative of the conditional value of evac wrt 'wait' coefficients
	vector<double> stay_der_cv_evac;//Derivative of the conditional value of evac wrt 'stay' coefficients
	double dfact_der_cv_evac;//Derivative of the conditional value of evac wrt the discount factor

	vector<double> evac_der_cv_wait;//Derivative of the conditional value of wait wrt 'evac' coefficients
	vector<double> wait_der_cv_wait;//Derivative of the conditional value of wait wrt 'wait' coefficients
	vector<double> stay_der_cv_wait;//Derivative of the conditional value of wait wrt 'stay' coefficients
	double dfact_der_cv_wait;//Derivative of the conditional value of wait wrt the discount factor*/

							 //Intensity forecasts
	//double** intensity_forecast;//[i][j] - i is the future time step and j is the intensity (0 to 6)	

								// Forecast values
	vector<int> order_forecast;//order perveivedness for future time steps
	vector<int> manorder_forecast;//Mandatory order status for future time steps
	vector<int> volorder_forecast;//Voluntary order status for future time steps
	vector<double> distlandfall_forecast;//Distance from household to lanfall area for future time steps
	vector<int> intensity_forecast;// Intensity values for future time steps
	vector<double> floodlevel_forecast; // Potential level of flooding for future time steps
	vector<double> distance_forecast;//Distance from household to predicted location for future time steps
	vector<double> prob74_forecast;// Prob of wind speed greater than 74mph for future time steps
	

	

									 //------------Probability variables (are calculated in calculate_probabilities())--------------------
	double prob_wait;//probability of waiting

					 //Derivatives of probabilities (are calculated in calculate_prob_partials())
	vector<double> evac_der_prob_wait;//Derivative of the probability of wait wrt to 'evac' coefficients
	vector<double> wait_der_prob_wait;
	vector<double> stay_der_prob_wait;
	double dfact_der_prob_wait;

	//----------------------------------Elasticity variables----------------------------------------------
	double** covar_der_exante_value;//[u][s][j] - u is future time step, s is future state, j is covariate index

	vector<double> covar_der_cv_evac;//Derivative of the conditional value of evac wrt covariates
	vector<double> covar_der_cv_wait;//Derivative of the conditional value of wait wrt covariates

	vector<double> covar_der_prob;//This is the probability of wait
};

class households
{
	friend void set_household_data(vector<households> &, vector<string> &, string);
	//friend void set_intensity_forecasts(vector<households> &);
	friend void set_distance_forecasts(vector<households> &, string, string);
	friend void calculate_value_functions(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
	friend void calculate_ev_partials(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
	friend void calculate_cv_partials(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
	friend void calculate_probabilities(vector<households> &);
	friend void calculate_prob_partials(vector<households> &);
	friend void calculate_log_likelihood(vector<households> &);
	friend void calculate_likelihood_partials(vector<households> &);
	friend void calculate_score(vector<households> &, vector<int> &, vector<int> &, vector<int> &);
	friend void calculate_outer_prod(vector<households> &, int);
	friend void estimate_parameters(vector<households> &, vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, vector<double> &, double &);
	friend void display_estimation_results(vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, vector<double> &, vector<string> &, double &);
	friend void calculate_elasticity_probabilities(vector<households> &);
	friend void display_value_functions(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);

public:
	households();
	~households();
	void display_data(vector<string> &) const;
	void display_value_functions() const;
	void display_derivatives() const;
	void display_probabilities() const;
	void display_prob_partials() const;
	void display_elasticities(vector<string> &) const;
	void output_elasticities(vector<string> &, int) const;
	/*void display_cv_partials() const;

	void free_memory();*/

public:
	int id;//household id
	vector<time_steps> time_step;

	double latitude;
	double longitude;

	double prob_policy;//Probability of the policy (Sequence of actions) of the household
	double log_likelihood;//Contribution to the log-likelihood function

	vector<double> evac_der_lnprob;//Derivative of the log of probability of policy wrt to evac coefficients
	vector<double> wait_der_lnprob;
	vector<double> stay_der_lnprob;
	double dfact_der_lnprob;

	vector<double> score;//A collection of parital derivatives wrt all parameters that need to be optimized
	double** outer_product;
};

#endif