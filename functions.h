#ifndef functions_h
#define functions_h

//Forward declaration of classes
class households;

//Function Prototypes
vector<string> split_string(const char *, char, char);
double compute_distance(double, double, double, double);
double product_rule(vector<double> &, vector<double> &);
void set_household_data(vector<households> &, vector<string> &, string);
//void set_intensity_forecasts(vector<households> &);
void set_distance_forecasts(vector<households> &, string, string);
void set_include_indicators(vector<int> &, vector<int> &, vector<int> &, vector<string> &);
void calculate_value_functions(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
void calculate_ev_partials(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
void calculate_cv_partials(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
void calculate_probabilities(vector<households> &);
void calculate_prob_partials(vector<households> &);
void calculate_log_likelihood(vector<households> &);
void calculate_likelihood_partials(vector<households> &);
void calculate_score(vector<households> &, vector<int> &, vector<int> &, vector<int> &);
void calculate_outer_prod(vector<households> &, int);
void calculate_direction(vector<households> &, vector<double> &, vector<double> &, double &, int);
void estimate_parameters(vector<households> &, vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, vector<double> &, double &);
void display_estimation_results(vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, vector<double> &, vector<string> &, double &);
void calculate_elasticities(vector<households> &, vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, double &);
void calculate_elasticity_probabilities(vector<households> &);
void display_value_functions(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
void output_finite_difference(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &, vector<string> &, int, int);
void output_parameter_values(vector<double> &, vector<double> &, vector<double> &, vector<int> &, vector<int> &, vector<int> &, double &);
void set_initial_values(vector<double> &, vector<double> &, vector<double> &, double &);
void predict_choice_probabilities(vector<households> &, vector<double> &, vector<double> &, vector<double> &, double &);
void validate_model(vector<int> &, vector<int> &, vector<int> &, vector<string> &);

//Templates
template <class T> void display_matrix(T **matrix, int size, string message)
{
	cout << "Printing: " << message << endl;
	for (int i = 0; i<size; i++)
	{
		for (int j = 0; j<size; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

template <class T> void display_vector(T &vector, string message)
{
	cout << "Printing: " << message << endl;
	for (int i = 0; i<vector.size(); i++)
		cout << vector[i] << " ";
	cout << endl;
}
#endif
