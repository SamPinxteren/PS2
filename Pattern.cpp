#include "Pattern.h"

BigInt Pattern::g(int b, int e)
{
	if (b < e){
		int temp = e;
		e = b;
		b = temp;
	}

	BigInt f;
	auto search = g_dict.find(std::make_pair(b, e));
	if (search == g_dict.end()){
		for(int i = b + 1; i <= b + e; ++i){
			f *= i;
		}

		for(int i = 2; i <= e; ++i){
			f /= i;
		}
		g_dict[std::make_pair(b, e)] = f;
	} else {
		f = search->second;
	}

	return f;
}

double Pattern::C(std::vector<unsigned int> X)
{
	if (X.size() == 1) return 1.0;
	if (C_dict.find(X) != C_dict.end()) return C_dict[X];

	unsigned int edge_sum = 0;
	for (unsigned int i: X){
		edge_sum += i;
	}

	double* Ve = new double[edge_sum]();
	double* Vo = new double[edge_sum]();
	double* V;
	double* temp;
	Ve[0] = 1;
	Vo[0] = 1;

	unsigned int l = 0;
	for (unsigned int n = 1; n < X.size(); ++n){
		l += X[n - 1];
		BigInt norm_term = g(l, X[n]);

		if (n % 2){
			V = Vo;
			temp = Ve;
		} else {
			V = Ve;
			temp = Vo;
		}

		for (unsigned int j = 0; j < l; ++j){
			temp[j] = 0;
		}

		for (unsigned int j = 0; j < l; ++j){
			for (unsigned int p = j + 1; p <= l; ++p){
				for (unsigned int t = 0; t < X[n]; ++t){
					double f1 = V[j];
					BigInt f2 = g(j, t);
					f2 *= g(l - p, X[n] - t - 1);
					f2 /= norm_term;
					temp[p + t] += f1 * f2.to_double();
				}
			}
		}
	}

	double result = 0;
	for (unsigned int i = 0; i < edge_sum; ++i){
		result += temp[i];
	}

	delete [] Ve;
	delete [] Vo;

	C_dict[X] = result;
	return result;
}

double Pattern::occurs_probability()
{
	std::vector<unsigned int> X;
	for (auto const& e: symbol_counts){
		if (e.second == 0) return 0;
		X.push_back(e.second);
	}
	std::sort(X.begin(), X.end());
	return C(X);
}

Pattern::Pattern(std::vector<std::string> pattern_symbols, unsigned int verbosity):
	Pattern{pattern_symbols}
{
	verbose = verbosity;
}

Pattern::Pattern(std::vector<std::string> pattern_symbols):
	expected_value(0),
	variance(0),
	real_value(0),
	active_symbol(0),
	symbols{pattern_symbols}
{
	// Ensure input does not contain duplicates
	std::vector<std::string> unique_list = pattern_symbols;
	std::sort(unique_list.begin(), unique_list.end());
	auto u = std::unique(unique_list.begin(), unique_list.end());
	if (std::distance(unique_list.begin(), u) != pattern_symbols.size()){
		std::stringstream error_msg;
		error_msg << "Duplicates in pattern: ";
		for (auto i: pattern_symbols){
			error_msg << i << " ";
		}
		throw std::domain_error(error_msg.str());
	}

	// Initialize counters
	for (auto const& e: symbols){
		symbol_counts[e] = 0;
		total_symbol_counts[e] = 0;
	}
}

double Pattern::get_standard_deviation() const
{
	return sqrt(variance);
}

double Pattern::get_expected_value() const
{
	return expected_value;
}

unsigned int Pattern::get_support() const
{
	return real_value;
}

double Pattern::get_p_normal() const
{
	double z = (get_support() - 0.5 - get_expected_value()) / get_standard_deviation();
	return erfc(z / sqrt(2.0)) / 2.0;
}

double Pattern::get_p_exact() const
{
	double* Q = new double[P.size()+1]();
	Q[0] = 1;
	for (const double p: P){
		for (int i = P.size(); i >= 1; --i){
			Q[i] = Q[i]*(1-p) + Q[i-1]*p;
		}
		Q[0] = Q[0]*(1-p);
	}
	double p = 0;
	for (int i = get_support(); i <= P.size(); ++i){
		p += Q[i];
	}

	delete [] Q;
	return p;
}

double Pattern::get_p_poisson() const
{
	double lambda = get_expected_value();
	unsigned int val = get_support();
	unsigned int max_val = get_non_zero_sequences();
	bool reverse = true;
	if (max_val - val < val){
		val = max_val - val;
		lambda = max_val - lambda;
		reverse = false;
	}
	double s = 0;
	double f = 1;
	for (unsigned int i = 1; i < val; ++i){
		f *= i;
		s += pow(lambda, i) / f;
	}

	if (reverse){
		return 1 - (exp(-lambda) * (1 + s));
	} else {
		return exp(-lambda) * (1 + s);
	}
}

unsigned int Pattern::get_non_zero_sequences() const
{
	return P.size();
}

void Pattern::process()
{
	double occurs_prob = occurs_probability();
	double var = occurs_prob * (1.0 - occurs_prob);
	int occuring = (active_symbol == symbols.size());

	if (verbose >= 2){
		std::cout
			<< to_string()
			<< " occuring=" << occuring
			<< " var=" << var
			<< " p=" << occurs_prob
			<< std::endl;
	}

	expected_value += occurs_prob;
	real_value += occuring;
	variance += var;
	if (occurs_prob > 0){
		P.push_back(occurs_prob);
	}
}

void Pattern::clear()
{
	active_symbol = 0;
	for (auto const& e: symbols){
		symbol_counts[e] = 0;
	}
}

bool Pattern::symbol_seen(std::string symbol)
{
	if (std::find(symbols.begin(), symbols.end(), symbol) != symbols.end()){
		if (active_symbol < symbols.size() && symbols.at(active_symbol) == symbol){
			active_symbol++;
		}
		symbol_counts[symbol] += 1;
		total_symbol_counts[symbol] += 1;
		return true;
	} else {
		return false;
	}
}

std::string Pattern::to_string() const
{
	std::stringstream result;
	for (auto i: symbols){
		result << i << " ";
	}
	return result.str();
}

#ifdef SIGSPAN
double* Pattern::sigspan(double* probabilities, unsigned int pattern_length, unsigned int sequence_length) const{
	double* Qx_ = new double[sequence_length]();
	double* Qx = new double[sequence_length]();

	for (unsigned int i = 0; i < sequence_length; ++i){
		Qx_[i] = probabilities[0] * pow(1-probabilities[0], i);
	}

	for (unsigned int k = 2; k <= pattern_length; ++k){
		for (unsigned int i = 0; i < sequence_length; ++i){
			Qx[i] = 0;
		}
		for (unsigned int i = k; i <= sequence_length; ++i){
			for (unsigned int j = k - 1; j <= i - 1; ++j){
				Qx[i-1] += Qx_[j-1] * probabilities[k - 1] * pow(1-probabilities[k-1], i-j-1);
			}
		}
		double* tmp = Qx_;
		Qx_ = Qx;
		Qx = tmp;
	}
	delete [] Qx;

	for (unsigned int i = 1; i < sequence_length; ++i){
		Qx_[i] += Qx_[i-1];
	}

	return Qx_;
}


double Pattern::get_expected_value_sigspan(std::map<unsigned int, unsigned int> dataset_shape) const
{
	unsigned int max_sequence_length = 0;
	unsigned int dataset_size = 0;
	for (auto const& x: dataset_shape){
		max_sequence_length = std::max(x.first, max_sequence_length);
		dataset_size += x.first * x.second;
	}

	double* probabilities = new double[symbols.size()]();

	for (unsigned int i = 0; i < symbols.size(); ++i){
		probabilities[i] = (double) total_symbol_counts.at(symbols[i]) / dataset_size;
	}

	double* result = sigspan(probabilities, symbols.size(), max_sequence_length);

	if (verbose >= 2){
		// Output item probabilities
		std::cout << "|SigSpan| " << to_string();
		for (unsigned int i = 0; i < symbols.size(); ++i){
			std::cout << " " << probabilities[i];
		}
		std::cout << std::endl;

		// Output probabilities per length
		std::cout << "|SigSpan| " << to_string();
		for (unsigned int i = 0; i < max_sequence_length; ++i){
			std::cout << " " << result[i];
		}
		std::cout << std::endl;

		// Output length counts
		std::cout << "|SigSpan| " << to_string();
		for (unsigned int i = symbols.size(); i <= max_sequence_length; ++i){
			if (dataset_shape.find(i) == dataset_shape.end()){
				std::cout << "0 ";
			} else {
				std::cout << dataset_shape.at(i) << " ";
			}
			
		}
		std::cout << std::endl;
	}

	double expected_support = 0;
	for (unsigned int i = symbols.size(); i <= max_sequence_length; ++i){
		if (dataset_shape.find(i) != dataset_shape.end()){
			expected_support += dataset_shape.at(i) * result[i-1];
		}
	}

	delete [] probabilities;
	delete [] result;

	return expected_support;
}

double Pattern::get_p_sigspan(std::map<unsigned int, unsigned int> dataset_shape) const
{
	double n = 0;
	for (auto const& x: dataset_shape){
		n += x.second;
	}

	return exp((-2.0/n) * pow(get_support() - get_expected_value_sigspan(dataset_shape), 2));
}
#endif

// Memoization variables
std::map<std::pair<int, int>, BigInt> Pattern::g_dict;
std::map<std::vector<unsigned int>, double> Pattern::C_dict;
