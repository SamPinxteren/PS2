#include "Pattern.h"

BigInt Pattern::G(int b, int e)
{
	if (b < e){
		int temp = e;
		e = b;
		b = temp;
	}

	BigInt f;
	auto search = m_G.find(std::make_pair(b, e));
	if (search == m_G.end()){
		for(int i = b + 1; i <= b + e; ++i){
			f *= i;
		}

		for(int i = 2; i <= e; ++i){
			f /= i;
		}
		m_G[std::make_pair(b, e)] = f;
	} else {
		f = search->second;
	}

	return f;
}

double Pattern::C(std::vector<unsigned int> X)
{
	if (X.size() == 1) return 1.0;
	if (m_C.find(X) != m_C.end()) return m_C[X];

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
		BigInt norm_term = G(l, X[n]);

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
					BigInt f2 = G(j, t);
					f2 *= G(l - p, X[n] - t - 1);
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

	m_C[X] = result;
	return result;
}

double Pattern::OccursProbability()
{
	std::vector<unsigned int> X;
	for (auto const& e: m_SymbolCounts){
		if (e.second == 0) return 0;
		X.push_back(e.second);
	}
	std::sort(X.begin(), X.end());
	return C(X);
}

Pattern::Pattern(std::vector<std::string> patternSymbols, unsigned int verbosity):
	Pattern{patternSymbols}
{
	m_Verbose = verbosity;
}

Pattern::Pattern(std::vector<std::string> patternSymbols):
	m_ExpectedValue(0),
	m_Variance(0),
	m_Verbose(0),
	m_ActiveSymbol(0),
	m_Symbols{patternSymbols}
{
	// Ensure input does not contain duplicates
	std::vector<std::string> unique_list = patternSymbols;
	std::sort(unique_list.begin(), unique_list.end());
	auto u = std::unique(unique_list.begin(), unique_list.end());
	if (std::distance(unique_list.begin(), u) != patternSymbols.size()){
		std::stringstream error_msg;
		error_msg << "Duplicates in pattern: ";
		for (auto i: patternSymbols){
			error_msg << i << " ";
		}
		throw std::domain_error(error_msg.str());
	}

	// Initialize counters
	for (auto const& e: m_Symbols){
		m_SymbolCounts[e] = 0;
		m_TotalSymbolCounts[e] = 0;
	}
}

double Pattern::StandardDeviation() const
{
	return sqrt(m_Variance);
}

double Pattern::ExpectedValue() const
{
	return m_ExpectedValue;
}

unsigned int Pattern::Support() const
{
	return m_RealValue;
}

double Pattern::PNormal() const
{
	double z = (Support() - 0.5 - ExpectedValue()) / StandardDeviation();
	return erfc(z / sqrt(2.0)) / 2.0;
}

double Pattern::PExact() const
{
	double* Q = new double[m_P.size()+1]();
	Q[0] = 1;
	for (const double p: m_P){
		for (int i = m_P.size(); i >= 1; --i){
			Q[i] = Q[i]*(1-p) + Q[i-1]*p;
		}
		Q[0] = Q[0]*(1-p);
	}
	double p = 0;
	for (int i = Support(); i <= m_P.size(); ++i){
		p += Q[i];
	}

	delete [] Q;
	return p;
}

double Pattern::PPoisson() const
{
	double lambda = ExpectedValue();
	unsigned int val = Support();
	unsigned int max_val = NonZeroSequences();
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

unsigned int Pattern::NonZeroSequences() const
{
	return m_P.size();
}

void Pattern::Process(bool onlyCount)
{
	if (onlyCount == false){
		Process();
	} else {
		int occuring = (m_ActiveSymbol == m_Symbols.size());
		if (m_Verbose >= 2){
			std::cout
				<< ToString()
				<< " occuring=" << occuring
				<< std::endl;
		}
		m_RealValue += occuring;
	}
}

void Pattern::Process()
{
	double occurs_prob = OccursProbability();
	double var = occurs_prob * (1.0 - occurs_prob);
	int occuring = (m_ActiveSymbol == m_Symbols.size());

	if (m_Verbose >= 2){
		std::cout
			<< ToString()
			<< " occuring=" << occuring
			<< " var=" << var
			<< " p=" << occurs_prob
			<< std::endl;
	}

	m_ExpectedValue += occurs_prob;
	m_RealValue += occuring;
	m_Variance += var;
	if (occurs_prob > 0){
		m_P.push_back(occurs_prob);
	}
}

bool Pattern::SymbolSeen(std::string symbol)
{
	return SymbolSeen(symbol, false);
}

bool Pattern::SymbolSeen(std::string symbol, bool onlyCount)
{
	if (std::find(m_Symbols.begin(), m_Symbols.end(), symbol) != m_Symbols.end()){
		if (m_ActiveSymbol < m_Symbols.size() && m_Symbols.at(m_ActiveSymbol) == symbol){
			m_ActiveSymbol++;
		}
		m_SymbolCounts[symbol] += 1;
		if (!onlyCount)
			m_TotalSymbolCounts[symbol] += 1;
		return true;
	} else {
		return false;
	}
}

void Pattern::Reset()
{
	Clear();
	m_RealValue = 0;
}

void Pattern::Clear()
{
	m_ActiveSymbol = 0;
	for (auto const& e: m_Symbols){
		m_SymbolCounts[e] = 0;
	}
}

std::string Pattern::ToString() const
{
	std::stringstream result;
	for (auto i: m_Symbols){
		result << i << " ";
	}
	return result.str();
}

#ifdef SIGSPAN
double* Pattern::Sigspan(double* probabilities, unsigned int pattern_length, unsigned int sequence_length) const{
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


double Pattern::ExpectedValueSigspan(std::map<unsigned int, unsigned int> dataset_shape) const
{
	unsigned int max_sequence_length = 0;
	unsigned int dataset_size = 0;
	for (auto const& x: dataset_shape){
		max_sequence_length = std::max(x.first, max_sequence_length);
		dataset_size += x.first * x.second;
	}

	double* probabilities = new double[m_Symbols.size()]();

	for (unsigned int i = 0; i < m_Symbols.size(); ++i){
		probabilities[i] = (double) m_TotalSymbolCounts.at(m_Symbols[i]) / dataset_size;
	}

	double* result = Sigspan(probabilities, m_Symbols.size(), max_sequence_length);

	if (m_Verbose >= 2){
		// Output item probabilities
		std::cout << "|SigSpan| " << ToString();
		for (unsigned int i = 0; i < m_Symbols.size(); ++i){
			std::cout << " " << probabilities[i];
		}
		std::cout << std::endl;

		// Output probabilities per length
		std::cout << "|SigSpan| " << ToString();
		for (unsigned int i = 0; i < max_sequence_length; ++i){
			std::cout << " " << result[i];
		}
		std::cout << std::endl;

		// Output length counts
		std::cout << "|SigSpan| " << ToString();
		for (unsigned int i = m_Symbols.size(); i <= max_sequence_length; ++i){
			if (dataset_shape.find(i) == dataset_shape.end()){
				std::cout << "0 ";
			} else {
				std::cout << dataset_shape.at(i) << " ";
			}
			
		}
		std::cout << std::endl;
	}

	double expected_support = 0;
	for (unsigned int i = m_Symbols.size(); i <= max_sequence_length; ++i){
		if (dataset_shape.find(i) != dataset_shape.end()){
			expected_support += dataset_shape.at(i) * result[i-1];
		}
	}

	delete [] probabilities;
	delete [] result;

	return expected_support;
}

double Pattern::PSigspan(std::map<unsigned int, unsigned int> dataset_shape) const
{
	double n = 0;
	for (auto const& x: dataset_shape){
		n += x.second;
	}

	return exp((-2.0/n) * pow(Support() - ExpectedValueSigspan(dataset_shape), 2));
}
#endif

// Memoization variables
std::map<std::pair<int, int>, BigInt> Pattern::m_G;
std::map<std::vector<unsigned int>, double> Pattern::m_C;
