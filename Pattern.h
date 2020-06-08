#ifndef PATTERN_H
#define PATTERN_H

#include "BigInt.h"

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

class Pattern
{
	private:
		// Pattern data
		std::vector<std::string> symbols;

		// Sequence state
		unsigned int active_symbol;
		std::map<std::string, unsigned int> symbol_counts;
		std::map<std::string, unsigned int> total_symbol_counts;

		// Probability statistics
		std::vector<double> P;
		double expected_value;
		double variance;
		unsigned int real_value;

		// verbosity level
		unsigned int verbose;

		static std::map<std::pair<int, int>, BigInt> g_dict;
		BigInt g(int b, int e);

		// Compute the permutation probability exactly
		static std::map<std::vector<unsigned int>, double> C_dict;
		double C(std::vector<unsigned int> X);

		double occurs_probability();

		#ifdef SIGSPAN
		void sigspan_recursive(double* result, double* probabilities, unsigned int pattern_length, unsigned int sequence_length, double p) const;
		#endif

	public:
		Pattern(std::vector<std::string> pattern_symbols, unsigned int verbosity);
		Pattern(std::vector<std::string> pattern_symbols);

		// Get data for currently processed sequences
		double get_standard_deviation() const;
		double get_expected_value() const;
		unsigned int get_support() const;
		double get_p_normal() const;
		double get_p_exact() const;
		double get_p_poisson() const;
		#ifdef SIGSPAN
		double get_p_sigspan(std::map<unsigned int, unsigned int> dataset_shape) const;
		#endif
		unsigned int get_non_zero_sequences() const;

		// Process the last symbols seen
		void process();
		// Reset for new sequence
		void clear();
		// Handle a new symbol for the current sequence
		bool symbol_seen(std::string symbol);

		// Create a string describing this pattern
		std::string to_string() const;
};
#endif
