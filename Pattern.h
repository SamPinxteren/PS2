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
		std::vector<std::string> m_Symbols;

		// Sequence state
		unsigned int m_ActiveSymbol;
		std::map<std::string, unsigned int> m_SymbolCounts;
		std::map<std::string, unsigned int> m_TotalSymbolCounts;

		// Probability statistics
		std::vector<double> m_P;
		double m_ExpectedValue;
		double m_Variance;
		unsigned int m_RealValue;

		// verbosity level
		unsigned int m_Verbose;

		static std::map<std::pair<int, int>, BigInt> m_G;
		BigInt G(int b, int e);

		// Compute the permutation probability exactly
		static std::map<std::vector<unsigned int>, double> m_C;
		double C(std::vector<unsigned int> X);

		double OccursProbability();

		#ifdef SIGSPAN
		double* Sigspan(double* probabilities, unsigned int pattern_length, unsigned int sequence_length) const;
		#endif

	public:
		Pattern(std::vector<std::string> patternSymbols, unsigned int verbosity);
		Pattern(std::vector<std::string> patternSymbols);

		// Get data for currently processed sequences
		double StandardDeviation() const;
		double ExpectedValue() const;
		unsigned int Support() const;
		double PNormal() const;
		double PExact() const;
		double PPoisson() const;
		#ifdef SIGSPAN
		double ExpectedValueSigspan(std::map<unsigned int, unsigned int> dataset_shape) const;
		double PSigspan(std::map<unsigned int, unsigned int> dataset_shape) const;
		#endif
		unsigned int NonZeroSequences() const;

		// Process the last symbols seen
		void Process(bool onlyCount);
		void Process();
		// Handle a new symbol for the current sequence
		bool SymbolSeen(std::string symbol, bool onlyCount);
		bool SymbolSeen(std::string symbol);
		// Clear for new sequence
		void Reset();
		void Clear();

		// Create a string describing this pattern
		std::string ToString() const;
};
#endif
