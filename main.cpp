#include "FileReader.h"
#include "Pattern.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::map<unsigned int, unsigned int> applyFileToPatterns(std::vector<Pattern>* patterns, FileReader* sequenceFile, bool onlyCount, unsigned int verbose){
	// Iterate sequences
	std::map<unsigned int, unsigned int> databaseShape;
	#ifdef SIGSPAN
	unsigned int sequenceLength = 0;
	#endif

	unsigned int sequenceCounter = 0;
	std::string newItem;
	while (sequenceFile->Item(newItem)){
		if (newItem == "\n"){
			if (verbose >= 2) std::cout << std::endl;
			sequenceCounter++;

			#ifdef SIGSPAN
			if (databaseShape.find(sequenceLength) == databaseShape.end()){
				databaseShape[sequenceLength] = 0;
			}
			databaseShape[sequenceLength]++;
			sequenceLength = 0;
			#endif

			for (auto& p: *patterns){
				p.Process(onlyCount);
				p.Clear();
			}
			if (verbose == 1) std::cout << "\r" << sequenceCounter << " sequences processed." << std::flush;
		} else {
			#ifdef SIGSPAN
			sequenceLength++;
			#endif

			if (verbose >= 2) std::cout << newItem << " ";
			for (auto& p: *patterns){
				p.SymbolSeen(newItem, onlyCount);
			}
		}
	}
	if (verbose >= 1) std::cout << std::endl;

	return databaseShape;
}

bool toDouble(char* s, double &result) {
	char* end;
	result = std::strtod(s, &end);
	if (end == s || *end != '\0')
		return false;
	return true;
}

int main(int argc, char** argv)
{
	if (argc < 3) {
		std::cout << argv[0] << " [options] <data> <patterns>" << std::endl;
		std::cout << "output options:" << std::endl;
		std::cout << " -o <filename> output result to file instead of stdio" << std::endl;
		std::cout << " -v Verbose" << std::endl;
		std::cout << " -V Extra verbose" << std::endl;
		std::cout << "Pattern Statistics:" << std::endl;
		std::cout << " -s Output support" << std::endl;
		std::cout << " -c Output number of sequences with non-zero probability" << std::endl;
		std::cout << "PS² options:" << std::endl;
		std::cout << " -e Output expected value" << std::endl;
		std::cout << " -d Output standard deviation" << std::endl;
		std::cout << " -p Output p-value (exact)" << std::endl;
		std::cout << " -P Output -log(p-value) (exact)" << std::endl;
		std::cout << " -n Output p-value (normal approximation)" << std::endl;
		std::cout << " -N Output -log(p-value) (normal approximation)" << std::endl;
		std::cout << " -l Output p-value (Poisson approximation)" << std::endl;
		std::cout << " -L Output -log(p-value) (Poisson approximation)" << std::endl;
		std::cout << "Significance options:" << std::endl;
		std::cout << " -B <alpha> Bonferroni significance threshold" << std::endl;
		std::cout << " -W <alpha> Westfall-Young significance threshold (PS²)" << std::endl;
		#ifdef SIGSPAN
		std::cout << "SigSpan options:" << std::endl;
		std::cout << " -b Output expected value" <<std::endl;
		std::cout << " -i Output p-value" <<std::endl;
		std::cout << " -I Output -log(p-value)" <<std::endl;
		#endif
		return 0;
	}

	unsigned int verbose = 0;
	std::ofstream outputFile;
	double tBonferroni = 0;
	double tWestfallYoung = 0;

	for (unsigned int i = 1; i <= argc-3; ++i){
		if (std::strlen(argv[i]) != 2 or argv[i][0] != '-') continue;
		switch(argv[i][1]){
				case 'v':
					verbose = 1;
					break;
				case 'V':
					verbose = 2;
					break;
				case 'o':
					outputFile.open(std::string(argv[i+1]));
					i += 1;
					break;
				case 'B':
					if (!toDouble(argv[i+1], tBonferroni)){
						std::cout << "-B " << argv[i+1] << " does not define a valid significance threshold, use e.g. -B 0.05" << std::endl;
						return 0;
					}
					if (tBonferroni <= 0 || tBonferroni >= 1){
						std::cout << "-B <alpha> needs to be within range (0,1), e.g. -B 0.05" << std::endl;
						return 0;
					}
					i += 1;
					break;
				case 'W':
					if (!toDouble(argv[i+1], tWestfallYoung)){
						std::cout << "-W " << argv[i+1] << " does not define a valid significance threshold, use e.g. -W 0.05" << std::endl;
						return 0;
					}
					if (tWestfallYoung <= 0 || tWestfallYoung >= 1){
						std::cout << "-W <alpha> needs to be within range (0,1), e.g. -W 0.05" << std::endl;
						return 0;
					}
					i += 1;
					break;
		}
	}

	// Load Patterns
	FileReader patternFile = FileReader(argv[argc - 1], ' ', '\n', false);
	std::vector<Pattern> patterns;
	std::vector<std::string> newSymbol;
	while (patternFile.Line(newSymbol)){
		patterns.push_back(Pattern(newSymbol, verbose));
	}
	if (verbose >= 1) std::cout << patterns.size() << " patterns loaded." << std::endl;

	// Iterate sequences
	FileReader sequenceFile = FileReader(argv[argc - 2], ' ', '\n', false);
	std::map<unsigned int, unsigned int> databaseShape = applyFileToPatterns(&patterns, &sequenceFile, false, verbose);

	// Perform significance tests if requested
	if (tBonferroni != 0){
		std::cout << "Bonferroni significance:" << std::endl;
		std::cout << "  B(" << tBonferroni << ") = " << tBonferroni / patterns.size() << std::endl;
		std::cout << "  -log(B(" << tBonferroni << ")) = " << -log(tBonferroni / patterns.size()) << std::endl;
		std::cout << std::endl;
	}

	if (tWestfallYoung != 0){
		std::cout << "Westfall-Young significance:" << std::endl;

		std::vector<double> ps;
		sequenceFile.SetShuffled(true);
		for (int i = 0; i < 100; ++i){
			std::cout << "\r" << "(" << i+1 << "/100";
			sequenceFile.Clear();
			for (auto& p: patterns){
				p.Reset();
			}
			applyFileToPatterns(&patterns, &sequenceFile, true, verbose);
			double minP = std::numeric_limits<double>::infinity();
			for (auto const& p: patterns){
				minP = std::min(minP, p.PExact());
			}
			ps.push_back(minP);
			std::cout << "\r" << "Sample " << i+1 << "/100" << std::flush;
		}
		double threshold = ps[4];

		std::cout << "\r  W(" << tWestfallYoung << ") = " << threshold << std::endl;
		std::cout << "  -log(W(" << tWestfallYoung << ")) = " << -log(threshold) << std::endl;
		std::cout << std::endl;
	}

	// Output results per pattern
	std::ostream& out_stream = (outputFile.is_open() ? outputFile : std::cout);
	for (auto const& p: patterns){
		std::ostringstream resultString;
		for (unsigned int i = 1; i <= argc-3; ++i){
			if (std::strlen(argv[i]) != 2 or argv[i][0] != '-') continue;
			switch(argv[i][1]){
				case 's':
					resultString << p.Support() << " ";
					break;
				case 'e':
					resultString << p.ExpectedValue() << " ";
					break;
				case 'd':
					resultString << p.StandardDeviation() << " ";
					break;
				case 'c':
					resultString << p.NonZeroSequences() << " ";
					break;
				case 'n':
					resultString << p.PNormal() << " ";
					break;
				case 'N':
					resultString << -log(p.PNormal()) << " ";
					break;
				case 'p':
					resultString << p.PExact() << " ";
					break;
				case 'P':
					resultString << -log(p.PExact()) << " ";
					break;
				case 'l':
					resultString << p.PPoisson() << " ";
					break;
				case 'L':
					resultString << -log(p.PPoisson()) << " ";
					break;
			#ifdef SIGSPAN
				case 'b':
					resultString << p.ExpectedValueSigspan(databaseShape) << " ";
					break;
				case 'i':
					resultString << p.PSigspan(databaseShape) << " ";
					break;
				case 'I':
					resultString << -log(p.PSigspan(databaseShape)) << " ";
					break;
			#endif
				case 'o':
					i +=1;
					break;
			}
		}
		if (resultString.str() != "")
			out_stream << resultString.str() << p.ToString() << std::endl;
	}
	if (outputFile.is_open()){
		outputFile.close();
	}
}
