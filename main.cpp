#include "FileReader.h"
#include "Pattern.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 3) {
		std::cout << argv[0] << " [options] <data> <patterns>" << std::endl;
		std::cout << "output options:" << std::endl;
		std::cout << " -o <filename> output result to file instead of stdio" << std::endl;
		std::cout << " -v verbose" << std::endl;
		std::cout << " -V extra verbose" << std::endl;
		std::cout << "score options:" << std::endl;
		std::cout << " -s output support" << std::endl;
		std::cout << " -e output expected value" << std::endl;
		std::cout << " -d output standard deviation" << std::endl;
		std::cout << " -c output number of sequences with non-zero probability" << std::endl;
		std::cout << " -p output p-value (exact)" << std::endl;
		std::cout << " -P output -log(p-value) (exact)" << std::endl;
		std::cout << " -n output p-value (normal approximation)" << std::endl;
		std::cout << " -N output -log(p-value) (normal approximation)" << std::endl;
		std::cout << " -l output p-value (Poisson approximation)" << std::endl;
		std::cout << " -L output -log(p-value) (Poisson approximation)" << std::endl;
		#ifdef SIGSPAN
		std::cout << "sigspan options:" << std::endl;
		std::cout << " -i output p-value" <<std::endl;
		std::cout << " -I output -log(p-value)" <<std::endl;
		#endif
		return 0;
	}

	unsigned int verbose = 0;
	std::ofstream output_file;

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
					output_file.open(std::string(argv[i+1]));
					i += 1;
					break;
		}
	}

	// Load Patterns
	FileReader pattern_file = FileReader(argv[argc-1], ' ', '\n');

	std::vector<Pattern> patterns;
	std::vector<std::string> new_symbol;
	while (pattern_file.getline(new_symbol)){
		patterns.push_back(Pattern(new_symbol, verbose));
	}

	if (verbose >= 1) std::cout << patterns.size() << " patterns loaded." << std::endl;

	// Iterate sequences
	#ifdef SIGSPAN
	std::map<unsigned int, unsigned int> database_shape;
	unsigned int sequence_length = 0;
	#endif

	unsigned int sequence_counter = 0;
	FileReader sequence_file = FileReader(argv[argc-2], ' ', '\n');
	std::string new_item;
	while (sequence_file.getitem(new_item)){
		if (new_item == "\n"){
			if (verbose >= 2) std::cout << std::endl;
			sequence_counter++;

			#ifdef SIGSPAN
			if (database_shape.find(sequence_length) == database_shape.end()){
				database_shape[sequence_length] = 0;
			}
			database_shape[sequence_length]++;
			sequence_length = 0;
			#endif

			for (auto& p: patterns){
				p.process();
				p.clear();
			}
			if (verbose == 1) std::cout << "\r" << sequence_counter << " sequences processed." << std::flush;
		} else {
			#ifdef SIGSPAN
			sequence_length++;
			#endif

			if (verbose >= 2) std::cout << new_item << " ";
			for (auto& p: patterns){
				p.symbol_seen(new_item);
			}
		}
	}
	if (verbose >= 1) std::cout << std::endl;

	std::ostream& out_stream = (output_file.is_open() ? output_file : std::cout);
	for (auto const& p: patterns){
		for (unsigned int i = 1; i <= argc-3; ++i){
			if (std::strlen(argv[i]) != 2 or argv[i][0] != '-') continue;
			switch(argv[i][1]){
				case 's':
					out_stream << p.get_support() << " ";
					break;
				case 'e':
					out_stream << p.get_expected_value() << " ";
					break;
				case 'd':
					out_stream << p.get_standard_deviation() << " ";
					break;
				case 'c':
					out_stream << p.get_non_zero_sequences() << " ";
					break;
				case 'n':
					out_stream << p.get_p_normal() << " ";
					break;
				case 'N':
					out_stream << -log(p.get_p_normal()) << " ";
					break;
				case 'p':
					out_stream << p.get_p_exact() << " ";
					break;
				case 'P':
					out_stream << -log(p.get_p_exact()) << " ";
					break;
				case 'l':
					out_stream << p.get_p_poisson() << " ";
					break;
				case 'L':
					out_stream << -log(p.get_p_poisson()) << " ";
					break;
			#ifdef SIGSPAN
				case 'i':
					out_stream << p.get_p_sigspan(database_shape) << " ";
					break;
				case 'I':
					out_stream << -log(p.get_p_sigspan(database_shape)) << " ";
					break;
			#endif
				case 'o':
					i +=1;
					break;
			}
		}
		out_stream << p.to_string() << std::endl;
	}
	if (output_file.is_open()){
		output_file.close();
	}
}
