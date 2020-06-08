#include <string>
#include <fstream>
#include <sstream>
#include <vector>

class FileReader{
	private:
		std::ifstream file;
		std::stringstream line_stream;
		char sym_sep;
		char line_sep;

	public:
		FileReader(std::string filename, char symbol_separator, char line_separator);

		bool getitem(std::string& item_data);
		bool getline(std::vector<std::string>& line_data);
};
