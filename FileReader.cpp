#include "FileReader.h"

FileReader::FileReader(std::string filename, char symbol_separator, char line_separator){
	file.open(filename);
	sym_sep = symbol_separator;
	line_sep = line_separator;
	std::string line;
	std::getline(file, line, line_sep);
	line_stream.str(line);
}

bool FileReader::getitem(std::string& item_data){
	item_data.clear();

	if (line_stream.str().size() == 0){
		std::string line;
		if (!std::getline(file, line, line_sep)){
			return false;
		}
		line_stream.str(line);
		line_stream.clear();
	}

	if (!std::getline(line_stream, item_data, sym_sep)){
		item_data = "\n";
		line_stream.str("");
		line_stream.clear();
		return true;
	}
	return true;
}

bool FileReader::getline(std::vector<std::string>& line_data){
	line_data.clear();
	std::string symbol = "";
	while (getitem(symbol) && symbol != "\n"){
		line_data.push_back(symbol);
	}
	return symbol == "\n";
}
