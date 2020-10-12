#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <iostream>

class FileReader{
	private:
		std::ifstream m_File;
		std::vector<std::string> m_Line;
		char m_SymbolSeparator;
		char m_LineSeparator;
		bool m_Shuffled;
		bool m_FileFinished;

		bool EnsureLine();

	public:
		FileReader(std::string filename, char symbol_separator, char line_separator, bool shuffled);

		bool Item(std::string& item_data);
		bool Line(std::vector<std::string>& line_data);

		void SetShuffled(bool shuffled);
		void Clear();
};
