#include "FileReader.h" 

FileReader::FileReader(std::string filename, char symbolSeparator, char lineSeparator, bool shuffled) :
	m_SymbolSeparator(symbolSeparator),
	m_LineSeparator(lineSeparator),
	m_Shuffled(shuffled) {
	m_File.open(filename);
	Clear();
}

bool FileReader::EnsureLine(){
	std::string line;
	while (m_Line.size() == 0){
		if (!std::getline(m_File, line, m_LineSeparator))
			return false;
		std::string newElement;
		std::stringstream lineStream = std::stringstream(line);
		while (std::getline(lineStream, newElement, m_SymbolSeparator))
			m_Line.push_back(newElement);

		if (m_Shuffled)
			std::random_shuffle(m_Line.begin(), m_Line.end());
	}
	return true;
}

bool FileReader::Item(std::string& itemData){
	if (m_FileFinished){
		return false;
	} else if (m_Line.size() == 0){
		m_FileFinished = !EnsureLine();
		itemData = "\n";
	} else {
		itemData = m_Line[0];
		m_Line.erase(m_Line.begin());
	}
	return true;
}

bool FileReader::Line(std::vector<std::string>& lineData){
	if (!EnsureLine()) return false;

	lineData = m_Line;
	m_Line.clear();
}

void FileReader::SetShuffled(bool shuffled){
	m_Shuffled = shuffled;
}

void FileReader::Clear(){
	m_File.clear();
	m_File.seekg(0);
	m_FileFinished = false;
	EnsureLine();
}
