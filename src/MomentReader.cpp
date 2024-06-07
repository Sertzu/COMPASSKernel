#include "MomentReader.h"

MomentReader::MomentReader(std::string momentPath)
{
    readMagmomFile(momentPath);
}

std::vector<Vec3> MomentReader::get_magmoms()
{
	return m_magmoms;
}

void MomentReader::readMagmomFile(const std::string& fileName)
{
    // Open the file
    std::ifstream file(fileName);

    if (!file) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<Vec3> temp_magmoms;
    std::string line;

    while (std::getline(file, line)) {
        // Skip if line is empty or starts with '!' or '#'
        if (line.empty() || line[0] == '!' || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double value;
        Vec3 temp = { 0.f,0.f,1.f };
        int pos = 0;

        for (int i = 0; i < 6; ++i) { // assuming there are always 6 values in a line
            if (!(iss >> value)) {
                throw std::runtime_error("File format error");
            }

            // Only append 4th, 5th and 6th values to temporary vector
            if (i >= 3 && i <= 5) {
                temp[pos] = value;
                pos++;
            }
        }

        temp_magmoms.push_back(temp); // Append temporary vector to final vector
    }
    m_magmoms = temp_magmoms;
}
