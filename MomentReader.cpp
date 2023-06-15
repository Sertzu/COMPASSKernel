#include "MomentReader.h"

MomentReader::MomentReader()
{
}

MomentReader::MomentReader(std::string momentPath)
{
    magmoms = readMagmomFile(momentPath);
}

std::vector<std::vector<double>> MomentReader::get_magmoms()
{
	return magmoms;
}

std::vector<std::tuple<int, std::vector<double>>> MomentReader::get_atomCoordinates()
{
	return std::vector<std::tuple<int, std::vector<double>>>();
}

std::vector<std::vector<double>> MomentReader::readMagmomFile(const std::string& fileName)
{
    // Open the file
    std::ifstream file(fileName);

    if (!file) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<std::vector<double>> vectors;
    std::string line;

    while (std::getline(file, line)) {
        // Skip if line is empty or starts with '!' or '#'
        if (line.empty() || line[0] == '!' || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double value;
        std::vector<double> temp;

        for (int i = 0; i < 9; ++i) { // assuming there are always 9 values in a line
            if (!(iss >> value)) {
                throw std::runtime_error("File format error");
            }

            // Only append 4th, 5th and 6th values to temporary vector
            if (i >= 3 && i <= 5) {
                temp.push_back(value);
            }
        }

        vectors.push_back(temp); // Append temporary vector to final vector
    }

    return vectors;
}
