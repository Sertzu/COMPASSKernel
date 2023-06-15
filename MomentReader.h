#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <random>
#include <cmath>
#include <map>
#include <algorithm>

class MomentReader {
public:
	MomentReader();
	MomentReader(std::string momentPath);
	std::vector<std::vector<double>> get_magmoms();
	std::vector<std::tuple<int, std::vector<double>>> get_atomCoordinates();
private:
	std::vector<std::vector<double>> readMagmomFile(const std::string& fileName);

	std::vector<std::vector<double>> magmoms;
	std::vector<std::tuple<int, std::vector<double>>> atomCoordinates;
};