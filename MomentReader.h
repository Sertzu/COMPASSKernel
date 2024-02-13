#pragma once
#include <string>
#include <vector>
#include <array>
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

typedef std::array<float, 3> Vec3;

class MomentReader {
public:
	MomentReader(std::string momentPath);
	std::vector<Vec3> get_magmoms();
private:
	MomentReader();
	void readMagmomFile(const std::string& fileName);

	std::vector<Vec3> m_magmoms;
};