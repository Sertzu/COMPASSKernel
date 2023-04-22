#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <tuple>
#include <random>
#include <cmath>

#include "helpers.h"

class SimEnvironment
{
public:
	SimEnvironment(std::string inputPath, std::string outputPath, double temperature);
	void runSim(int steps, bool measurement);
	std::vector<double> getParameters();
private:
	SimEnvironment();

	double energy_calculator(int& index, std::vector<double>& oldMom, std::vector<double>& newMom);
	std::vector<double> generate_random_vec();

	std::vector<std::vector<std::tuple<int, double>>> links;
	std::vector<std::vector<double>> magmoms;
	std::vector<std::vector<std::vector<double>>> magmomsHistory;
	std::vector<double> meanmagmomsHistory;

	std::string inputPath;
	std::string outputPath;
	int atomnum;
	double temperature;

	double kB = 8.617330337217213e-05;
	std::random_device rd;
	std::mt19937 gen;
	

};