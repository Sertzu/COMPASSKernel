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

#include <thread>
#include <mutex>

class ThreadLogger {
public:
    ThreadLogger(int threadnum = 0, bool writeToConsole = false) : writeToConsole(writeToConsole), threadnum(threadnum) {
        if (!writeToConsole) {
            std::ostringstream filename;
            filename << "SimLog_" << threadnum << ".txt";
            file.open(filename.str());
        }
    }

    template <typename T>
    ThreadLogger& operator<<(const T& value) {
        if (writeToConsole) {
            std::cout << value;
        }
        else {
            file << value;
        }
        return *this;
    }

    // Handling manipulators like std::endl
    ThreadLogger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (writeToConsole) {
            std::cout << manip;
        }
        else {
            file << manip;
        }
        return *this;
    }

private:
    int threadnum;
    bool writeToConsole;
    std::ofstream file;
};



class SimEnvironment
{
public:
	SimEnvironment(std::string inputPath, std::string outputPath, double temperature);
	SimEnvironment(std::string inputPath, std::string outputPath, double temperature, unsigned int threadnum);
	void runSim(int steps, bool measurement);
	void setTemperature(double temp);
	std::string getOutputPath();
	std::vector<double> getParameters();
private:
	SimEnvironment();

	double energy_diff_calculator(int& index, std::vector<double>& oldMom, std::vector<double>& newMom);
	double energy_calculator();
	std::vector<double> generate_random_vec();

    ThreadLogger tlog;
	unsigned int threadnum;

	std::vector<std::vector<std::tuple<int, double>>> links;
	std::vector<std::vector<double>> magmoms;
	std::vector<std::vector<std::vector<double>>> magmomsHistory;
	std::vector<double> meanmagmomsHistory;
	std::vector<double> energyHistory;

	std::string inputPath;
	std::string outputPath;
	int atomnum;
	double temperature;

	double kB = 8.617330337217213e-05;
	std::random_device rd;
	std::mt19937 gen;
	

};