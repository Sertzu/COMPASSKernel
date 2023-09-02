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
#include <algorithm>
#include <array> 

#include "helpers.h"
#include "timeUtility.h"
#include "BS_thread_pool.hpp"
#include <thread>
#include <mutex>
#include <wincrypt.h>

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
    SimEnvironment(std::string inputPath, std::string outputPath, double temperature, unsigned int threadnum, const std::vector<std::vector<std::tuple<int, double>>>& links, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinates, const std::vector<std::vector<std::vector<int>>>& linksNN, std::vector<double> interactionKZ = { 0.0, 0.0, 0.0 }, std::vector<double> interactionC = { 0.0, 0.0, 0.0 }, std::vector<std::vector<double>> magmoms = {}, std::vector<double> magneticFieldHIn = { 0.0, 0.0, 0.0 });

    std::vector<std::vector<std::tuple<int, double>>> getLinks();
    std::vector<std::vector<std::vector<int>>> getLinksNN();
    std::vector<std::tuple<int, std::vector<double>>> getAtomCoordinates();
    std::vector<std::vector<double>> getMagmoms();

	void runSim(int steps, bool measurement, bool approach = 0, double initialtemp = 300);
	void setTemperature(double temp);

	std::string getOutputPath();
	std::vector<double> getParameters();
    void writeMagneticMomentsToFile();
private:
	SimEnvironment();

	double energy_diff_calculator(int& index, std::vector<double>& oldMom, std::vector<double>& newMom);
    double rate_calculator(int& index,double& beta, std::vector<double>& oldMom, std::vector<double>& newMom);
	double energy_calculator();
	std::vector<double> generateRandomVecSingle();
    void generateRandomVecArray(std::vector<std::vector<double>>& vecIn);

    ThreadLogger tlog;
	unsigned int threadnum;

	std::vector<std::vector<std::tuple<int, double>>> links;
    std::vector<std::vector<std::vector<int>>> linksNN;
	std::vector<std::vector<double>> magmoms;
	std::vector<std::vector<std::vector<double>>> magmomsHistory;
    std::vector<std::tuple<int, std::vector<double>>> atomCoordinates;
	std::vector<double> meanmagmomsHistory;
	std::vector<double> energyHistory;

	std::string inputPath;
	std::string outputPath;
	int atomnum;
	double temperature;
    bool enableCompassAnisotropy;

	const double kB = 8.617330337217213e-02; // meV/K
    const int MAXMAGMOMHISTSIZE = 300;
    std::vector<double> interactionK;
    std::vector<double> interactionC;
    std::vector<double> magneticFieldH;
	std::random_device rd;
	std::mt19937 gen;
	

};