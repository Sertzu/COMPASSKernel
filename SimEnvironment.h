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
#include <barrier>

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
    SimEnvironment();
	SimEnvironment(std::string inputPath, std::string outputPath, double temperature);
	SimEnvironment(std::string inputPath, std::string outputPath, double temperature, unsigned int threadnum);
    SimEnvironment(std::string inputPath, std::string outputPath, double temperature, unsigned int threadnum, const std::vector<std::vector<std::tuple<int, double>>>& links, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinates, const std::vector<std::vector<std::vector<int>>>& linksNN, std::vector<double> interactionComp = { 0.0, 0.0, 0.0 }, std::vector<double> interactionC = { 0.0, 0.0, 0.0 }, std::vector<std::vector<double>> magmoms = {}, std::vector<double> magneticFieldHIn = { 0.0, 0.0, 0.0 });

    std::vector<std::vector<std::tuple<int, double>>> getLinks();
    std::vector<std::vector<std::vector<int>>> getLinksNN();
    std::vector<std::tuple<int, std::vector<double>>> getAtomCoordinates();
    std::vector<std::vector<double>> getMagmoms();

	void runSim(int steps, bool measurement, bool approachTemp = 0, double initialtemp = 300, bool approachMag = 0, double initialH = 0.0);
	void setTemperature(double temp);

    void setMagneticField(double xH, double yH, double zH);
    void setSingleIonAnisotropy(double xC, double yC, double zC);
    void setCompassAnisotropy(double xComp, double yComp, double zComp);
    void setOutputPath(std::string out);
    void setMagPattern(std::vector<double> mag_pattern);

    void setStatusSteps(int steps);

	std::string getOutputPath();
	std::vector<double> getParameters();
    void writeMagneticMomentsToFile(std::string path);
private:

	double energy_diff_calculator(const int& index,const std::vector<double>& oldMom,const std::vector<double>& newMom, const std::vector<std::vector<double>>& magmoms, const std::vector<std::vector<std::tuple<int, double>>>& atomlinks);
    double rate_calculator(int& index,double& beta, std::vector<double>& oldMom, std::vector<double>& newMom, const std::vector<std::vector<double>>& magmoms, const std::vector<std::vector<std::tuple<int, double>>>& atomlinks);
	double energy_calculator();
	std::vector<double> generateRandomVecSingle();

    void generateAcceptanceVec(std::vector<double>& vecIn);
    void generateRandomVecArray(std::vector<std::vector<double>>& vecIn, std::uniform_real_distribution<>& randomizer, std::mt19937& generator);

    ThreadLogger m_tlog;
	unsigned int m_threadnum;
    int m_workerCount = 6;

	std::vector<std::vector<std::tuple<int, double>>> m_links;
    std::vector<std::vector<std::vector<int>>> m_linksNN;
	std::vector<std::vector<double>> m_magmoms;
    std::vector<std::tuple<int, std::vector<double>>> m_atomCoordinates;
    std::vector<int> m_atomTypes;
    std::vector<std::string> m_atomNames;

    std::vector<std::vector<std::vector<double>>> m_magmomsHistory;
	std::vector<double> m_meanmagmomsHistory;
	std::vector<double> m_energyHistory;

	std::string m_inputPath;
	std::string m_outputPath;
	int m_atomnum;
	double m_temperature;
    bool m_enableCompassAnisotropy;
    bool m_useStaggeredMagnetization = false;
    std::vector<double> m_magnetizationPattern;

	const double m_kB = 8.617330337217213e-02; // meV/K
    const int m_MAXMAGMOMHISTSIZE = 300;
    std::vector<double> m_singleIonAnisotropyTerm;
    std::vector<double> m_compassAnisotropyTerm;
    std::vector<double> m_zeemanTerm;
    int m_magDir;
    int m_singleIonDir;
    int m_compassDir;

	std::random_device rd;
	std::mt19937 gen;
	
    int m_statusSteps = 1000;
};