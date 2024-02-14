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
#include <set>

#include "helpers.h"
#include "MomentReader.h"
#include "timeUtility.h"
#include "BS_thread_pool.hpp"
#include <thread>
#include <mutex>
#include <wincrypt.h>
#include <barrier>

typedef std::array<float, 3> Vec3;

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
	SimEnvironment(std::string inputPath, std::string outputPath, float temperature);

    std::vector<std::vector<std::tuple<int, float>>> getLinks();
    std::vector<std::vector<std::vector<int>>> getLinksNN();
    std::vector<std::tuple<int, std::vector<float>>> getAtomCoordinates();
    std::vector<Vec3> getMagmoms();

	void runSim(int steps, bool measurement, bool approachTemp = 0, float initialtemp = 300, bool approachMag = 0, float initialH = 0.0);
	void setTemperature(float temp);

    void setMagneticField(float xH, float yH, float zH);
    void setSingleIonAnisotropy(float xC, float yC, float zC);
    void setCompassAnisotropy(float xComp, float yComp, float zComp);
    void setOutputPath(std::string out);
    void setMagPattern(std::vector<double> mag_pattern);

    void setStatusSteps(int steps);

	std::string getOutputPath();
	std::vector<double> getParameters();
    void writeMagneticMomentsToFile(std::string path);
    void readMagneticMomentsFromFile(std::string path);
private:

    float energy_diff_calculator(const int& index,const Vec3& oldMom,const Vec3& newMom, const std::vector<Vec3>& magmoms, const std::vector<std::vector<std::tuple<int, float>>>& atomlinks);
    float rate_calculator(int& index,float& beta, Vec3& oldMom, Vec3& newMom, const std::vector<Vec3>& magmoms, const std::vector<std::vector<std::tuple<int, float>>>& atomlinks);
    float energy_calculator();
    Vec3 generateRandomVecSingle();

    void generateAcceptanceVec(std::vector<float>& vecIn);
    void generateRandomVecArray(std::vector<Vec3>& vecIn, std::uniform_real_distribution<>& randomizer, std::mt19937& generator);

    ThreadLogger m_tlog;
	unsigned int m_threadnum;

	std::vector<std::vector<std::tuple<int, float>>> m_links;
    std::vector<std::vector<std::vector<int>>> m_linksNN;
	std::vector<Vec3> m_magmoms;
    std::vector<std::tuple<int, std::vector<float>>> m_atomCoordinates;

    std::vector<int> m_atomTypes;
    std::vector<std::string> m_atomNames;
    std::set<int> m_uniqueAtomTypes;
    std::unordered_map<int, std::string> m_atomTypeToName;

    std::vector<std::vector<Vec3>> m_magmomsHistory;
	std::vector<float> m_meanmagmomsHistory;
    std::vector<Vec3> m_meanRawMagmomsHistory;
	std::vector<float> m_energyHistory;

	std::string m_inputPath;
	std::string m_outputPath;
	int m_atomnum;
	float m_temperature;
    bool m_enableCompassAnisotropy;
    bool m_useStaggeredMagnetization = false;
    std::vector<float> m_magnetizationPattern;

	const float m_kB = 8.617330337217213e-02; // meV/K
    const int m_MAXMAGMOMHISTSIZE = 300;
    std::vector<float> m_singleIonAnisotropyTerm;
    std::vector<float> m_compassAnisotropyTerm;
    std::vector<float> m_zeemanTerm;
    int m_magDir;
    int m_singleIonDir;
    int m_compassDir;

	std::random_device rd;
	std::mt19937 gen;
	
    int m_statusSteps = 1000;
    int m_workerCount = 8;
};