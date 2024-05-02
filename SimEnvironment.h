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
#include <memory>

#include "helpers.h"
#include "MomentReader.h"
#include "timeUtility.h"
#include "BS_thread_pool.hpp"
#include "3DGrid.h"
#include <thread>
#include <mutex>
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

struct IndivdualParameters
{
    std::string atomName;
    int atomType = 0;

    std::vector<double> parameters;
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
    Vec3 getSingeIonAnisotropy();
    void setCompassAnisotropy(float xComp, float yComp, float zComp);
    void setOutputPath(std::string out);
    void setMagPattern(std::vector<double> mag_pattern);
    void setWorkerCount(int workerCount);


    void setStatusSteps(int steps);

	std::string getOutputPath();
	std::vector<double> getParameters();
    std::vector<IndivdualParameters> getIndivdualParameters();

    void writeMagneticMomentsToFile(const std::string& path);
    void readMagneticMomentsFromFile(const std::string& path);
    void writeCorrelationFunction(const std::string& path);
private:

    float energy_diff_calculator(const int& index,const Vec3& oldMom,const Vec3& newMom, const std::vector<Vec3>& magmoms, const std::vector<std::vector<std::tuple<int, float>>>& atomlinks);
    float rate_calculator(int& index,float& beta, Vec3& oldMom, Vec3& newMom, const std::vector<Vec3>& magmoms, const std::vector<std::vector<std::tuple<int, float>>>& atomlinks);
    float energy_calculator(std::array<float, 20>& individualEnergy);
    Vec3 generateRandomVecSingle();

    void generateAcceptanceVec(std::vector<float>& vecIn);
    void generateRandomVecArray(std::vector<Vec3>& vecIn, std::uniform_real_distribution<>& randomizer, std::mt19937& generator);

    std::array<std::vector<float>, 3> calculateCorrelationFunction(int atomType);

    ThreadLogger m_tlog;
	unsigned int m_threadnum;

	std::vector<std::vector<std::tuple<int, float>>> m_links;
    std::vector<std::vector<std::vector<int>>> m_linksNN;
	std::vector<Vec3> m_magmoms;
    std::vector<std::array<int, 3>> m_atomGridPos;
    std::vector<std::tuple<int, std::vector<float>>> m_atomCoordinates;

    std::vector<int> m_atomTypes;
    std::vector<std::string> m_atomNames;
    std::set<int> m_uniqueAtomTypes;
    std::unordered_map<int, std::string> m_atomTypeToName;

    std::vector<std::vector<std::unique_ptr<Array3D>>> m_gridHistory;
    std::vector<std::vector<Vec3>> m_magmomsHistory;
	std::vector<float> m_meanmagmomsHistory;
    std::vector<Vec3> m_meanRawMagmomsHistory;
	std::vector<float> m_energyHistory;

    std::vector<std::array<float, 20>> m_individualMeanmagmomsHistory;
    std::vector<std::array<float, 20>> m_individualEnergyHistory;

	std::string m_inputPath;
	std::string m_outputPath;
	int m_atomnum;
    std::array<int, 20> m_individualAtomnum;
    std::array<int, 3> m_gridSize;

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
    int m_gridSteps = 100;
    int m_workerCount = 1;
};