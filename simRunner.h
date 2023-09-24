#pragma once
#include "helpers.h"
#include "MomentReader.h"
#include "SettingsClass.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include <memory>

enum class XYZ {X = 0, Y = 1, Z = 2};

struct SimulationVariables
{
	bool isInit = true;

	std::string inputPath = "";
	std::string outputPath = "";
	int atomnum = 0;
	double temperature = -1;
	double HStrength = 0.0;
	XYZ HDir = XYZ::X;

	std::vector<std::vector<std::tuple<int, double>>> links;
	std::vector<std::vector<std::vector<int>>> linksNN;
	std::vector<std::vector<double>> magmoms;
	std::vector<std::tuple<int, std::vector<double>>> atomCoordinates;
};

class SimRunner
{
public:
	SimRunner(std::string settingsPath);

	void LOAD(std::string inputFile = "");
	void SETOUT(std::string outFile = "");
	void SETTEMP(double temp);
	void SETMAGFIELD(double HVal, XYZ dir);

	void APPROACHTEMP(double temperature, int steps);
	void APPROACHMAG(double HVal, int steps);

	void EQUILIB(int steps);
	void MEASUREMENT(int steps);
private:
	SimRunner();

	void print(std::string toPrint);

	void saveMeasurement(std::vector<double> outVals);

	std::unique_ptr<SimEnvironment> simEnv;
	SettingsClass settings;
	SimulationVariables simVars;

	ThreadLogger tlogRunner;
	bool loadFlag = false;
};