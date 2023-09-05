#pragma once
#include "helpers.h"
#include "MomentReader.h"
#include "SettingsClass.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include <memory>

struct SimulationVariables
{
	std::string inputPath = "";
	std::string outputPath = "";
	int atomnum = 0;
	double temperature = -1;

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
private:
	SimRunner();

	void print(std::string toPrint);

	std::unique_ptr<SimEnvironment> simEnv;
	SettingsClass settings;
	SimulationVariables simVars;

	ThreadLogger tlogRunner;
	bool loadFlag = false;
};