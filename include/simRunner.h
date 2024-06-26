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

	std::string momentsPath = "";
	std::string individualAtomPath = "";
	std::string correlationFunctionPath = "";

	bool autoSaveMoments = false;
	int workerCount = 1;

	double HStrength = 0.0;
	XYZ HDir = XYZ::X;

	double CStrength = 0.0;
	XYZ CDir = XYZ::X;

	double CompStrength = 0.0;
	XYZ CompDir = XYZ::X;

	std::vector<double> mag_pattern;
	std::vector<std::vector<std::tuple<int, double>>> links;
	std::vector<std::vector<std::vector<int>>> linksNN;
	std::vector<std::vector<double>> magmoms;
	std::vector<std::tuple<int, std::vector<double>>> atomCoordinates;
};

class SimRunner
{
public:
	SimRunner();

	void LOAD(std::string inputFile = "");
	void SETOUT(std::string outFile = "");
	void LOADCHECKPOINT(std::string outPath = "");
	void SAVECHECKPOINT(std::string outPath = "");
	void SETMAGMOMAUTOSAVE(std::string outPath = "");
	void SETMAGPATTERN(std::vector<double> mag_pattern);
	void SETWORKERCOUNT(int workerCount);

	void SETTEMP(double temp);
	void SETMAGFIELD(double HVal, XYZ dir);

	void SETEASYPLANE(double CVal, XYZ dir);
	void SETCOMPASSANISO(double CompVal, XYZ dir);

	void APPROACHTEMP(double temperature, int steps);

	void APPROACHMAG(double HVal, int steps);

	void EQUILIB(int steps);
	void MEASUREMENT(int steps);

	void SWEEPTEMP(double targetTemp, double tempSteps, int approachSteps, int equilibSteps, int measurementSteps);

	void SWEEPMAG(double targetH, double HSteps, int approachSteps, int equilibSteps, int measurementSteps);
private:
	void print(const std::string toPrint);

	void saveMeasurement(std::vector<double> outVals);
	void saveIndividualParameters(std::vector<IndivdualParameters> params);
	std::unique_ptr<SimEnvironment> simEnv;
	SimulationVariables simVars;

	ThreadLogger tlogRunner;
	bool loadFlag = false;
};