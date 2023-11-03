#include "simRunner.h"
#include <format>

SimRunner::SimRunner(const std::string settingsPath) : settings(settingsPath)
{
	simVars.inputPath = settings.get<std::string>("inputPath");
	simVars.outputPath = settings.get<std::string>("outputPath");
	simVars.temperature = -1;

	tlogRunner = ThreadLogger(0, true);
}

// WILL OVERWRITE EVERYTHING IN THE ENVIRONMENT
void SimRunner::LOAD(std::string inputFile)
{
	if (inputFile != "")
		simVars.inputPath = inputFile;
	std::unique_ptr<SimEnvironment> temp = std::make_unique<SimEnvironment>(simVars.inputPath, simVars.outputPath, simVars.temperature);
	simEnv.swap(temp);
	loadFlag = true;

	print("ENVIRONMENT LOADED FROM FILE: " + simVars.inputPath);
}

void SimRunner::SETOUT(std::string outFile)
{
	if (outFile != "")
		simVars.outputPath = outFile;
	simEnv->setOutputPath(simVars.outputPath);

	print("ENVIRONMENT OUTPUT SET TO FILE: " + simVars.outputPath);
}

void SimRunner::SAVECHECKPOINT(std::string outPath)
{
	simEnv->writeMagneticMomentsToFile(outPath);

	print("ENVIRONMENT SAVED MAGMOMS TO: " + outPath);
}

void SimRunner::SETMAGMOMAUTOSAVE(std::string outPath)
{
	simVars.autoSaveMoments = true;
	simVars.momentsPath = outPath;

	print("ENVIRONMENT WILL AUTOSAVE MAGMOMS TO: " + outPath);
}

void SimRunner::SETTEMP(double temp)
{
	simVars.temperature = temp;
	simEnv->setTemperature(simVars.temperature);

	print("ENVIRONMENT TEMPERATURE SET TO: " + std::to_string(simVars.temperature) + "K");
}

void SimRunner::SETMAGFIELD(double HVal, XYZ dir)
{
	std::string dirString = "";

	switch (dir) {
		case XYZ::X: simEnv->setMagneticField(HVal, 0.0, 0.0); dirString = "X"; break;
		case XYZ::Y: simEnv->setMagneticField(0.0, HVal, 0.0); dirString = "Y"; break;
		case XYZ::Z: simEnv->setMagneticField(0.0, 0.0, HVal); dirString = "Z"; break;
		default: throw std::runtime_error("Invalid direction for the H Field");
	}

	simVars.HDir = dir;
	simVars.HStrength = HVal;

	print("ENVIRONMENT MAGNETIC FIELD SET TO: " + std::to_string(simVars.HStrength) + "T, IN DIRECTION: " + dirString);
}

void SimRunner::SETEASYPLANE(double CVal, XYZ dir)
{
	std::string dirString = "";

	switch (dir) {
	case XYZ::X: simEnv->setSingleIonAnisotropy(CVal, 0.0, 0.0); dirString = "X"; break;
	case XYZ::Y: simEnv->setSingleIonAnisotropy(0.0, CVal, 0.0); dirString = "Y"; break;
	case XYZ::Z: simEnv->setSingleIonAnisotropy(0.0, 0.0, CVal); dirString = "Z"; break;
	default: throw std::runtime_error("Invalid direction for the easy plane anisotropy: ONLY X, Y, Z directions allowed");
	}

	simVars.CDir = dir;
	simVars.CStrength = CVal;

	print("ENVIRONMENT EASY PLANE INTERACTION SET TO: " + std::to_string(simVars.CStrength) + " IN DIRECTION: " + dirString);
}

void SimRunner::SETCOMPASSANISO(double CompVal, XYZ dir)
{
	std::string dirString = "";

	switch (dir) {
	case XYZ::X: simEnv->setCompassAnisotropy(0.0, CompVal, CompVal); dirString = "X"; break;
	case XYZ::Y: simEnv->setCompassAnisotropy(CompVal, 0.0, CompVal); dirString = "Y"; break;
	case XYZ::Z: simEnv->setCompassAnisotropy(CompVal, CompVal, 0.0); dirString = "Z"; break;
	default: throw std::runtime_error("Invalid direction for the compass anisotropy: ONLY X, Y, Z directions allowed");
	}

	simVars.CompDir = dir;
	simVars.CompStrength = CompVal;

	print("ENVIRONMENT COMPASS INTERACTION SET TO: " + std::to_string(simVars.CompStrength) + " IN NORMAL DIRECTION: " + dirString);
}

void SimRunner::APPROACHTEMP(double temperature, int steps)
{
	if (simVars.isInit)
		print("APPROACHING TEMPERATURE: " + std::to_string(temperature) + "K FOR " + std::to_string(steps) + " STEPS");
	else
		throw std::runtime_error("ABORTING SIMULATION, ENVIRONMENT IS NOT FULLY INITIALIZED!");

	simEnv->setTemperature(temperature);
	simEnv->runSim(steps, false, true, simVars.temperature);
	simVars.temperature = temperature;
}

void SimRunner::APPROACHMAG(double HVal, int steps)
{
	if (simVars.isInit)
		print("APPROACHING MAGNETIC FIELD STRENGTH: " + std::to_string(HVal) + "T FOR " + std::to_string(steps) + " STEPS");
	else
		throw std::runtime_error("ABORTING SIMULATION, ENVIRONMENT IS NOT FULLY INITIALIZED!");

	switch (simVars.HDir) {
		case XYZ::X: simEnv->setMagneticField(HVal, 0.0, 0.0); break;
		case XYZ::Y: simEnv->setMagneticField(0.0, HVal, 0.0); break;
		case XYZ::Z: simEnv->setMagneticField(0.0, 0.0, HVal); break;
		default: throw std::runtime_error("Invalid direction for the H Field");
	}
	simEnv->runSim(steps, false, false, simVars.temperature, true, simVars.HStrength + 0.00001);
	simVars.HStrength = HVal;
}

void SimRunner::EQUILIB(int steps)
{
	if (simVars.isInit)
		print("EQUILIB @ TEMPERATURE: " + std::to_string(simVars.temperature) + "K FOR " + std::to_string(steps) + " STEPS");
	else
		throw std::runtime_error("ABORTING SIMULATION, ENVIRONMENT IS NOT FULLY INITIALIZED!");

	simEnv->runSim(steps, false);
}

void SimRunner::MEASUREMENT(int steps)
{
	if (simVars.isInit)
		print("MEASURING @ TEMPERATURE: " + std::to_string(simVars.temperature) + "K FOR " + std::to_string(steps) + " STEPS");
	else
		throw std::runtime_error("ABORTING SIMULATION, ENVIRONMENT IS NOT FULLY INITIALIZED!");

	simEnv->runSim(steps, true);

	saveMeasurement(simEnv->getParameters());
	if (simVars.autoSaveMoments)
		simEnv->writeMagneticMomentsToFile(simVars.momentsPath);
}

void SimRunner::SWEEPTEMP(double targetTemp, double tempSteps, int approachSteps, int equilibSteps, int measurementSteps)
{
	if (targetTemp > simVars.temperature)
		throw std::logic_error("TARGET TEMP IS ABOVE CURRENT TEMP");

	double currentTemp = simVars.temperature;

	print("STARTING TEMPERATURE SWEEP FROM: " + std::to_string(currentTemp) + "K TO: " + std::to_string(targetTemp) + "K\n");

	while (currentTemp > targetTemp)
	{
		APPROACHTEMP(currentTemp, approachSteps);
		EQUILIB(equilibSteps);
		MEASUREMENT(measurementSteps);

		currentTemp -= tempSteps;
	}
	simVars.temperature = currentTemp;
	print("FINISHED TEMPERATURE SWEEP FROM: " + std::to_string(currentTemp) + "K TO: " + std::to_string(targetTemp) + "K\n");
}

void SimRunner::SWEEPMAG(double targetH, double HSteps, int approachSteps, int equilibSteps, int measurementSteps)
{
	double currentFieldStrength = simVars.HStrength;

	print("STARTING MAGNETIC FIELD SWEEP FROM: " + std::to_string(currentFieldStrength) + "T TO: " + std::to_string(targetH) + "T\n");
	bool targetHigher = (targetH - currentFieldStrength > 0.);

	if(targetHigher)
		while (currentFieldStrength < targetH)
		{
			if (HSteps <= 0.)
				throw std::runtime_error("ENDLESS LOOP IN SWEEPMAG");

			APPROACHMAG(currentFieldStrength, approachSteps);
			EQUILIB(equilibSteps);
			MEASUREMENT(measurementSteps);

			currentFieldStrength += HSteps;
		}
	else
		while (currentFieldStrength > targetH)
		{
			if (HSteps >= 0.)
				throw std::runtime_error("ENDLESS LOOP IN SWEEPMAG");

			APPROACHMAG(currentFieldStrength, approachSteps);
			EQUILIB(equilibSteps);
			MEASUREMENT(measurementSteps);

			currentFieldStrength += HSteps;
		}
	simVars.HStrength = currentFieldStrength;
	print("FINISHED MAGNETIC FIELD SWEEP FROM: " + std::to_string(currentFieldStrength) + "T TO: " + std::to_string(targetH) + "T\n");
}


void SimRunner::print(const std::string toPrint)
{
	auto time = getCurrentTime().time_string;
	auto msg = toPrint;
	tlogRunner << (time + " [SimRunner] " + msg + "\n");
}

void SimRunner::saveMeasurement(std::vector<double> outVals)
{
	print("SAVING MEASUREMENT TO FILE: " + simVars.outputPath);
	std::string out(formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5) + "  " + formatDouble(outVals[6], 9, 5) + "  " + formatDouble(outVals[7], 9, 5));

	appendToFile(simVars.outputPath, out);
}
