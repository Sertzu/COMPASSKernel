#include "simRunner.h"
#include <format>

SimRunner::SimRunner(std::string settingsPath) : settings(settingsPath)
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

void SimRunner::SETTEMP(double temp)
{
	simVars.temperature = temp;
	simEnv->setTemperature(simVars.temperature);

	print("ENVIRONMENT TEMPERATURE SET TO: " + std::to_string(simVars.temperature));
}

void SimRunner::print(std::string toPrint)
{
	auto time = getCurrentTime().time_string;
	auto msg = toPrint;
	tlogRunner << (time + " [SimRunner] " + msg + "\n");
}
