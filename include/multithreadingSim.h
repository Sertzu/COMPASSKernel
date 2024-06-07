#pragma once
#include "helpers.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include "SettingsClass.h"
#include "MomentReader.h"
#include <thread>
#include <mutex>


// void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec, int equilibsteps, int measuresteps, std::vector<double> KVec = {0.0}, std::vector<double> CVec = { 0.0 });
// void create_and_run_threads(SettingsClass simSettings);