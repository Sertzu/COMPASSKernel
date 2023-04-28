#pragma once
#include "helpers.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include <thread>
#include <mutex>

void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec,int equilibsteps, int measuresteps);