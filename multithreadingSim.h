#pragma once
#include "helpers.h"
#include "SimEnvironment.h"
#include <thread>
#include <mutex>

void print_thread_id(int id);
void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec,int equilibsteps, int measuresteps);