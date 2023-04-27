#include "multithreadingSim.h"

std::mutex mtx;

void print_thread_id(int id)
{
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << "Thread ID: " << id << std::endl;
    lock.unlock();
}

void run_sim_thread(unsigned int threadnum, std::string inputPath, std::string outputPath, double temp, int equilibsteps, int measuresteps)
{
    std::unique_lock<std::mutex> lock(mtx);
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temp, threadnum);
    myEnvironment.setTemperature(temp);
    lock.unlock();

    myEnvironment.runSim(equilibsteps, 0);
    myEnvironment.runSim(measuresteps, 1);
    std::vector<double> outVals = myEnvironment.getParameters();
    std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5);

    lock.lock();
    appendToFile(myEnvironment.getOutputPath() + "\\results.dat", out);
    lock.unlock();
    
}

void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec, int equilibsteps, int measurementsteps)
{
    std::vector<std::thread> threads;

    // Create and start threads
    for (unsigned int i = 0; i < tempVec.size(); ++i) {
        threads.emplace_back(run_sim_thread, i, inputPath, outputPath, tempVec[i], equilibsteps, measurementsteps);
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}
