#include "multithreadingSim.h"

std::mutex mtx;


void run_sim_thread(unsigned int threadnum, std::string inputPath, std::string outputPath, double temp, int equilibsteps, int measuresteps, const std::vector<std::vector<std::tuple<int, double>>>& links, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinates)
{
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup started for thread: " << threadnum << std::endl;
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temp, threadnum, links, atomCoordinates);
    myEnvironment.setTemperature(temp);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup finished for thread: " << threadnum << std::endl;
    std::cout << getCurrentTime().time_string << " [Thread Handler] Starting simulation in thread:  " << threadnum << std::endl;
    lock.unlock();

    myEnvironment.runSim(equilibsteps, 0);
    myEnvironment.runSim(measuresteps, 1);
    std::vector<double> outVals = myEnvironment.getParameters();
    std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5);

    lock.lock();
    appendToFile(myEnvironment.getOutputPath() + "\\results.dat", out);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Finished simulation in thread:  " << threadnum << std::endl;
    lock.unlock();
    std::cout << getCurrentTime().time_string << " [Thread Handler] Dumping magmoms from thread: " << threadnum << std::endl;
    myEnvironment.writeMagneticMomentsToFile();
    
}

void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec, int equilibsteps, int measurementsteps)
{
    std::vector<std::thread> threads;
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, 0);
    auto links = myEnvironment.getLinks();
    auto atomCoordinates = myEnvironment.getAtomCoordinates();
    // Create and start threads
    for (int i = tempVec.size() - 1; i >= 0; --i) {
        threads.emplace_back(run_sim_thread, i, inputPath, outputPath, tempVec[i], equilibsteps, measurementsteps, links, atomCoordinates);
    }
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << getCurrentTime().time_string << " [Kernel] All threads started!" << std::endl;
    lock.unlock();
    // Wait for all threads to finish
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}
