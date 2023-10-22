#include "multithreadingSim.h"

std::mutex mtx;


void run_sim_thread(unsigned int threadnum, std::string inputPath, std::string outputPath, double temp, int equilibsteps, int measuresteps, const std::vector<std::vector<std::tuple<int, double>>>& links, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinates, const std::vector<std::vector<std::vector<int>>>& linksNN, double KVec, double CVec)
{
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup started for thread: " << threadnum << std::endl;
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temp, threadnum, links, atomCoordinates, linksNN, { 0.0, 0.0, KVec }, {CVec, CVec, 0.0});
    myEnvironment.setTemperature(temp);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup finished for thread: " << threadnum << std::endl;
    std::cout << getCurrentTime().time_string << " [Thread Handler] Starting simulation in thread:  " << threadnum << std::endl;
    lock.unlock();

    //TESTING
    myEnvironment.runSim(30000, 0, 1, 3000);

    myEnvironment.runSim(equilibsteps, 0);
    myEnvironment.runSim(measuresteps, 1);
    std::vector<double> outVals = myEnvironment.getParameters();
    std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5);

    lock.lock();
    appendToFile(myEnvironment.getOutputPath() + "\\results.dat", out);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Finished simulation in thread:  " << threadnum << std::endl;
    lock.unlock();
    std::cout << getCurrentTime().time_string << " [Thread Handler] Dumping magmoms from thread: " << threadnum << std::endl;
    myEnvironment.writeMagneticMomentsToFile("resultsMom");
    
}

void create_and_run_threads(std::string inputPath, std::string outputPath, std::vector<double> tempVec, int equilibsteps, int measurementsteps, std::vector<double> KVec, std::vector<double> CVec)
{
    std::vector<std::thread> threads;
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, 0);
    auto links = myEnvironment.getLinks();
    auto atomCoordinates = myEnvironment.getAtomCoordinates();
    auto linksNN = myEnvironment.getLinksNN();
    // Create and start threads
    for (int i = tempVec.size() - 1; i >= 0; --i) {
        for(int j = KVec.size() - 1; j >= 0; --j)
            for(int k = CVec.size() - 1; k >= 0; --k)
                threads.emplace_back(run_sim_thread, k + j * CVec.size() + i * KVec.size(), inputPath, outputPath, tempVec[i], equilibsteps, measurementsteps, links, atomCoordinates, linksNN, KVec[j], CVec[k]);
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


void run_sim_thread_wset(unsigned int threadnum, SettingsClass &simSettings, SimEnvironment &initEnvironment, double temperature, double KVec, double CVec, double HVec)
{
    std::unique_lock<std::mutex> lock(mtx);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup started for thread: " << threadnum << std::endl;

    auto inputPath = simSettings.get<std::string>("inputPath");
    auto outputPath = simSettings.get<std::string>("outputPath");

    int equilibsteps = simSettings.get<int>("equilibsteps");
    int measuresteps = simSettings.get<int>("measuresteps");

    std::vector<double> HVecIn = { 0.0, HVec, 0.0 };


    MomentReader magmomsFromFile;
    if (simSettings.get<int>("restartFromCheckpoint"))
    {
        magmomsFromFile = MomentReader(simSettings.get<std::string>("restartCheckpoint"));
    }
    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temperature, threadnum, initEnvironment.getLinks(), initEnvironment.getAtomCoordinates(), initEnvironment.getLinksNN(), {0.0, 0.0, KVec}, {CVec, CVec, 0.0}, magmomsFromFile.get_magmoms(), HVecIn);
    myEnvironment.setTemperature(temperature);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Setup finished for thread: " << threadnum << std::endl;
    std::cout << getCurrentTime().time_string << " [Thread Handler] Starting simulation in thread:  " << threadnum << std::endl;
    lock.unlock();

    myEnvironment.runSim(simSettings.get<int>("approachsteps"), 0, 1, simSettings.get<double>("initialtemp"));
    myEnvironment.runSim(equilibsteps, 0);
    myEnvironment.runSim(measuresteps, 1);

    std::vector<double> outVals = myEnvironment.getParameters();
    std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5) + "  " + formatDouble(outVals[6], 9, 5) + "  "  + formatDouble(outVals[7], 9, 5);

    lock.lock();
    appendToFile(myEnvironment.getOutputPath() + "\\results.dat", out);
    std::cout << getCurrentTime().time_string << " [Thread Handler] Finished simulation in thread:  " << threadnum << std::endl;
    lock.unlock();
    std::cout << getCurrentTime().time_string << " [Thread Handler] Dumping magmoms from thread: " << threadnum << std::endl;
    myEnvironment.writeMagneticMomentsToFile("resultsMom");

}


void create_and_run_threads(SettingsClass simSettings)
{
    std::vector<std::thread> threads;

    std::vector<double> tempVec = linspace(simSettings.get<double>("mintemp"), simSettings.get<double>("maxtemp"), simSettings.get<int>("tempspacing"));
    std::vector<double> KVec = linspace(simSettings.get<double>("SIAstart"), simSettings.get<double>("SIAend"), simSettings.get<int>("SIAspacing"));
    std::vector<double> CVec = linspace(simSettings.get<double>("CAstart"), simSettings.get<double>("CAend"), simSettings.get<int>("CAspacing"));
    std::vector<double> HVec = linspace(simSettings.get<double>("Hstart"), simSettings.get<double>("Hend"), simSettings.get<int>("Hspacing"));
    SimEnvironment initEnvironment = SimEnvironment(simSettings.get<std::string>("inputPath"), simSettings.get<std::string>("outputPath"), 0);
    // Create and start threads
    for (int i = tempVec.size() - 1; i >= 0; --i) {
        for (int j = KVec.size() - 1; j >= 0; --j)
            for (int k = CVec.size() - 1; k >= 0; --k)
                for (int l = HVec.size() - 1; l >= 0; --l)
                    threads.emplace_back(run_sim_thread_wset, l + k * HVec.size() + j * CVec.size() * HVec.size() + i * KVec.size() * CVec.size() * HVec.size(),
                                        std::ref(simSettings), std::ref(initEnvironment), tempVec[i], KVec[j], CVec[k], HVec[l]);
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
