#include <iostream>
#include <string>
#include <vector>

#include "helpers.h"
#include "SimEnvironment.h"

int main(int argc, char* argv[])
{
    std::string outputPath;
    std::string inputPath;
    double temperature;
    int equilibsteps, measuresteps;

    if (read_args(argc, argv, inputPath, outputPath, temperature, equilibsteps, measuresteps))
    {
        std::cerr << "Couldn't read arguments correctly, exiting!";
        return 1;
    }

    //SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temperature);
    //myEnvironment.runSim(equilibsteps, 0);
    //myEnvironment.runSim(measuresteps, 1);

    SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, 250);
    myEnvironment.runSim(20000, 0);
    myEnvironment.runSim(60000, 1);

    auto parameters = myEnvironment.getParameters();
    return 0;
}