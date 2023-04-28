#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#include "helpers.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include "multithreadingSim.h"


// -e 2000 -m 2000 -t 50 -o C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\results -i C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\inputs\Jij_my_test.dat
int main(int argc, char* argv[])
{
    std::string outputPath;
    std::string inputPath;
    double temperature;
    int equilibsteps, measuresteps;
    double mintemp = 200;
    double tempstep = 20;
    bool multithreading = 1;
    auto time = getCurrentTime();
    if(argc > 1)
    {
        if (read_args(argc, argv, inputPath, outputPath, temperature, equilibsteps, measuresteps))
        {
            std::cerr << "Couldn't read arguments correctly, exiting!";
            return 1;
        }
    }
    else
    {
        outputPath = "results";
        //inputPath = "inputs\\Jij_26_17_cor.dat";
        inputPath = "inputs\\Jij_my_test_10.dat";
        temperature = 400;
        equilibsteps = 20000;
        measuresteps = 20000;
    }
    if (argc > 1)
    {
        SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temperature);
        myEnvironment.runSim(equilibsteps, 0);
        myEnvironment.runSim(measuresteps, 1);
        std::vector<double> outVals = myEnvironment.getParameters();
        std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5);

        appendToFile(outputPath + "\\results.dat", out);
    }
    else if(!multithreading)
    {
        SimEnvironment myEnvironment = SimEnvironment(inputPath, outputPath, temperature);
        while(temperature >= mintemp)
        {
            myEnvironment.setTemperature(temperature);
            std::cout << getCurrentTime().time_string << " [Kernel] Doing calculation for temperature: " << temperature << "K" << std::endl;
            myEnvironment.runSim(equilibsteps, 0);
            myEnvironment.runSim(measuresteps, 1);
            std::vector<double> outVals = myEnvironment.getParameters();
            std::string out = formatDouble(outVals[0], 9, 5) + "  " + formatDouble(outVals[1], 9, 5) + "  " + formatDouble(outVals[2], 9, 5) + "  " + formatDouble(outVals[3], 9, 5) + "  " + formatDouble(outVals[4], 9, 5) + "  " + formatDouble(outVals[5], 9, 5);

            appendToFile(outputPath + "\\results.dat", out);
            std::cout << getCurrentTime().time_string << " [Kernel] Reducing temperature by " << tempstep << "K";
            temperature -= tempstep;
            std::cout << " || New temperature is: " << temperature << "K" << std::endl;
        }
    }
    else
    {
        std::vector<double> tempVec;
        while (temperature >= mintemp)
        {
            tempVec.push_back(temperature);
            temperature -= tempstep;
        }
        create_and_run_threads(inputPath, outputPath, tempVec, equilibsteps, measuresteps);
        std::cout << getCurrentTime().time_string << " [Kernel] All threads finished. Sorting file: "  << outputPath + "\\results.dat" << std::endl;
        sort_file_by_first_entry(outputPath + "\\results.dat");
    }

    return 0;
}