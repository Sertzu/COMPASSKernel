#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#include "helpers.h"
#include "MomentReader.h"
#include "SettingsClass.h"
#include "timeUtility.h"
#include "SimEnvironment.h"
#include "multithreadingSim.h"


// -e 2000 -m 2000 -t 50 -o C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\results -i C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\inputs\Jij_my_test.dat
int main(int argc, char* argv[])
{

    SettingsClass simSettings("inputs\\SETTINGS.cfg");
    auto time = getCurrentTime();

    if(!simSettings.get<int>("multithreading"))
    {
        /*
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
        */
    }
    else
    {
        auto outputPath = simSettings.get<std::string>("outputPath");
        create_and_run_threads(simSettings);
        std::cout << getCurrentTime().time_string << " [Kernel] All threads finished. Sorting file: "  << outputPath + "\\results.dat" << std::endl;
        sort_file_by_first_entry(outputPath + "\\results.dat");
    }

    return 0;
}