
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
#include "simRunner.h"
#include "multithreadingSim.h"


// -e 2000 -m 2000 -t 50 -o C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\results -i C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\inputs\Jij_my_test.dat
int main(int argc, char* argv[])
{

    SimRunner runner("inputs/SETTINGS.cfg");

    // RUNNER TESTING
    runner.LOAD();
    runner.SETOUT("results/testSrCoO3.dat");
    runner.SETTEMP(1000.);
    runner.SETMAGFIELD(0.0, XYZ::X);

    double temp = 300;
    int steps = 10000;
    while (temp > 250)
    {
        runner.APPROACHTEMP(temp, steps);
        runner.EQUILIB(steps);
        runner.MEASUREMENT(steps);

        temp -= 5;
    }







    // END RUNNER TEST
    auto time = getCurrentTime();

    SettingsClass simSettings("inputs/SETTINGS.cfg");
    auto outputPath = simSettings.get<std::string>("outputPath");
    create_and_run_threads(simSettings);
    std::cout << getCurrentTime().time_string << " [Kernel] All threads finished. Sorting file: "  << outputPath + "\\results.dat" << std::endl;
    sort_file_by_first_entry(outputPath + "\\results.dat");

    return 0;
}