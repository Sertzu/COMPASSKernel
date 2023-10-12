
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
#include "scriptHandler.h"

// -e 2000 -m 2000 -t 50 -o C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\results -i C:\Users\Sertzu\source\repos\COMPASSKernel\COMPASSKernel\inputs\Jij_my_test.dat
int main(int argc, char* argv[])
{
    ScriptHandler scriptExecutor("inputs/SETTINGS.cfg");
    scriptExecutor.loadScript("inputs/testScriptBase.comps");
    scriptExecutor.showScript();
    scriptExecutor.runScript();


    exit(1);
    SimRunner simRunner("inputs/SETTINGS.cfg");
    
    int stepsInit = 40000;
    int steps = 1000;
    // RUNNER TESTING
    simRunner.LOAD("inputs/GdMn6Sn6_b20_c0.1_mJValue.dat");
    simRunner.SETOUT("results/GdMn6Sn6_b20_c0.1_mJValue.dat");
    simRunner.SETTEMP(1000.);
    simRunner.SETMAGFIELD(0.1, XYZ::Y);

    simRunner.APPROACHTEMP(1., stepsInit);
    simRunner.EQUILIB(stepsInit);
    simRunner.MEASUREMENT(stepsInit);

    simRunner.SWEEPMAG(100., 10., steps, steps, steps);
    simRunner.SWEEPMAG(-100., -10., steps, steps, steps);
    simRunner.SWEEPMAG(100., 10., steps, steps, steps);

    //simRunner.SWEEPTEMP(250., 5., 20000, 20000, 20000);
    
    exit(1);

    // END RUNNER TEST
    auto time = getCurrentTime();

    SettingsClass simSettings("inputs/SETTINGS.cfg");
    auto outputPath = simSettings.get<std::string>("outputPath");
    create_and_run_threads(simSettings);
    std::cout << getCurrentTime().time_string << " [Kernel] All threads finished. Sorting file: "  << outputPath + "\\results.dat" << std::endl;
    sort_file_by_first_entry(outputPath + "\\results.dat");

    return 0;
}