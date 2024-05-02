
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
    //scriptExecutor.loadScript("inputs/ScriptEuRh2Si2.comps");
    scriptExecutor.loadScript("FeTeInputs/ScriptFeTe.comps");
    //scriptExecutor.loadScript("inputs/ScriptSFMO.comps");
    //scriptExecutor.loadScript("inputs/ScriptGdRu2Si2.comps");
    scriptExecutor.showScript();
    scriptExecutor.runScript();
    return 0;
}