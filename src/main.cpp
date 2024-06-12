#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <argparse.hpp>

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
    argparse::ArgumentParser program("compasskernel");

    program.add_argument("Script Path")
        .help("Specify the path to the simulation script!")
        .required();  

    std::string path;
    try {
        program.parse_args(argc, argv);

        path = program.get<std::string>("Script Path");

        std::cout << "Provided script path: " << path << std::endl;
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    ScriptHandler scriptExecutor;

    scriptExecutor.loadScript(path);
    scriptExecutor.showScript();
    scriptExecutor.runScript();

    return 0;
}