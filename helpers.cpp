#include "helpers.h"


static void show_usage(std::string name)
{
    std::cerr << "Usage of this program:" << std::endl
        << "Options:\n"
        << "\t-h,--help\t\tShow this help message.\n"
        << "\t-o,--output\t\tSpecify the output path.\n"
        << "\t-i,--input\t\tSpecify the input file.\n"
        << "\t-t,--temperature\t\tSpecify the temperature.\n"
        << "\t-e,--equilibsteps\t\tSpecify the equilibsteps.\n"
        << "\t-m,--measuresteps\t\tSpecify the measuresteps.\n\n"
        << "NOTE: ALL OF THESE EXCEPT --help HAVE TO BE SPECIFIED!!!!\n"
        << std::endl;
}

int read_args(int argc, char* argv[], std::string& input, std::string& output, double& temperature, int& equilibsteps, int& measuresteps)
{
    if (argc < 3) {
        show_usage(argv[0]);
        return 1;
    }
    temperature = -1;
    equilibsteps = -1;
    measuresteps = -1;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        }
        else if ((arg == "-i") || (arg == "--input")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                std::cout << "Input given is: " << argv[++i] << std::endl;
                input = argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
            }
            else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--output option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-o") || (arg == "--output")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                std::cout << "Destination given is: " << argv[++i] << std::endl;
                output = argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
            }
            else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--output option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-t") || (arg == "--temperature")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                std::cout << "Temperature given is: " << argv[++i] << std::endl;
                temperature = std::stod(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
            }
            else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--temperature option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-e") || (arg == "--equilibsteps")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                std::cout << "Equilibsteps given is: " << argv[++i] << std::endl;
                equilibsteps = std::stoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
            }
            else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--equilibsteps option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-m") || (arg == "--measuresteps")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                std::cout << "Measuresteps given is: " << argv[++i] << std::endl;
                measuresteps = std::stoi(argv[i]); // Increment 'i' so we don't get the argument as the next argv[i].
            }
            else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--measuresteps option requires one argument." << std::endl;
                return 1;
            }
        }
    }
    if (equilibsteps < 0)
    {
        std::cerr << "--equilibsteps option not given or invalid." << std::endl;
        return 2;
    }
    if (measuresteps < 0)
    {
        std::cerr << "--measuresteps option not given or invalid." << std::endl;
        return 2;
    }
    if (temperature < 0)
    {
        std::cerr << "--temperature option not given or invalid." << std::endl;
        return 2;
    }
    if (!input.size())
    {
        std::cerr << "--input option not given or invalid." << std::endl;
        return 2;
    }
    if (!output.size())
    {
        std::cerr << "--output option not given or invalid." << std::endl;
        return 2;
    }
    return 0;
}


double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size() || vec1.size() != 3) {
        throw std::invalid_argument("Both vectors must be of size 3.");
    }

    double result = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }

    return result;
}