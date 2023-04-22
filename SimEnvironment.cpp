#include "SimEnvironment.h"

std::vector<std::string> read_file_without_comments(const std::string& file_path) {
    std::vector<std::string> lines;
    std::ifstream input_file(file_path);

    if (input_file.is_open()) {
        std::string line;
        while (std::getline(input_file, line)) {
            if (line.empty() || line[0] != '#') {
                lines.push_back(line);
            }
        }
        input_file.close();
    }
    else {
        std::cerr << "Unable to open file: " << file_path << std::endl;
    }

    return lines;
}

std::tuple<int, int, double> extract_values(const std::string& line) {
    std::istringstream iss(line);
    int first, third;
    double last;
    std::vector<double> values;
    double value;

    while (iss >> value) {
        values.push_back(value);
    }

    first = static_cast<int>(values[0]);
    third = static_cast<int>(values[2]);
    last = values.back();

    return std::make_tuple(first, third, last);
}

SimEnvironment::SimEnvironment(std::string input, std::string output, double temp) : gen(rd())
{
	inputPath = input;
	outputPath = output;
	temperature = temp;
    std::cout << "[SimEnvironment] Setting up SimEnvironment!" << std::endl;
    std::cout << "[SimEnvironment] Reading input file: " << input << std::endl;
    auto lines = read_file_without_comments(inputPath);
    std::cout << "[SimEnvironment] Read input file, converting data!" << std::endl;
    for (int i = 0; i < lines.size(); i++)
    {
        auto values = extract_values(lines[i]);
        if (std::get<0>(values) > links.size())
            links.resize(std::get<0>(values));
        links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values) / 1000));
    }
    atomnum = links.size();
    std::cout << "[SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generate_random_vec());
    std::cout << "[SimEnvironment] Setup complete!" << std::endl << std::endl;
}

void SimEnvironment::runSim(int steps, bool measurement)
{
    if (measurement)
        std::cout << "[Kernel] Starting MC simulation in measurement mode for " << steps << " steps!" << std::endl;
    else
        std::cout << "[Kernel] Starting MC simulation in equilib mode for " << steps << " steps!" << std::endl;
    for (int step = 1; step <= steps; step++)
    {
        if (step % 100 == 0)
        {
            if (measurement)
                std::cout << "[Kernel] Currently at step " << step << " of " << steps << " in measurement mode!" << std::endl;
            else
                std::cout << "[Kernel] Currently at step " << step << " / " << steps << " in equilib mode!" << std::endl;
        }
        if (measurement)
        {
            magmomsHistory.push_back(magmoms);
            double M1 = 0.0;
            double M2 = 0.0;
            double M3 = 0.0;
            for (int i = 0; i < atomnum; i++)
            {
                M1 += magmoms[i][0];
                M2 += magmoms[i][1];
                M3 += magmoms[i][2];
            }
            M1 /= atomnum;
            M2 /= atomnum;
            M3 /= atomnum;
            meanmagmomsHistory.push_back(std::sqrt(M1 * M1 + M2 * M2 + M3 * M3));
        }

        double diff;
        double rate;
        std::uniform_real_distribution<> ran(0.0, 1.0);

        for (int i = 0; i < atomnum; i++)
        {
            auto randomVec = generate_random_vec();
            diff = energy_calculator(i, magmoms[i], randomVec);
            rate = std::exp(-diff / (kB * temperature));
            if (ran(gen) < rate)
                magmoms[i] = randomVec;
        }

    }
    if (measurement)
        std::cout << "[Kernel] Finished MC simulation in measurement mode!" << std::endl << std::endl;
    else
        std::cout << "[Kernel] Finished MC simulation in equilib mode!" << std::endl << std::endl;
}

double SimEnvironment::energy_calculator(int& index, std::vector<double>& oldMom, std::vector<double>& newMom)
{
    double Hold = 0.0;
    double Hnew = 0.0;
    int link;
    double param;
    for (auto linkforatom : links[index])
    {
        link = std::get<0>(linkforatom);
        param = std::get<1>(linkforatom);
        Hold -= dotProduct(oldMom, magmoms[link]) * param;
        Hnew -= dotProduct(newMom, magmoms[link]) * param;
    }

    return Hnew - Hold;
}

std::vector<double> SimEnvironment::generate_random_vec()
{
    double C1, C2, Csq;

    std::uniform_real_distribution<> ran(0.0, 1.0);

    Csq = 2;
    while (Csq >= 1)
    {
        C1 = 1.0 - 2.0 * ran(gen);
        C2 = 1.0 - 2.0 * ran(gen);
        Csq = C1 * C1 + C2 * C2;
    }
    std::vector<double> ranvec = { 2 * C1 * std::sqrt((1 - Csq)), 2 * C2 * std::sqrt((1 - Csq)), 1 - 2 * Csq };
    return ranvec;

}

std::vector<double> SimEnvironment::getParameters()
{
    double firstorder = 0.0;
    double secondorder = 0.0;
    double fourthorder = 0.0;
    for (double value : meanmagmomsHistory)
    {
        firstorder += value;
        secondorder += value * value;
        fourthorder += value * value * value * value;
    }
    firstorder /= meanmagmomsHistory.size();
    secondorder /= meanmagmomsHistory.size();
    fourthorder /= meanmagmomsHistory.size();

    double MagMom = firstorder;
    double Chi = atomnum / (kB * temperature) * (secondorder - firstorder * firstorder);
    double U4 = 1 - fourthorder / (3 * secondorder * secondorder);

    std::cout << "[SimEnvironment] Parameters are:" << std::endl;
    std::cout << "[SimEnvironment] Temperature: " << temperature << std::endl;
    std::cout << "[SimEnvironment] MagMom: " << MagMom << std::endl;
    std::cout << "[SimEnvironment] Chi: " << Chi << std::endl;
    std::cout << "[SimEnvironment] Kumulante of the fourth order: " << U4 << std::endl << std::endl;

    std::vector<double> returnVals = { temperature,MagMom,Chi,U4 };
    return returnVals;
}