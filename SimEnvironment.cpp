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
    threadnum = 0;
    tlog = ThreadLogger(threadnum, true);
    inputPath = input;
    outputPath = output;
    temperature = temp;
    tlog << getCurrentTime().time_string <<  " [SimEnvironment] Setting up SimEnvironment!" << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Reading input file: " << input << std::endl;
    auto lines = read_file_without_comments(inputPath);
    tlog << getCurrentTime().time_string << " [SimEnvironment] Read input file, converting data!" << std::endl;
    for (int i = 0; i < lines.size(); i++)
    {
        auto values = extract_values(lines[i]);
        if (std::get<0>(values) > links.size())
            links.resize(std::get<0>(values));
        links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values) / 1000));
    }
    atomnum = links.size();
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generate_random_vec());
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

SimEnvironment::SimEnvironment(std::string input, std::string output, double temp, unsigned int threadnumber) : gen(rd())
{
    threadnum = threadnumber;
    if(threadnum)
        tlog = ThreadLogger(threadnum, false);
    else
        tlog = ThreadLogger(threadnum, true);
    inputPath = input;
    outputPath = output;
    temperature = temp;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setting up SimEnvironment!" << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Reading input file: " << input << std::endl;
    auto lines = read_file_without_comments(inputPath);
    tlog << getCurrentTime().time_string << " [SimEnvironment] Read input file, converting data!" << std::endl;
    for (int i = 0; i < lines.size(); i++)
    {
        auto values = extract_values(lines[i]);
        if (std::get<0>(values) > links.size())
            links.resize(std::get<0>(values));
        links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values) / 1000));
    }
    atomnum = links.size();
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generate_random_vec());
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

void SimEnvironment::runSim(int steps, bool measurement)
{
    magmomsHistory.clear();
    meanmagmomsHistory.clear();
    energyHistory.clear();
    if (measurement)
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in measurement mode for " << steps << " steps!" << std::endl;
    else
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in equilib mode for " << steps << " steps!" << std::endl;
    auto lasttime = getCurrentTime();
    for (int step = 1; step <= steps; step++)
    {
        if (step % 100 == 0)
        {
            auto currenttime = getCurrentTime();
            if (measurement)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " of " << steps << " in measurement mode. Last 100 steps took: " 
                << currenttime.unix_time - lasttime.unix_time << " seconds!" << std::endl;
            else
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in equilib mode. Last 100 steps took: " 
                << currenttime.unix_time - lasttime.unix_time << " seconds!" << std::endl;
            lasttime = currenttime;
        }
        if (measurement)
        {
            //magmomsHistory.push_back(magmoms);
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
            energyHistory.push_back(energy_calculator());
        }

        double diff;
        double rate;
        std::uniform_real_distribution<> ran(0.0, 1.0);

        for (int i = 0; i < atomnum; i++)
        {
            auto randomVec = generate_random_vec();
            diff = energy_diff_calculator(i, magmoms[i], randomVec);
            rate = std::exp(-diff / (kB * temperature));
            if (ran(gen) < rate)
                magmoms[i] = randomVec;
        }

    }
    if (measurement)
        tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in measurement mode!" << std::endl << std::endl;
    else
        tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in equilib mode!" << std::endl << std::endl;
}

void SimEnvironment::setTemperature(double temp)
{
    temperature = temp;
}

std::string SimEnvironment::getOutputPath()
{
    return outputPath;
}

double SimEnvironment::energy_diff_calculator(int& index, std::vector<double>& oldMom, std::vector<double>& newMom)
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

double SimEnvironment::energy_calculator()
{
    double energy = 0.0;
    int link;
    double param;
    for(int i = 0; i < atomnum; i++)
        for (auto linkforatom : links[i])
        {
            link = std::get<0>(linkforatom);
            param = std::get<1>(linkforatom);
            energy -= dotProduct(magmoms[i], magmoms[link]) * param;
        }

    return energy / 2;
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
    double magmomFirstorder = 0.0;
    double magmomSecondorder = 0.0;
    double magmomFourthorder = 0.0;
    for (double value : meanmagmomsHistory)
    {
        magmomFirstorder += value;
        magmomSecondorder += value * value;
        magmomFourthorder += value * value * value * value;
    }
    magmomFirstorder /= meanmagmomsHistory.size();
    magmomSecondorder /= meanmagmomsHistory.size();
    magmomFourthorder /= meanmagmomsHistory.size();

    double energyFirstorder = 0.0;
    double energySecondorder = 0.0;
    for (double value : energyHistory)
    {
        energyFirstorder += value / atomnum;
        energySecondorder += value * value / (atomnum * atomnum);
    }
    energyFirstorder /= energyHistory.size();
    energySecondorder /= energyHistory.size();

    double MagMom = magmomFirstorder;
    double Chi = atomnum / (kB * temperature) * (magmomSecondorder - magmomFirstorder * magmomFirstorder);
    double U4 = 1 - magmomFourthorder / (3 * magmomSecondorder * magmomSecondorder);
    double E = energyFirstorder;
    double HeatCapacity = atomnum * atomnum / (kB * temperature * temperature) * (energySecondorder - energyFirstorder * energyFirstorder);

    tlog << getCurrentTime().time_string << " [SimEnvironment] Parameters are:" << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Temperature: " << temperature << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] MagMom: " << MagMom << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Chi: " << Chi << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Energy: " << E << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] HeatCapacity: " << HeatCapacity << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Kumulante of the fourth order: " << U4 << std::endl << std::endl;

    std::vector<double> returnVals = { temperature,MagMom,Chi,U4 ,E, HeatCapacity };
    return returnVals;
}
