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

std::tuple<int, int, double, bool, bool, bool> extractLinks(const std::string& line) {
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
    last = values[8];

    bool isNearestNeighbourX = false;
    bool isNearestNeighbourY = false;
    bool isNearestNeighbourZ = false;
    /*
    if((abs(static_cast<int>(values[12])) + abs(static_cast<int>(values[13])) + abs(static_cast<int>(values[14]))) < 2)
    {
        if (abs(static_cast<int>(values[12])) == 1)
            isNearestNeighbourX = true;
        if (abs(static_cast<int>(values[13])) == 1)
            isNearestNeighbourY = true;
        if (abs(static_cast<int>(values[14])) == 1)
            isNearestNeighbourZ = true;
    }
    */
    return std::make_tuple(first, third, last, isNearestNeighbourX, isNearestNeighbourY, isNearestNeighbourZ);
}


std::tuple<int, std::vector<double>> extractCoordinates(const std::string& line)
{
    std::istringstream iss(line);
    int first, third;
    double last;
    std::vector<double> values;
    double value;
    std::vector<double> coords;
    while (iss >> value) {
        values.push_back(value);
    }

    first = static_cast<int>(values[0]);
    coords.push_back(values[9]);
    coords.push_back(values[10]);
    coords.push_back(values[11]);

    return std::make_tuple(first, coords);
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

    atomnum = std::get<0>(extractLinks(lines.back()));
    linksNN.resize(3);
    links.resize(atomnum);
    linksNN[0].resize(atomnum);
    linksNN[1].resize(atomnum);
    linksNN[2].resize(atomnum);
    atomCoordinates.resize(atomnum);

    for (int i = 0; i < lines.size(); i++)
    {
        auto values = extractLinks(lines[i]);
        auto coords = extractCoordinates(lines[i]);
        links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values)));
        if(std::get<3>(values))
        {
            linksNN[0][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<4>(values))
        {
            linksNN[1][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<5>(values))
        {
            linksNN[2][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }

        atomCoordinates[std::get<0>(coords) - 1] = coords;

    }
    interactionK = { 0.0, 0.0, 0.0 };
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generateRandomVecSingle());
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

    atomnum = std::get<0>(extractLinks(lines.back()));
    linksNN.resize(3);
    links.resize(atomnum);
    linksNN[0].resize(atomnum);
    linksNN[1].resize(atomnum);
    linksNN[2].resize(atomnum);
    atomCoordinates.resize(atomnum);

    for (int i = 0; i < lines.size(); i++)
    {
        auto values = extractLinks(lines[i]);
        auto coords = extractCoordinates(lines[i]);
        links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values)));
        if (std::get<3>(values))
        {
            linksNN[0][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<4>(values))
        {
            linksNN[1][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<5>(values))
        {
            linksNN[2][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }

        atomCoordinates[std::get<0>(coords) - 1] = coords;
    }
    interactionK = { 0.0, 0.0, 0.0 };
    atomnum = links.size();
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generateRandomVecSingle());
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

SimEnvironment::SimEnvironment(std::string input, std::string output, double temp, unsigned int threadnumber, const std::vector<std::vector<std::tuple<int, double>>> &linksIn, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinatesIn,const std::vector<std::vector<std::vector<int>>>& linksNNIn, std::vector<double> interactionKZIn, std::vector<double> interactionCIn, std::vector<std::vector<double>> magmomsIn, std::vector<double> magneticFieldHIn)
{
    threadnum = threadnumber;
    if (threadnum)
        tlog = ThreadLogger(threadnum, false);
    else
        tlog = ThreadLogger(threadnum, true);
    inputPath = input;
    outputPath = output;
    temperature = temp;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setting up SimEnvironment!" << std::endl;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Using links from reference!" << input << std::endl;
    links = linksIn;
    linksNN = linksNNIn;
    interactionK = interactionKZIn;
    interactionC = interactionCIn;
    magneticFieldH = magneticFieldHIn;
    atomCoordinates = atomCoordinatesIn;
    atomnum = links.size();
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;

    if (magmomsIn.size() == 0)
    {
        for (int i = 0; i < atomnum; i++)
            magmoms.push_back(generateRandomVecSingle());
    }
    else
        magmoms = magmomsIn;
    if ((abs(interactionC[0]) + abs(interactionC[0]) + abs(interactionC[0])) < 1e-8)
        enableCompassAnisotropy = true;
    else
        enableCompassAnisotropy = false;
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

std::vector<std::vector<std::tuple<int, double>>> SimEnvironment::getLinks()
{
    return links;
}

std::vector<std::vector<std::vector<int>>> SimEnvironment::getLinksNN()
{
    return linksNN;
}

std::vector<std::tuple<int, std::vector<double>>> SimEnvironment::getAtomCoordinates()
{
    return atomCoordinates;
}

std::vector<std::vector<double>> SimEnvironment::getMagmoms()
{
    return magmoms;
}

void SimEnvironment::runSim(int steps, bool measurement, bool approach, double initialtemp)
{
    magmomsHistory.clear();
    meanmagmomsHistory.clear();
    energyHistory.clear();
    if (measurement)
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in measurement mode for " << steps << " steps!" << std::endl;
    else if(approach)
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in approach mode for " << steps << " steps!" << std::endl;
    else
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in equilib mode for " << steps << " steps!" << std::endl;
    auto lasttime = getCurrentTime();
    double approachFactor;
    double runningTemperature = temperature;
    if (approach)
    {
        runningTemperature = initialtemp;
        approachFactor = std::pow(temperature / initialtemp, 1.0 / steps);
        tlog << getCurrentTime().time_string << " [Kernel] Initial temperature is " << initialtemp << "K and approachFactor is " << approachFactor << std::endl;
    }

    std::vector<std::vector<double>> randomNumberVector3D(atomnum, { 0.0 ,0.0 , 0.0});
    std::vector<double> rateVector(atomnum, 0.0);

    BS::thread_pool pool;

    for (int step = 1; step <= steps; step++)
    {
        if (step % 1000 == 0)
        {
            auto currenttime = getCurrentTime();
            if (measurement)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " of " << steps << " in measurement mode. Last 100 steps took: " 
                << currenttime.unix_time - lasttime.unix_time << " seconds!" << std::endl;
            else if (approach)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in approach mode. Current temperature is " 
                << runningTemperature << "K !" << std::endl;
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
            if (steps - step < MAXMAGMOMHISTSIZE)
                magmomsHistory.push_back(magmoms);
        }

        double diff;
        double rate;
        std::uniform_real_distribution<> ran(0.0, 1.0);
        double beta = 1 / (kB * runningTemperature);
        generateRandomVecArray(randomNumberVector3D);
        /*
        //Normal Way
        for (int i = 0; i < atomnum; i++)
        {
            diff = energy_diff_calculator(i, magmoms[i], randomNumberVector[i]);
            rate = std::exp(-diff * beta);
            if (ran(gen) < rate)
            {
                magmoms[i] = randomNumberVector[i];
            }
        }
        */
        
        //Parallel Way
        pool.push_loop(atomnum,
            [&](const int a, const int b)
            {
                for (int i = a; i < b; ++i)
                    rateVector[i] = rate_calculator(i, beta, magmoms[i], randomNumberVector3D[i]);
            });

        pool.wait_for_tasks();

        pool.push_loop(atomnum,
            [&](const int a, const int b)
            {
                for (int i = a; i < b; ++i)
                    if (ran(gen) < rateVector[i])
                    {
                        magmoms[i] = randomNumberVector3D[i];
                    }
            });

        pool.wait_for_tasks();
        if (approach)
            runningTemperature *= approachFactor;

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
    double param, sum;
    for (auto linkforatom : links[index])
    {
        link = std::get<0>(linkforatom);
        param = std::get<1>(linkforatom);
        dotProduct(sum, oldMom, magmoms[link]);
        Hold -= sum * param;
        dotProduct(sum, newMom, magmoms[link]);
        Hnew -= sum * param;
    }
    // Single Ion anisotropy + Zeeman Term
    for (int i = 0; i < 3; i++)
    {
        Hold -= oldMom[i] * oldMom[i] * interactionK[i];
        Hnew -= newMom[i] * newMom[i] * interactionK[i];

        Hold -= oldMom[i] * magneticFieldH[i];
        Hnew -= newMom[i] * magneticFieldH[i];
    }
    // Compass anisotropy
    if(enableCompassAnisotropy)
        for(int i = 0; i < 3; i++)
            for (auto linkforatom : linksNN[i][index])
            {
                Hold -= oldMom[i] * magmoms[linkforatom][i] * interactionC[i];
                Hnew -= newMom[i] * magmoms[linkforatom][i] * interactionC[i];
            }

    return Hnew - Hold;
}

double SimEnvironment::rate_calculator(int& index, double& beta, std::vector<double>& oldMom, std::vector<double>& newMom)
{
    return std::exp(-energy_diff_calculator(index, oldMom, newMom) * beta);
}

double SimEnvironment::energy_calculator()
{
    double energy = 0.0;
    int link;
    double param, sum;
    for(int i = 0; i < atomnum; i++)
        for (auto linkforatom : links[i])
        {
            link = std::get<0>(linkforatom);
            param = std::get<1>(linkforatom);
            dotProduct(sum, magmoms[i], magmoms[link]);
            energy -= sum * param;
        }

    return energy / 2;
}

std::vector<double> SimEnvironment::generateRandomVecSingle()
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
    double A, B, C;
    A = 2 * C1 * std::sqrt((1 - Csq));
    B = 2 * C2 * std::sqrt((1 - Csq));
    C = 1 - 2 * Csq;
    std::vector<double> ranvec = { A, B, C};
    return ranvec;

}

void SimEnvironment::generateRandomVecArray(std::vector<std::vector<double>>& vecIn)
{
    double C1, C2, Csq;

    std::uniform_real_distribution<> ran(0.0, 1.0);
    for (int i = 0; i < vecIn.size(); i++)
    {
        Csq = 2;
        while (Csq >= 1)
        {
            C1 = 1.0 - 2.0 * ran(gen);
            C2 = 1.0 - 2.0 * ran(gen);
            Csq = C1 * C1 + C2 * C2;
        }
        vecIn[i][0] = 2 * C1 * std::sqrt((1 - Csq));
        vecIn[i][1] = 2 * C2 * std::sqrt((1 - Csq));
        vecIn[i][2] = 1 - 2 * Csq;
    }
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

    auto magtotal = magneticFieldH[0] + magneticFieldH[1] + magneticFieldH[2];

    std::vector<double> returnVals = { temperature, MagMom, Chi, U4, E, HeatCapacity, interactionK[2], magtotal };
    return returnVals;
}

void SimEnvironment::writeMagneticMomentsToFile()
{
    // Open the output file
    std::ostringstream KString;
    KString << std::scientific << std::setprecision(1) << interactionK[2];
    std::string interactionKString = KString.str();

    std::ostringstream CString;
    CString << std::scientific << std::setprecision(1) << interactionC[0];
    std::string interactionCString = CString.str();

    std::ostringstream HString;
    HString << std::scientific << std::setprecision(1) << magneticFieldH[0] + magneticFieldH[1] + magneticFieldH[2];
    std::string magneticFieldHString = HString.str();

    std::ofstream outFile(joinPaths(outputPath, "endMagmoms_temp" + std::to_string(temperature) + "_K" + interactionKString + "_C" + interactionCString + "_H" + magneticFieldHString + ".magres"));
    if (!outFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    // Iterate through both vectors and write the entries to the file
    for (size_t i = 0; i < atomCoordinates.size(); ++i) {
        for (const auto& value : std::get<1>(atomCoordinates[i])) {
            outFile << formatDouble(value, 9, 5) << ' ';
        }
        double sumx, sumy, sumz;
        double sumxnorm, sumynorm, sumznorm;
        sumx = sumy = sumz = 0;
        for (int j = 0; j < magmomsHistory.size(); j++) {
            sumx += magmomsHistory[j][i][0];
            sumy += magmomsHistory[j][i][1];
            sumz += magmomsHistory[j][i][2];
        }
        double norm = std::sqrt(sumx * sumx + sumy * sumy + sumz * sumz);

        outFile << formatDouble(sumx / norm, 9, 5) << ' ';
        outFile << formatDouble(sumy / norm, 9, 5) << ' ';
        outFile << formatDouble(sumz / norm, 9, 5) << ' ';

        outFile << formatDouble(sumx, 11, 5) << ' ';
        outFile << formatDouble(sumy, 11, 5) << ' ';
        outFile << formatDouble(sumz, 11, 5) << ' ';
        outFile << '\n';
    }

    // Close the output file
    outFile.close();
}
