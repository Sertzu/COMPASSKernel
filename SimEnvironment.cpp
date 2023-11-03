#include "SimEnvironment.h"
#include <filesystem> // Requires C++17

namespace fs = std::filesystem;

inline void dotProduct(double& sum, const std::vector<double>& vec1, const std::vector<double>& vec2) {
    sum = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

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
    
    if(values.size() == 15)
        if((abs(static_cast<int>(values[12])) + abs(static_cast<int>(values[13])) + abs(static_cast<int>(values[14]))) < 2)
        {
            if (abs(static_cast<int>(values[12])) == 1)
                isNearestNeighbourX = true;
            if (abs(static_cast<int>(values[13])) == 1)
                isNearestNeighbourY = true;
            if (abs(static_cast<int>(values[14])) == 1)
                isNearestNeighbourZ = true;
        }
    
    return std::make_tuple(first, third, last, isNearestNeighbourX, isNearestNeighbourY, isNearestNeighbourZ);
}


std::tuple<int, std::vector<double>> extractCoordinates(const std::string& line)
{
    std::istringstream iss(line);
    int first, third, counter;
    double last;
    static std::vector<double> values(20, 0.0);
    double value;
    static std::vector<double> coords{0.0,0.0,0.0};
    counter = 0;
    while (iss >> value && counter < 12) {
        values[counter] = value;
        counter++;
    }

    first = static_cast<int>(values[0]);
    coords[0] = values[9];
    coords[1] = values[10];
    coords[2] = values[11];

    return std::make_tuple(first, coords);
}

SimEnvironment::SimEnvironment()
{
    inputPath = "";
    outputPath = "";
    temperature = 0;
}

SimEnvironment::SimEnvironment(std::string input, std::string output, double temp) : gen(rd())
{
    threadnum = 0;
    tlog = ThreadLogger(threadnum, true);
    inputPath = input;
    outputPath = output;
    temperature = temp;

    if (inputPath == "")
        return;
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
    singleIonAnisotropyTerm = { 0.0, 0.0, 0.0 };
    compassAnisotropyTerm = { 0.0, 0.0, 0.0 };
    zeemanTerm = { 0.0, 0.0, 0.0 };
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    magmoms.reserve(atomnum);
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
    singleIonAnisotropyTerm = { 0.0, 0.0, 0.0 };
    atomnum = links.size();
    tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    for (int i = 0; i < atomnum; i++)
        magmoms.push_back(generateRandomVecSingle());
    tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

SimEnvironment::SimEnvironment(std::string input, std::string output, double temp, unsigned int threadnumber, const std::vector<std::vector<std::tuple<int, double>>> &linksIn, const std::vector<std::tuple<int, std::vector<double>>>& atomCoordinatesIn,const std::vector<std::vector<std::vector<int>>>& linksNNIn, std::vector<double> interactionCompIn, std::vector<double> interactionCIn, std::vector<std::vector<double>> magmomsIn, std::vector<double> magneticFieldHIn)
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
    singleIonAnisotropyTerm = interactionCompIn;
    compassAnisotropyTerm = interactionCIn;
    zeemanTerm = magneticFieldHIn;
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
    if ((abs(compassAnisotropyTerm[0]) + abs(compassAnisotropyTerm[0]) + abs(compassAnisotropyTerm[0])) < 1e-8)
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

void SimEnvironment::runSim(int steps, bool measurement, bool approachTemp, double initialtemp, bool approachMag, double initialH)
{
    magmomsHistory.clear();
    meanmagmomsHistory.clear();
    energyHistory.clear();
    if (measurement)
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in measurement mode for " << steps << " steps!" << std::endl;
    else if(approachTemp)
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in approach mode for " << steps << " steps!" << std::endl;
    else
        tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in equilib mode for " << steps << " steps!" << std::endl;
    auto lasttime = std::chrono::high_resolution_clock::now();
    double approachValue;
    double runningTemperature = temperature;
    double runningMag = zeemanTerm[magDir];
    if (approachTemp)
    {
        runningTemperature = initialtemp;
        approachValue = std::pow(temperature / initialtemp, 1.0 / steps);
        tlog << getCurrentTime().time_string << " [Kernel] Initial temperature is " << initialtemp << "K and Approachfactor is " << approachValue << std::endl;
    }
    else if (approachMag)
    {
        runningMag = initialH;
        approachValue = (zeemanTerm[magDir] - initialH) / steps;
        zeemanTerm[magDir] = runningMag;
        tlog << getCurrentTime().time_string << " [Kernel] Initial Magnetic Field is " << initialH << "T and Approachvalue is " << approachValue << std::endl;
    }

    std::vector<std::vector<double>> randomNumberVector3D(atomnum, { 0.0 ,0.0 , 0.0});
    std::vector<double> rateVector(atomnum, 0.0);
    std::vector<double> acceptanceVector(atomnum, 0.0);

    //generateRandomVecArray(randomNumberVector3D);
    //generateAcceptanceVec(acceptanceVector);

    std::vector<std::vector<double>> randomNumberVector3DSwap(atomnum, { 0.0 ,0.0 , 0.0 });
    std::vector<double> acceptanceVectorSwap(atomnum, 0.0);

    BS::thread_pool pool;
    magmomsHistory.clear();
    meanmagmomsHistory.clear();
    energyHistory.clear();

    meanmagmomsHistory.resize(steps);
    energyHistory.resize(steps);
    std::future<double> energy;

    int workerCount = 12;

    std::vector<std::vector<std::vector<double>>> magmomCopy;
    for (int i = 0; i < workerCount; i++)
    {
        magmomCopy.push_back(magmoms);
    }
    auto startIndices = getStartIndices(magmoms.size(), workerCount);

    std::barrier syncPointInit(workerCount + 1);
    std::barrier syncPointRun(workerCount + 1);

    auto mcWorker = [&](std::stop_token stopToken, double& beta, std::vector<std::vector<double>> &magmoms, std::vector<std::vector<double>>& randomNumberVector3D, std::vector<double> &acceptanceVector, int start, int end)
    {
        std::uniform_real_distribution<> randomizer(0.0, 1.0);
        std::random_device randDevice;
        std::mt19937 generator(randDevice());
        auto atomlinks = links;
        double rate = 0.0;
        double acceptor = 0.0;
        std::vector<std::vector<double>> randomVec = { {0.0,0.0,0.0} };
        generateRandomVecArray(randomVec, randomizer, generator);

        syncPointInit.arrive_and_wait();
        while (!stopToken.stop_requested())
        {
            for (int i = start; i < end; ++i)
            {
                generateRandomVecArray(randomVec, randomizer, generator);
                acceptor = randomizer(generator);
                rate = rate_calculator(i, beta, magmoms[i], randomVec[0], magmoms, atomlinks);
                if (acceptor < rate)
                {
                    magmoms[i][0] = randomVec[0][0];
                    magmoms[i][1] = randomVec[0][1];
                    magmoms[i][2] = randomVec[0][2];
                }
            }
            syncPointRun.arrive_and_wait();
            for (auto& magmomEntry : magmomCopy)
                for (int i = start; i < end; ++i)
                {
                    magmomEntry[i][0] = magmoms[i][0];
                    magmomEntry[i][1] = magmoms[i][1];
                    magmomEntry[i][2] = magmoms[i][2];
                }
            syncPointRun.arrive_and_wait();
        }
    };

    std::vector<std::jthread> mcWorkers;
    double beta = 1 / (kB * runningTemperature);

    mcWorkers.reserve(workerCount);
    for (int i = 0; i < workerCount; ++i) {
        int start = startIndices[i];
        int end = (i == workerCount - 1) ? magmoms.size() : startIndices[i + 1];
        mcWorkers.emplace_back(mcWorker, std::ref(beta), std::ref(magmomCopy[i]), std::ref(randomNumberVector3D),std::ref(acceptanceVector), start, end);
    }

    for (int step = 1; step <= steps; step++)
    {
        if (step % statusSteps == 0)
        {
            auto currenttime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currenttime - lasttime);
            if (measurement)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " of " << steps << " in measurement mode. Last " << statusSteps << " steps took: "
                << duration.count() << " ms!" << std::endl;
            else if (approachTemp)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in approachTemp mode. Current temperature is " 
                << runningTemperature << "K !" << " Last " << statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            else if (approachMag)
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in approachMag mode. Current Magnetic Field is "
                << zeemanTerm[magDir] << "T !" << " Last " << statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            else
                tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in equilib mode. Last " << statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            lasttime = currenttime;
        }
        if (measurement)
        {
            
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
            meanmagmomsHistory[step - 1] = std::sqrt(M1* M1 + M2 * M2 + M3 * M3);
            energy = pool.submit(&SimEnvironment::energy_calculator, this);
            //energyHistory[step - 1] = energy_calculator();
            if (steps - step < MAXMAGMOMHISTSIZE)
                magmomsHistory.push_back(magmoms);
        }

        if (step == 1)
        {
            if (measurement)
                energyHistory[step - 1] = energy.get();
            syncPointInit.arrive_and_wait();
        }

        if (measurement && step != 1)
            energyHistory[step - 1] = energy.get();

        //generateRandomVecArray(randomNumberVector3DSwap);
        //generateAcceptanceVec(acceptanceVectorSwap);

        syncPointRun.arrive_and_wait();

        //randomNumberVector3D.swap(randomNumberVector3DSwap);
        //acceptanceVector.swap(acceptanceVectorSwap);

        /*
        //Parallel Way
        pool.push_loop(atomnum,
            [&](const int a, const int b)
            {
                for (int i = a; i < b; ++i)
                    rateVector[i] = rate_calculator(i, beta, magmoms[i], randomNumberVector3D[i], magmoms);
            });

        generateAcceptanceVec(acceptanceVector);
        pool.wait_for_tasks();

        if (measurement)
            energyHistory[step - 1] = energy.get();


        pool.push_loop(atomnum,
            [&](const int a, const int b)
            {
                for (int i = a; i < b; ++i)
                    if (acceptanceVector[i] < rateVector[i])
                    {
                        magmoms[i] = randomNumberVector3D[i];
                    }
            });

        pool.wait_for_tasks();
        */

        if (approachTemp)
            runningTemperature *= approachValue;
        else if (approachMag)
        {
            runningMag += approachValue;
            zeemanTerm[magDir] = runningMag;
        }

        if(step == steps)
            for(auto & worker : mcWorkers)
                worker.request_stop();

        syncPointRun.arrive_and_wait();
        magmoms = magmomCopy[0];

    }
    for (auto& worker : mcWorkers)
        worker.join();
    if (measurement)
        tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in measurement mode!" << std::endl << std::endl;
    else
        tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in equilib mode!" << std::endl << std::endl;
}

void SimEnvironment::setTemperature(double temp)
{
    temperature = temp;
}

void SimEnvironment::setMagneticField(double xH, double yH, double zH)
{
    if (abs(xH) > 0.00001)
        magDir = 0;
    else if (abs(yH) > 0.00001)
        magDir = 1;
    else if (abs(zH) > 0.00001)
        magDir = 2;
    else
        magDir = -1;

    zeemanTerm[0] = xH;
    zeemanTerm[1] = yH;
    zeemanTerm[2] = zH;
}

void SimEnvironment::setSingleIonAnisotropy(double xC, double yC, double zC)
{
    if (abs(xC) > 0.00001)
        singleIonDir = 0;
    else if (abs(yC) > 0.00001)
        singleIonDir = 1;
    else if (abs(zC) > 0.00001)
        singleIonDir = 2;
    else
        singleIonDir = -1;

    singleIonAnisotropyTerm[0] = xC;
    singleIonAnisotropyTerm[1] = yC;
    singleIonAnisotropyTerm[2] = zC;
}

void SimEnvironment::setCompassAnisotropy(double xComp, double yComp, double zComp)
{
    if (abs(xComp) > 0.00001)
        compassDir = 0;
    else if (abs(yComp) > 0.00001)
        compassDir = 1;
    else if (abs(zComp) > 0.00001)
        compassDir = 2;
    else
        compassDir = -1;

    compassAnisotropyTerm[0] = xComp;
    compassAnisotropyTerm[1] = yComp;
    compassAnisotropyTerm[2] = zComp;
}

void SimEnvironment::setOutputPath(std::string out)
{
    outputPath = out;
}

void SimEnvironment::setStatusSteps(int steps)
{
    statusSteps = steps;
}

std::string SimEnvironment::getOutputPath()
{
    return outputPath;
}

double SimEnvironment::energy_diff_calculator(const int& index, const std::vector<double>& oldMom, const std::vector<double>& newMom, const std::vector<std::vector<double>>& magmoms, const std::vector<std::vector<std::tuple<int, double>>>& atomlinks)
{
    double Hold = 0.0;
    double Hnew = 0.0;
    int link;
    double param, sum;
    for (const auto& linkforatom : atomlinks[index])
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
        Hold -= oldMom[i] * oldMom[i] * singleIonAnisotropyTerm[i];
        Hnew -= newMom[i] * newMom[i] * singleIonAnisotropyTerm[i];

        // 0.672 -> unity conversion factor mu_B/k_B https://pubs.aip.org/aip/adv/article/5/12/127124/661186/Modeling-of-hysteresis-loops-by-Monte-Carlo

        Hold -= oldMom[i] * zeemanTerm[i] * 0.672;
        Hnew -= newMom[i] * zeemanTerm[i] * 0.672;
    }
    // Compass anisotropy
    if(enableCompassAnisotropy)
        for(int i = 0; i < 3; i++)
            for (auto linkforatom : linksNN[i][index])
            {
                Hold -= oldMom[i] * magmoms[linkforatom][i] * compassAnisotropyTerm[i];
                Hnew -= newMom[i] * magmoms[linkforatom][i] * compassAnisotropyTerm[i];
            }

    return Hnew - Hold;
}

double SimEnvironment::rate_calculator(int& index, double& beta, std::vector<double>& oldMom, std::vector<double>& newMom, const std::vector<std::vector<double>>& magmoms,const std::vector<std::vector<std::tuple<int, double>>>& atomlinks)
{
    return std::exp(-energy_diff_calculator(index, oldMom, newMom, magmoms, atomlinks) * beta);
}

double SimEnvironment::energy_calculator()
{
    double energy = 0.0;
    static BS::thread_pool calcpool;
    static std::mutex energyLock;

    calcpool.push_loop(atomnum,
        [&](const int a, const int b)
        {
            double tempEnergy = 0.0;
            int link;
            double param, sum;
            for (int i = a; i < b; i++)
            {
                for (auto linkforatom : links[i])
                {
                    link = std::get<0>(linkforatom);
                    param = std::get<1>(linkforatom);
                    dotProduct(sum, magmoms[i], magmoms[link]);
                    tempEnergy -= sum * param * 0.5;
                }
                // Single Ion anisotropy + Zeeman Term
                for (int j = 0; j < 3; j++)
                {
                    tempEnergy -= magmoms[i][j] * magmoms[i][j] * singleIonAnisotropyTerm[j];

                    // 0.672 -> unity conversion factor mu_B/k_B https://pubs.aip.org/aip/adv/article/5/12/127124/661186/Modeling-of-hysteresis-loops-by-Monte-Carlo
                    tempEnergy -= magmoms[i][j] * zeemanTerm[j] * 0.672; 
                }
            }
            energyLock.lock();
            energy += tempEnergy;
            energyLock.unlock();
        });

    calcpool.wait_for_tasks();

    return energy;
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

void SimEnvironment::generateRandomVecArray(std::vector<std::vector<double>>& vecIn, std::uniform_real_distribution<>& randomizer, std::mt19937& generator)
{
    double C1, C2, Csq;

    for (int i = 0; i < vecIn.size(); i++)
    {
        Csq = 2;
        while (Csq >= 1)
        {
            C1 = 1.0 - 2.0 * randomizer(generator);
            C2 = 1.0 - 2.0 * randomizer(generator);
            Csq = C1 * C1 + C2 * C2;
        }
        vecIn[i][0] = 2 * C1 * std::sqrt((1 - Csq));
        vecIn[i][1] = 2 * C2 * std::sqrt((1 - Csq));
        vecIn[i][2] = 1 - 2 * Csq;
    }
}

void SimEnvironment::generateAcceptanceVec(std::vector<double>& vecIn)
{
    std::uniform_real_distribution<> ran(0.0, 1.0);
    for (int i = 0; i < vecIn.size(); i++)
    {
        vecIn[i] = ran(gen);
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

    auto magtotal = zeemanTerm[0] + zeemanTerm[1] + zeemanTerm[2];

    std::vector<double> returnVals = { temperature, MagMom, Chi, U4, E, HeatCapacity, singleIonAnisotropyTerm[2], magtotal };
    return returnVals;
}

void SimEnvironment::writeMagneticMomentsToFile(std::string path)
{
    // Open the output file
    std::ostringstream KString;
    KString << std::scientific << std::setprecision(1) << singleIonAnisotropyTerm[0] + singleIonAnisotropyTerm[1] + singleIonAnisotropyTerm[2];
    std::string interactionKString = KString.str();

    std::ostringstream CString;
    CString << std::scientific << std::setprecision(1) << compassAnisotropyTerm[0];
    std::string interactionCString = CString.str();

    std::ostringstream HString;
    HString << std::scientific << std::setprecision(1) << zeemanTerm[0] + zeemanTerm[1] + zeemanTerm[2];
    std::string magneticFieldHString = HString.str();


    if (!fs::exists(path)) {
        // The directory does not exist, create it
        fs::create_directories(path);
    }

    std::ofstream outFile(joinPaths(path, "endMagmoms_temp" + std::to_string(temperature) + "_K" + interactionKString + "_C" + interactionCString + "_H" + magneticFieldHString + ".magres"));
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
