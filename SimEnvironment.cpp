#include "SimEnvironment.h"
#include <filesystem> // Requires C++17
#include <charconv>

namespace fs = std::filesystem;

inline void dotProduct(float& sum, const std::vector<float>& vec1, const std::vector<float>& vec2) {
    sum = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

std::vector<std::string> read_file_without_comments(const std::string& file_path) {
    std::ifstream input_file(file_path);
    std::vector<std::string> lines;

    if (!input_file.is_open()) {
        std::cerr << "Unable to open file: " << file_path << std::endl;
        return lines; // Return an empty vector if file can't be opened.
    }

    // Count the total number of lines in a single pass.
    size_t num_lines = std::count(std::istreambuf_iterator<char>(input_file),
        std::istreambuf_iterator<char>(), '\n');

    // Reserve space in the vector based on the total number of lines.
    lines.reserve(num_lines);

    // Seek back to the beginning of the file to read it.
    input_file.clear(); // Clear eof and any other flags.
    input_file.seekg(0, std::ios::beg);

    // Read the file and store the lines, skipping those not needed.
    std::string line;
    while (std::getline(input_file, line)) {
        if (!(line.empty() || line[0] == '#' || line[0] == '!')) {
            lines.push_back(line);
        }
    }

    return lines;
}

void extractLinksAndCoordinates(const std::string& line, std::tuple<int, int, float, bool, bool, bool, int, std::string, std::tuple<int, std::vector<float>>>& values_in) {
    std::istringstream iss(line);

    float values[30]; 
    std::fill_n(values, 30, std::numeric_limits<float>::quiet_NaN());
    int value_count = 0; 

    std::string atomName;
    std::string token;

    while (iss >> token) {
        // Check if the token can be converted to a float and if we have not exceeded the array size.
        if (!token.empty() && value_count < 30 && std::all_of(token.begin(), token.end(), [](char c) { return std::isdigit(c) || c == '.' || c == '-'; })) {
            float value;
            auto result = std::from_chars(token.data(), token.data() + token.size(), value);
            if (result.ec == std::errc()) {
                values[value_count++] = value;
            }
            else {
                atomName = token;
                break;
            }
        }
        else {
            atomName = token;
            break;
        }
    }

    std::get<0>(values_in) = static_cast<int>(values[0]);
    std::get<1>(values_in) = static_cast<int>(values[2]);
    std::get<2>(values_in) = values[8];

    bool isNearestNeighbourX = false;
    bool isNearestNeighbourY = false;
    bool isNearestNeighbourZ = false;
    
    if(value_count == 19)
        if(((abs(static_cast<int>(values[12])) + abs(static_cast<int>(values[13])) + abs(static_cast<int>(values[14]))) < 2) 
            && static_cast<int>(values[1]) == static_cast<int>(values[3]))
        {
            if (abs(static_cast<int>(values[12])) == 1)
                isNearestNeighbourX = true;
            if (abs(static_cast<int>(values[13])) == 1)
                isNearestNeighbourY = true;
            if (abs(static_cast<int>(values[14])) == 1)
                isNearestNeighbourZ = true;
        }

    std::get<3>(values_in) = isNearestNeighbourX;
    std::get<4>(values_in) = isNearestNeighbourY;
    std::get<5>(values_in) = isNearestNeighbourZ;

    std::get<6>(values_in) = static_cast<int>(values[18]);
    std::get<7>(values_in) = atomName;

    std::get<0>(std::get<8>(values_in)) = static_cast<int>(values[0]);

    std::get<1>(std::get<8>(values_in)).resize(3);
    std::get<1>(std::get<8>(values_in))[0] = values[9];
    std::get<1>(std::get<8>(values_in))[1] = values[10];
    std::get<1>(std::get<8>(values_in))[2] = values[11];
}


std::tuple<int, std::vector<float>> extractCoordinates(const std::string& line)
{
    std::istringstream iss(line);
    int first, third, counter;
    float last;
    static std::vector<float> values(20, 0.0);
    float value;
    static std::vector<float> coords{0.0,0.0,0.0};
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
    m_inputPath = "";
    m_outputPath = "";
    m_temperature = 0;
}

SimEnvironment::SimEnvironment(std::string input, std::string output, float temp) : gen(rd())
{
    m_threadnum = 0;
    m_tlog = ThreadLogger(m_threadnum, true);
    m_inputPath = input;
    m_outputPath = output;
    m_temperature = temp;

    if (m_inputPath == "")
        return;
    m_tlog << getCurrentTime().time_string <<  " [SimEnvironment] Setting up SimEnvironment!" << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Reading input file: " << input << std::endl;
    auto lines = read_file_without_comments(m_inputPath);
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Read input file, converting data!" << std::endl;

    std::tuple<int, int, float, bool, bool, bool, int, std::string, std::tuple<int, std::vector<float>>> values;
    extractLinksAndCoordinates(lines.back(), values);

    m_atomnum = std::get<0>(values);
    m_linksNN.resize(3);
    m_links.resize(m_atomnum);
    m_linksNN[0].resize(m_atomnum);
    m_linksNN[1].resize(m_atomnum);
    m_linksNN[2].resize(m_atomnum);
    m_atomCoordinates.resize(m_atomnum);
    m_atomTypes.resize(m_atomnum);
    m_atomNames.resize(m_atomnum);

    auto coords = std::get<8>(values);

    for (int i = 0; i < lines.size(); i++)
    {
        extractLinksAndCoordinates(lines[i], values);
        coords = std::get<8>(values);
        m_links[std::get<0>(values) - 1].push_back(std::make_tuple(std::get<1>(values) - 1, std::get<2>(values)));
        if(std::get<3>(values))
        {
            m_linksNN[0][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<4>(values))
        {
            m_linksNN[1][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }
        if (std::get<5>(values))
        {
            m_linksNN[2][std::get<0>(values) - 1].push_back(std::get<1>(values) - 1);
        }

        m_atomCoordinates[std::get<0>(coords) - 1] = coords;
        m_atomTypes[std::get<0>(values) - 1] = std::get<6>(values);
        m_atomNames[std::get<0>(values) - 1] = std::get<7>(values);

    }
    m_singleIonAnisotropyTerm = { 0.0, 0.0, 0.0 };
    m_compassAnisotropyTerm = { 0.0, 0.0, 0.0 };
    m_zeemanTerm = { 0.0, 0.0, 0.0 };
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Initializing magnetic moments!" << std::endl;
    m_magmoms.reserve(m_atomnum);
    for (int i = 0; i < m_atomnum; i++)
        m_magmoms.push_back(generateRandomVecSingle());
    if ((abs(m_compassAnisotropyTerm[0]) + abs(m_compassAnisotropyTerm[1]) + abs(m_compassAnisotropyTerm[2])) > 1e-8)
        m_enableCompassAnisotropy = true;
    else
        m_enableCompassAnisotropy = false;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Setup complete!" << std::endl << std::endl;
}

std::vector<std::vector<std::tuple<int, float>>> SimEnvironment::getLinks()
{
    return m_links;
}

std::vector<std::vector<std::vector<int>>> SimEnvironment::getLinksNN()
{
    return m_linksNN;
}

std::vector<std::tuple<int, std::vector<float>>> SimEnvironment::getAtomCoordinates()
{
    return m_atomCoordinates;
}

std::vector<std::vector<float>> SimEnvironment::getMagmoms()
{
    return m_magmoms;
}

void SimEnvironment::runSim(int steps, bool measurement, bool approachTemp, float initialtemp, bool approachMag, float initialH)
{
    m_magmomsHistory.clear();
    m_meanmagmomsHistory.clear();
    m_energyHistory.clear();
    if (measurement)
        m_tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in measurement mode for " << steps << " steps!" << std::endl;
    else if(approachTemp)
        m_tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in approach mode for " << steps << " steps!" << std::endl;
    else
        m_tlog << getCurrentTime().time_string << " [Kernel] Starting MC simulation in equilib mode for " << steps << " steps!" << std::endl;
    auto lasttime = std::chrono::high_resolution_clock::now();
    float approachValue;
    float runningTemperature = m_temperature;
    float runningMag = m_zeemanTerm[m_magDir];
    if (approachTemp)
    {
        runningTemperature = initialtemp;
        approachValue = std::pow(m_temperature / initialtemp, 1.0 / steps);
        m_tlog << getCurrentTime().time_string << " [Kernel] Initial temperature is " << initialtemp << "K and Approachfactor is " << approachValue << std::endl;
    }
    else if (approachMag)
    {
        runningMag = initialH;
        approachValue = (m_zeemanTerm[m_magDir] - initialH) / steps;
        m_zeemanTerm[m_magDir] = runningMag;
        m_tlog << getCurrentTime().time_string << " [Kernel] Initial Magnetic Field is " << initialH << "T and Approachvalue is " << approachValue << std::endl;
    }

    std::vector<std::vector<float>> randomNumberVector3D(m_atomnum, { 0.0 ,0.0 , 0.0});
    std::vector<float> rateVector(m_atomnum, 0.0);
    std::vector<float> acceptanceVector(m_atomnum, 0.0);

    //generateRandomVecArray(randomNumberVector3D);
    //generateAcceptanceVec(acceptanceVector);

    std::vector<std::vector<float>> randomNumberVector3DSwap(m_atomnum, { 0.0 ,0.0 , 0.0 });
    std::vector<float> acceptanceVectorSwap(m_atomnum, 0.0);

    BS::thread_pool pool;
    m_magmomsHistory.clear();
    m_meanmagmomsHistory.clear();
    m_meanRawMagmomsHistory.clear();

    m_energyHistory.clear();

    m_meanmagmomsHistory.resize(steps);
    m_meanRawMagmomsHistory.resize(steps);
    m_energyHistory.resize(steps);
    std::future<float> energy;

    std::vector<std::vector<std::vector<float>>> magmomCopy;
    for (int i = 0; i < m_workerCount; i++)
    {
        magmomCopy.push_back(m_magmoms);
    }
    auto startIndices = getStartIndices(m_magmoms.size(), m_workerCount);

    std::barrier syncPointInit(m_workerCount + 1);
    std::barrier syncPointRun(m_workerCount + 1);

    auto mcWorker = [&](std::stop_token stopToken, float& beta, std::vector<std::vector<float>> &magmoms, std::vector<std::vector<float>>& randomNumberVector3D, std::vector<float> &acceptanceVector, int start, int end)
    {
        std::uniform_real_distribution<> randomizer(0.0, 1.0);
        std::random_device randDevice;
        std::mt19937 generator(randDevice());
        auto atomlinks = m_links;
        float rate = 0.0;
        float acceptor = 0.0;
        std::vector<std::vector<float>> randomVec = { {0.0,0.0,0.0} };
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
    float beta = 1 / (m_kB * runningTemperature);

    mcWorkers.reserve(m_workerCount);
    for (int i = 0; i < m_workerCount; ++i) {
        int start = startIndices[i];
        int end = (i == m_workerCount - 1) ? m_magmoms.size() : startIndices[i + 1];
        mcWorkers.emplace_back(mcWorker, std::ref(beta), std::ref(magmomCopy[i]), std::ref(randomNumberVector3D),std::ref(acceptanceVector), start, end);
    }

    for (int step = 1; step <= steps; step++)
    {
        if (step % m_statusSteps == 0)
        {
            auto currenttime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currenttime - lasttime);
            if (measurement)
                m_tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " of " << steps << " in measurement mode. Last " << m_statusSteps << " steps took: "
                << duration.count() << " ms!" << std::endl;
            else if (approachTemp)
                m_tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in approachTemp mode. Current temperature is " 
                << runningTemperature << "K !" << " Last " << m_statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            else if (approachMag)
                m_tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in approachMag mode. Current Magnetic Field is "
                << m_zeemanTerm[m_magDir] << "T !" << " Last " << m_statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            else
                m_tlog << getCurrentTime().time_string << " [Kernel] Currently at step " << step << " / " << steps << " in equilib mode. Last " << m_statusSteps << " steps took : "
                << duration.count() << " ms!" << std::endl;
            lasttime = currenttime;
        }
        if (measurement)
        {
            
            float M1 = 0.0;
            float M2 = 0.0;
            float M3 = 0.0;

            float M1Stag = 0.0;
            float M2Stag = 0.0;
            float M3Stag = 0.0;

            if(m_useStaggeredMagnetization)
            {
                for (int i = 0; i < m_atomnum; i++)
                {
                    M1 += m_magmoms[i][0];
                    M2 += m_magmoms[i][1];
                    M3 += m_magmoms[i][2];

                    float pattern = m_magnetizationPattern[m_atomTypes[i] - 1];
                    M1Stag += m_magmoms[i][0] * pattern;
                    M2Stag += m_magmoms[i][1] * pattern;
                    M3Stag += m_magmoms[i][2] * pattern;

                }
            }
            else
            {
                for (int i = 0; i < m_atomnum; i++)
                {
                    M1 += m_magmoms[i][0];
                    M2 += m_magmoms[i][1];
                    M3 += m_magmoms[i][2];
                }

                M1Stag = M1;
                M2Stag = M2;
                M3Stag = M3;
            }

            M1 /= m_atomnum;
            M2 /= m_atomnum;
            M3 /= m_atomnum;

            M1Stag /= m_atomnum;
            M2Stag /= m_atomnum;
            M3Stag /= m_atomnum;

            m_meanmagmomsHistory[step - 1] = std::sqrt(M1Stag * M1Stag + M2Stag * M2Stag + M3Stag * M3Stag);

            m_meanRawMagmomsHistory[step - 1][0] = M1;
            m_meanRawMagmomsHistory[step - 1][1] = M2;
            m_meanRawMagmomsHistory[step - 1][2] = M3;

            energy = pool.submit(&SimEnvironment::energy_calculator, this);
            //energyHistory[step - 1] = energy_calculator();
            if (steps - step < m_MAXMAGMOMHISTSIZE)
                m_magmomsHistory.push_back(m_magmoms);
        }

        if (step == 1)
        {
            if (measurement)
                m_energyHistory[step - 1] = energy.get();
            syncPointInit.arrive_and_wait();
        }

        if (measurement && step != 1)
            m_energyHistory[step - 1] = energy.get();

        syncPointRun.arrive_and_wait();

        if (approachTemp)
            runningTemperature *= approachValue;
        else if (approachMag)
        {
            runningMag += approachValue;
            m_zeemanTerm[m_magDir] = runningMag;
        }

        if(step == steps)
            for(auto & worker : mcWorkers)
                worker.request_stop();

        syncPointRun.arrive_and_wait();
        m_magmoms = magmomCopy[0];

    }
    for (auto& worker : mcWorkers)
        worker.join();
    if (measurement)
        m_tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in measurement mode!" << std::endl << std::endl;
    else
        m_tlog << getCurrentTime().time_string << " [Kernel] Finished MC simulation in equilib mode!" << std::endl << std::endl;
}

void SimEnvironment::setTemperature(float temp)
{
    m_temperature = temp;
}

void SimEnvironment::setMagneticField(float xH, float yH, float zH)
{
    if (abs(xH) > 0.00000001)
        m_magDir = 0;
    else if (abs(yH) > 0.00000001)
        m_magDir = 1;
    else if (abs(zH) > 0.00000001)
        m_magDir = 2;
    else
        m_magDir = -1;

    m_zeemanTerm[0] = xH;
    m_zeemanTerm[1] = yH;
    m_zeemanTerm[2] = zH;
}

void SimEnvironment::setSingleIonAnisotropy(float xC, float yC, float zC)
{
    if (abs(xC) > 0.00001)
        m_singleIonDir = 0;
    else if (abs(yC) > 0.00001)
        m_singleIonDir = 1;
    else if (abs(zC) > 0.00001)
        m_singleIonDir = 2;
    else
        m_singleIonDir = -1;

    m_singleIonAnisotropyTerm[0] = xC;
    m_singleIonAnisotropyTerm[1] = yC;
    m_singleIonAnisotropyTerm[2] = zC;
}

void SimEnvironment::setCompassAnisotropy(float xComp, float yComp, float zComp)
{
    if (abs(xComp) > 0.00001)
        m_compassDir = 0;
    else if (abs(yComp) > 0.00001)
        m_compassDir = 1;
    else if (abs(zComp) > 0.00001)
        m_compassDir = 2;
    else
        m_compassDir = -1;

    m_compassAnisotropyTerm[0] = xComp;
    m_compassAnisotropyTerm[1] = yComp;
    m_compassAnisotropyTerm[2] = zComp;

    if ((abs(m_compassAnisotropyTerm[0]) + abs(m_compassAnisotropyTerm[1]) + abs(m_compassAnisotropyTerm[2])) > 1e-8)
        m_enableCompassAnisotropy = true;
    else
        m_enableCompassAnisotropy = false;
}

void SimEnvironment::setOutputPath(std::string out)
{
    m_outputPath = out;
}

void SimEnvironment::setMagPattern(std::vector<double> mag_pattern)
{
    std::vector<float> tempVec;
    std::transform(mag_pattern.begin(), mag_pattern.end(), std::back_inserter(tempVec),
        [](double d) { return static_cast<float>(d); });

    m_magnetizationPattern = tempVec;
    m_useStaggeredMagnetization = true;
}

void SimEnvironment::setStatusSteps(int steps)
{
    m_statusSteps = steps;
}

std::string SimEnvironment::getOutputPath()
{
    return m_outputPath;
}

float SimEnvironment::energy_diff_calculator(const int& index, const std::vector<float>& oldMom, const std::vector<float>& newMom, const std::vector<std::vector<float>>& magmoms, const std::vector<std::vector<std::tuple<int, float>>>& atomlinks)
{
    float Hold = 0.0;
    float Hnew = 0.0;
    int link;
    float param, sum;
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
        Hold -= oldMom[i] * oldMom[i] * m_singleIonAnisotropyTerm[i];
        Hnew -= newMom[i] * newMom[i] * m_singleIonAnisotropyTerm[i];

        // 0.672 -> unity conversion factor mu_B/k_B https://pubs.aip.org/aip/adv/article/5/12/127124/661186/Modeling-of-hysteresis-loops-by-Monte-Carlo

        Hold -= oldMom[i] * m_zeemanTerm[i];
        Hnew -= newMom[i] * m_zeemanTerm[i];
    }
    // Compass anisotropy
    if(m_enableCompassAnisotropy)
        for(int i = 0; i < 3; i++)
            for (auto linkforatom : m_linksNN[i][index])
            {
                Hold -= oldMom[i] * magmoms[linkforatom][i] * m_compassAnisotropyTerm[i];
                Hnew -= newMom[i] * magmoms[linkforatom][i] * m_compassAnisotropyTerm[i];
            }

    return Hnew - Hold;
}

float SimEnvironment::rate_calculator(int& index, float& beta, std::vector<float>& oldMom, std::vector<float>& newMom, const std::vector<std::vector<float>>& magmoms,const std::vector<std::vector<std::tuple<int, float>>>& atomlinks)
{
    return std::exp(-energy_diff_calculator(index, oldMom, newMom, magmoms, atomlinks) * beta);
}

float SimEnvironment::energy_calculator()
{
    float energy = 0.0;
    static BS::thread_pool calcpool;
    static std::mutex energyLock;

    calcpool.push_loop(m_atomnum,
        [&](const int a, const int b)
        {
            float tempEnergy = 0.0;
            int link;
            float param, sum;
            for (int i = a; i < b; i++)
            {
                for (auto linkforatom : m_links[i])
                {
                    link = std::get<0>(linkforatom);
                    param = std::get<1>(linkforatom);
                    dotProduct(sum, m_magmoms[i], m_magmoms[link]);
                    tempEnergy -= sum * param * 0.5;
                }
                // Single Ion anisotropy + Zeeman Term
                for (int j = 0; j < 3; j++)
                {
                    tempEnergy -= m_magmoms[i][j] * m_magmoms[i][j] * m_singleIonAnisotropyTerm[j];

                    // 0.672 -> unity conversion factor mu_B/k_B https://pubs.aip.org/aip/adv/article/5/12/127124/661186/Modeling-of-hysteresis-loops-by-Monte-Carlo
                    tempEnergy -= m_magmoms[i][j] * m_zeemanTerm[j]; 
                }
                // Compass anisotropy
                if (m_enableCompassAnisotropy)
                    for (int j = 0; j < 3; j++)
                        for (auto linkforatom : m_linksNN[j][i])
                        {
                            tempEnergy -= m_magmoms[i][j] * m_magmoms[linkforatom][j] * m_compassAnisotropyTerm[j];
                        }
            }
            energyLock.lock();
            energy += tempEnergy;
            energyLock.unlock();
        });

    calcpool.wait_for_tasks();

    return energy;
}

std::vector<float> SimEnvironment::generateRandomVecSingle()
{
    float C1, C2, Csq;

    std::uniform_real_distribution<> ran(0.0, 1.0);

    Csq = 2;
    while (Csq >= 1)
    {
        C1 = 1.0 - 2.0 * ran(gen);
        C2 = 1.0 - 2.0 * ran(gen);
        Csq = C1 * C1 + C2 * C2;
    }
    float A, B, C;
    A = 2 * C1 * std::sqrt((1 - Csq));
    B = 2 * C2 * std::sqrt((1 - Csq));
    C = 1 - 2 * Csq;
    std::vector<float> ranvec = { A, B, C};
    return ranvec;

}

void SimEnvironment::generateRandomVecArray(std::vector<std::vector<float>>& vecIn, std::uniform_real_distribution<>& randomizer, std::mt19937& generator)
{
    float C1, C2, Csq;

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

void SimEnvironment::generateAcceptanceVec(std::vector<float>& vecIn)
{
    std::uniform_real_distribution<> ran(0.0, 1.0);
    for (int i = 0; i < vecIn.size(); i++)
    {
        vecIn[i] = ran(gen);
    }
}

std::vector<double> SimEnvironment::getParameters()
{
    float magmomFirstorder = 0.0;
    float magmomSecondorder = 0.0;
    float magmomFourthorder = 0.0;

    float rawMagmomXDir = 0.0;
    float rawMagmomYDir = 0.0;
    float rawMagmomZDir = 0.0;

    for (float value : m_meanmagmomsHistory)
    {
        magmomFirstorder += value;
        magmomSecondorder += value * value;
        magmomFourthorder += value * value * value * value;
    }

    for (const auto& magmomVal : m_meanRawMagmomsHistory)
    {
        rawMagmomXDir += magmomVal[0];
        rawMagmomYDir += magmomVal[1];
        rawMagmomZDir += magmomVal[2];
    }

    int historySize = m_meanmagmomsHistory.size();

    rawMagmomXDir /= historySize;
    rawMagmomYDir /= historySize;
    rawMagmomZDir /= historySize;

    magmomFirstorder /= historySize;
    magmomSecondorder /= historySize;
    magmomFourthorder /= historySize;

    float energyFirstorder = 0.0;
    float energySecondorder = 0.0;
    for (float value : m_energyHistory)
    {
        energyFirstorder += value / m_atomnum;
        energySecondorder += value * value / (m_atomnum * m_atomnum);
    }
    energyFirstorder /= m_energyHistory.size();
    energySecondorder /= m_energyHistory.size();    

    float MagMom = magmomFirstorder;
    float Chi = m_atomnum / (m_kB * m_temperature) * (magmomSecondorder - magmomFirstorder * magmomFirstorder);
    float U4 = 1 - magmomFourthorder / (3 * magmomSecondorder * magmomSecondorder);
    float E = energyFirstorder;
    float HeatCapacity = m_atomnum * m_atomnum / (m_kB * m_temperature * m_temperature) * (energySecondorder - energyFirstorder * energyFirstorder);

    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Parameters are:" << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Temperature: " << m_temperature << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] MagMom: " << MagMom << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] RawMagMomXDir: " << rawMagmomXDir << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] RawMagMomYDir: " << rawMagmomYDir << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] RawMagMomZDir: " << rawMagmomZDir << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Chi: " << Chi << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Energy: " << E << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] HeatCapacity: " << HeatCapacity << std::endl;
    m_tlog << getCurrentTime().time_string << " [SimEnvironment] Kumulante of the fourth order: " << U4 << std::endl << std::endl;

    auto magtotal = m_zeemanTerm[0] + m_zeemanTerm[1] + m_zeemanTerm[2];

    std::vector<double> returnVals = { m_temperature, MagMom, Chi, U4, E, HeatCapacity, m_singleIonAnisotropyTerm[2], magtotal, rawMagmomXDir, rawMagmomYDir, rawMagmomZDir };
    return returnVals;
}

void SimEnvironment::writeMagneticMomentsToFile(std::string path)
{
    // Open the output file
    std::ostringstream KString;
    KString << std::scientific << std::setprecision(1) << m_singleIonAnisotropyTerm[0] + m_singleIonAnisotropyTerm[1] + m_singleIonAnisotropyTerm[2];
    std::string interactionKString = KString.str();

    std::ostringstream CString;
    CString << std::scientific << std::setprecision(1) << m_compassAnisotropyTerm[0];
    std::string interactionCString = CString.str();

    std::ostringstream HString;
    HString << std::scientific << std::setprecision(1) << m_zeemanTerm[0] + m_zeemanTerm[1] + m_zeemanTerm[2];
    std::string magneticFieldHString = HString.str();


    if (!fs::exists(path)) {
        // The directory does not exist, create it
        fs::create_directories(path);
    }

    std::ofstream outFile(joinPaths(path, "endMagmoms_temp" + std::to_string(m_temperature) + "_K" + interactionKString + "_C" + interactionCString + "_H" + magneticFieldHString + ".magres"));
    if (!outFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    // Iterate through both vectors and write the entries to the file
    for (size_t i = 0; i < m_atomCoordinates.size(); ++i) {
        for (const auto& value : std::get<1>(m_atomCoordinates[i])) {
            outFile << formatDouble(value, 9, 5) << ' ';
        }
        float sumx, sumy, sumz;
        float sumxnorm, sumynorm, sumznorm;
        sumx = sumy = sumz = 0.;
        for (int j = 0; j < m_magmomsHistory.size(); j++) {
            sumx += m_magmomsHistory[j][i][0];
            sumy += m_magmomsHistory[j][i][1];
            sumz += m_magmomsHistory[j][i][2];
        }
        float norm = std::sqrt(sumx * sumx + sumy * sumy + sumz * sumz);

        outFile << formatDouble(sumx / norm, 9, 5) << ' ';
        outFile << formatDouble(sumy / norm, 9, 5) << ' ';
        outFile << formatDouble(sumz / norm, 9, 5) << ' ';

        outFile << m_atomTypes[i] << ' ';
        outFile << m_atomNames[i];
        outFile << '\n';
    }

    // Close the output file
    outFile.close();
}
