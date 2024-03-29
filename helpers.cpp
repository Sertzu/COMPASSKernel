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

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> linspaced;

    if (num == 0) {
        return linspaced;
    }

    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i) {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // ensure that end is included

    return linspaced;
}

std::vector<double> logspace(double start, double end, int num, double base)
{
    std::vector<double> logscaled;

    if (num == 0) {
        return logscaled;
    }

    if (num == 1) {
        logscaled.push_back(std::pow(base, start));
        return logscaled;
    }

    std::vector<double> linspace_vector = linspace(start, end, num);

    for (double num : linspace_vector) {
        logscaled.push_back(std::pow(base, num));
    }

    return logscaled;
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


void appendToFile(const std::string& filename, const std::string& data) {
    std::ofstream datFile;

    std::string fullFileName = filename + ".dat";
    // Check if the file exists
    bool datFileExists = std::filesystem::exists(fullFileName);

    if (datFileExists) {
        // Open the file in append mode
        datFile.open(fullFileName, std::ios::app);
    }
    else {
        // Create a new file and add the header
        datFile.open(fullFileName, std::ios::out);
        datFile << "# Temp         M          Chi        U4         E     HeatCapacity       K        H\n";
    }

    // Append the string
    datFile << data << '\n';

    // Close the file
    datFile.close();
}

void writeResultsToFile(const std::string& filename, const std::vector<double>& data) {
    std::ofstream datFile;

    std::string fullFileNameDat = filename + ".dat";
    // Check if the file exists
    bool datFileExists = std::filesystem::exists(fullFileNameDat);

    if (datFileExists) {
        // Open the file in append mode
        datFile.open(fullFileNameDat, std::ios::app);
    }
    else {
        // Create a new file and add the header
        datFile.open(fullFileNameDat, std::ios::out);
        datFile << "# Temp         M          Chi        U4         E     HeatCapacity       K        H          M_x        M_y        M_z\n";
    }
    std::string out(formatDouble(data[0], 9, 5) + "  " + formatDouble(data[1], 9, 5) + "  " + formatDouble(data[2], 9, 5) + "  " + formatDouble(data[3], 9, 5) + "  " + formatDouble(data[4], 9, 5) + "  " + formatDouble(data[5], 9, 5) + "  " + formatDouble(data[6], 9, 5) + "  " + formatDouble(data[7], 9, 5) + "  " + formatDouble(data[8], 9, 5) + "  " + formatDouble(data[9], 9, 5) + "  " + formatDouble(data[10], 9, 5));
    // Append the string
    datFile << out << '\n';

    // Close the file
    datFile.close();

    // Same for the csv file
    std::ofstream csvFile;
    std::string fullFileNameCsv = filename + ".csv";
    bool csvFileExists = std::filesystem::exists(fullFileNameCsv);

    if (csvFileExists) {
        // Open the file in append mode
        csvFile.open(fullFileNameCsv, std::ios::app);
    }
    else {
        // Create a new file and add the header
        csvFile.open(fullFileNameCsv, std::ios::out);
        csvFile << "temperature;magnetization;chi;kumulantU4;energy;heatCapacity;singleIonAnisotropyTerm;magneticFieldStrength;magnetizationXDir;magnetizationYDir;magnetizationZDir\n";
    }
    // Check if the file stream is open
    if (csvFile.is_open()) {
        // Iterate over the vector
        for (size_t i = 0; i < data.size(); ++i) 
        {
            csvFile << data[i];

            // Add a semicolon after each element except the last
            if (i < data.size() - 1) 
            {
                csvFile << ";";
            }
            else
            {
                csvFile << "\n";
            }
        }
    }
    else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
    // Close the file stream
    csvFile.close();
}

std::string formatDouble(double value, int width, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setw(width) << std::setprecision(precision) << value;
    return ss.str();
}


std::vector<double> generate_evenly_spaced_numbers(double a, double b, int n) {
    std::vector<double> result;

    if (n < 2) {
        std::cerr << "Error: n should be at least 2." << std::endl;
        return result;
    }

    double step = (b - a) / (n - 1);

    for (int i = 0; i < n; ++i) {
        double value = a + i * step;
        result.push_back(value);
    }

    return result;
}

bool compare_lines(const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
    double val_a, val_b;
    std::istringstream iss_a(a.second), iss_b(b.second);
    iss_a >> val_a;
    iss_b >> val_b;
    return val_a < val_b;
}

void sort_file_by_first_entry(const std::string& filename) {
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::vector<std::pair<int, std::string>> lines;
    std::map<int, std::string> ignored_lines;
    std::string line;
    int line_number = 0;

    while (std::getline(input_file, line)) {
        if (!line.empty() && line[0] != '#') {
            lines.push_back(std::make_pair(line_number, line));
        }
        else {
            ignored_lines[line_number] = line;
        }
        ++line_number;
    }

    input_file.close();
    std::sort(lines.begin(), lines.end(), compare_lines);

    std::vector<std::string> sorted_lines(line_number);
    for (const auto& ignored_line : ignored_lines) {
        sorted_lines[ignored_line.first] = ignored_line.second;
    }

    for (const auto& valid_line : lines) {
        auto pos = std::find_if(sorted_lines.begin(), sorted_lines.end(), [](const std::string& s) { return s.empty(); });
        *pos = valid_line.second;
    }

    std::ofstream output_file(filename);
    if (!output_file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (const std::string& line_to_write : sorted_lines) {
        output_file << line_to_write << std::endl;
    }
}

std::string joinPaths(const std::string& path1, const std::string& path2)
{
    std::string joinedPath = path1;

    // Ensure the first path ends with a separator
    if (joinedPath.back() != PATH_SEPARATOR) {
        joinedPath += PATH_SEPARATOR;
    }

    // If the second path starts with a separator, remove it
    if (path2.front() == PATH_SEPARATOR) {
        joinedPath += path2.substr(1);
    }
    else {
        joinedPath += path2;
    }

    return joinedPath;
}

std::vector<int> getStartIndices(int size, int threadcount)
{
    std::vector<int> startIndices;

    // Calculate the number of elements each thread will handle.
    int chunkSize = size / threadcount;
    int remainder = size % threadcount;

    int currentIndex = 0;
    for (int i = 0; i < threadcount; i++) {
        startIndices.push_back(currentIndex);

        // If there's a remainder, distribute it among the threads.
        if (remainder > 0) {
            currentIndex += chunkSize + 1;
            remainder--;
        }
        else {
            currentIndex += chunkSize;
        }
    }

    return startIndices;
}

std::vector<std::pair<int, int>> splitArray(int x, int y) {
    std::vector<std::pair<int, int>> parts;

    if (y <= 0) { // Ensure parts count is greater than zero.
        return parts;
    }

    int partSize = x / y;
    int remainder = x % y;

    int start = 0;
    int end = partSize - 1;

    for (int i = 0; i < y; i++) {
        if (i < remainder) {
            end += 1;
        }
        parts.push_back({ start, end });
        start = end + 1;
        end = start + partSize - 1;
    }

    return parts;
}
