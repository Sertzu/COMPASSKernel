#include "SettingsClass.h"

// Explicit specialization of the get method for std::string type
template<>
std::string SettingsClass::get<std::string>(const std::string& key) {
    return settings.at(key);
}

// Explicit specialization of the get method for int type
template<>
int SettingsClass::get<int>(const std::string& key) {
    try {
        return std::stoi(settings.at(key));
    }
    catch (...) {
        throw std::runtime_error("Cannot convert setting to int: " + key);
    }
}

// Explicit specialization of the get method for double type
template<>
double SettingsClass::get<double>(const std::string& key) {
    try {
        return std::stod(settings.at(key));
    }
    catch (...) {
        throw std::runtime_error("Cannot convert setting to double: " + key);
    }
}

SettingsClass::SettingsClass(std::string settingsFile)
{
    settings = readSettingsFile(settingsFile);
}

std::map<std::string, std::string> SettingsClass::get_settings()
{
    return settings;
}

std::map<std::string, std::string> SettingsClass::readSettingsFile(const std::string& fileName)
{
    std::map<std::string, std::string> settings;
    std::ifstream inputFile(fileName);

    if (inputFile.is_open()) {
        std::string line;

        while (std::getline(inputFile, line)) {
            // Trim leading and trailing spaces and ignore empty lines
            line.erase(0, line.find_first_not_of(' ')); // trim leading spaces
            line.erase(line.find_last_not_of(' ') + 1); // trim trailing spaces

            // Ignore lines that start with a "#" or are empty
            if (!line.empty() && line[0] != '#') {
                std::istringstream iss(line);
                std::string key, value;

                // Separate the line into key and value based on " = "
                if (std::getline(std::getline(iss, key, ' '), value, ' ') && std::getline(iss, value)) {
                    settings[key] = value;
                }
            }
        }

        inputFile.close();
    }
    else {
        std::cerr << "Unable to open file: " << fileName << std::endl;
    }

    return settings;
}