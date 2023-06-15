#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <random>
#include <cmath>
#include <map>
#include <algorithm>


class SettingsClass
{
public:
	SettingsClass(std::string settingsFile);

    template<typename T>
    T get(const std::string& key);

    // Specializations
    template<>
    std::string get<std::string>(const std::string& key) {
        return settings.at(key);
    }

    template<>
    int get<int>(const std::string& key) {
        try {
            return std::stoi(settings.at(key));
        }
        catch (...) {
            throw std::runtime_error("Cannot convert setting to int: " + key);
        }
    }

    template<>
    double get<double>(const std::string& key) {
        try {
            return std::stod(settings.at(key));
        }
        catch (...) {
            throw std::runtime_error("Cannot convert setting to double: " + key);
        }
    }

	std::map<std::string, std::string> get_settings();
private:
	std::map<std::string, std::string> readSettingsFile(const std::string& fileName);
	std::map<std::string, std::string> settings;
};