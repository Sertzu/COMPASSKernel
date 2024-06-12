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

	std::map<std::string, std::string> get_settings();
private:
	std::map<std::string, std::string> readSettingsFile(const std::string& fileName);
	std::map<std::string, std::string> settings;
};
