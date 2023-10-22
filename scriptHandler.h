#pragma once

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "simRunner.h"

enum class ScriptOption {
    LOAD,
    SETOUT,
    LOADCHECKPOINT,
    SAVECHECKPOINT,
    SETMAGMOMAUTOSAVE,
    SETTEMP,
    SETMAGFIELD,
    APPROACHTEMP,
    APPROACHMAG,
    EQUILIB,
    MEASUREMENT,
    SWEEPTEMP,
    SWEEPMAG,
    SETEASYPLANE,

    INVALID
};

ScriptOption stringToOptionEnum(const std::string& str);

std::vector<std::vector<std::string>> readScriptFile(const std::string& filename);

class ScriptHandler {
public:
    ScriptHandler(const std::string settingsPath);

    void loadScript(const std::string scriptPath);

    void showScript();

    void validateScript();

    void runScript();

private:
    SimRunner simRunner_;
    ThreadLogger tlogRunner_;
    int currentPos_ = 0;
    bool validcheck_ = false;

    void print(const std::string toPrint);

    std::vector<std::vector<std::string>> script_;
};