#include "scriptHandler.h"

bool canBeConvertedToInt(const std::string& s) {
    std::istringstream stream(s);
    int value;
    stream >> value;
    return stream.eof() && !stream.fail();
}

bool canBeConvertedToDouble(const std::string& s) {
    std::istringstream stream(s);
    double value;
    stream >> value;
    return stream.eof() && !stream.fail();
}

ScriptOption stringToOptionEnum(const std::string& str) {
    static const std::unordered_map<std::string, ScriptOption> strToEnumMap{
        {"LOAD", ScriptOption::LOAD},
        {"SETOUT", ScriptOption::SETOUT},
        {"LOADCHECKPOINT", ScriptOption::LOADCHECKPOINT},
        {"SETTEMP", ScriptOption::SETTEMP},
        {"SETMAGFIELD", ScriptOption::SETMAGFIELD},
        {"APPROACHTEMP", ScriptOption::APPROACHTEMP},
        {"APPROACHMAG", ScriptOption::APPROACHMAG},
        {"EQUILIB", ScriptOption::EQUILIB},
        {"MEASUREMENT", ScriptOption::MEASUREMENT},
        {"SWEEPTEMP", ScriptOption::SWEEPTEMP},
        {"SWEEPMAG", ScriptOption::SWEEPMAG}
    };

    auto it = strToEnumMap.find(str);
    if (it != strToEnumMap.end()) {
        return it->second;
    }
    else {
        return ScriptOption::INVALID;
    }
}

std::vector<std::vector<std::string>> readScriptFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file : " + filename);
    }

    std::vector<std::vector<std::string>> result;
    std::string line;
    while (getline(file, line)) {
        // Skip lines that are empty or start with '#' or '!'
        if (line.empty() || line[0] == '#' || line[0] == '!') {
            continue;
        }

        // Split the line by whitespace
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }

        result.push_back(tokens);
    }

    file.close();
    return result;
}

ScriptHandler::ScriptHandler(const std::string settingsPath) : simRunner_(settingsPath)
{
    tlogRunner_ = ThreadLogger(0, true);
}

void ScriptHandler::loadScript(const std::string scriptPath)
{
    auto tempScript = readScriptFile(scriptPath);
    validcheck_ = true;
    std::string invalidScriptCommand;

    for (const auto& element : tempScript)
    {
        if (stringToOptionEnum(element[0]) == ScriptOption::INVALID)
        {
            validcheck_ = false;
            invalidScriptCommand = element[0];
            break;
        }
    }

    if (!validcheck_)
        print(">>>WARNING: INVALID SCRIPTCOMMAND: " + invalidScriptCommand + "<<<");
    else
    {
        script_ = tempScript;
        currentPos_ = 0;
        print("LOADED SCRIPT FROM " + scriptPath + "\n");

        validateScript();
    }

}

void ScriptHandler::showScript()
{   
    if (!validcheck_)
    {
        print("NOT A VALID SCRIPT!");
        return;
    }
    std::string out;
    print("THE SCRIPT (CURRENT POS = ->)");
    for (int i = 0; i < script_.size(); i++)
    {
        out = std::to_string(i+1) + ": ";
        if (i == currentPos_)
            out += "-> ";
        for (const auto& element : script_[i])
        {
            out += element;
            out += " ";
        }
        if (i == script_.size() - 1)
            out += "\n";
        print(out);
    }
}

void ScriptHandler::validateScript()
{
    using namespace std::string_literals;

    if (!validcheck_)
    {
        print("COULDN'T VALIDATE DUE TO INVALID SCRIPT!");
        return;
    }

    for (const auto& command : script_)
    {
        auto commandType = stringToOptionEnum(command[0]);
        switch (commandType)
        {
            case(ScriptOption::LOAD):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND LOAD INVALID"); }
                break;
            case(ScriptOption::SETOUT):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SETOUT INVALID"); }
                break;
            case(ScriptOption::LOADCHECKPOINT):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND LOADCHECKPOINT INVALID"); }
                break;
            case(ScriptOption::SETTEMP):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SETTEMP INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETTEMP: TEMPERATURE IS NOT VALID "); }
                break;
            case(ScriptOption::SETMAGFIELD):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND SETMAGFIELD INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETMAGFIELD: FIELDSTRENGTH IS NOT VALID "); }
                if (!(command[2] == "X")) { print("HERE"); }
                if (command[2] != "X" && command[2] != "Y" && command[2] != "Z") { validcheck_ = false; print("COMMAND SETMAGFIELD: FIELDDIRECTION IS NOT VALID "); }
                break;
            case(ScriptOption::APPROACHTEMP):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND APPROACHTEMP INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND APPROACHTEMP: TEMPERATURE IS NOT VALID "); }
                if (!canBeConvertedToInt(command[2])) { validcheck_ = false; print("COMMAND APPROACHTEMP: STEPS IS NOT VALID "); }
                break;
            case(ScriptOption::APPROACHMAG):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND APPROACHMAG INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND APPROACHMAG: FIELD STRENGTH IS NOT VALID "); }
                if (!canBeConvertedToInt(command[2])) { validcheck_ = false; print("COMMAND APPROACHMAG: STEPS IS NOT VALID "); }
                break;
            case(ScriptOption::EQUILIB):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND EQUILIB INVALID"); }
                if (!canBeConvertedToInt(command[1])) { validcheck_ = false; print("COMMAND EQUILIB: STEPS IS NOT VALID "); }
                break;
            case(ScriptOption::MEASUREMENT):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND MEASUREMENT INVALID"); }
                if (!canBeConvertedToInt(command[1])) { validcheck_ = false; print("COMMAND MEASUREMENT: STEPS IS NOT VALID "); }
                break;
            case(ScriptOption::SWEEPTEMP):
                if (command.size() != 6) { validcheck_ = false; print("COMMAND SWEEPTEMP INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SWEEPTEMP: TEMPERATURE IS NOT VALID "); }
                if (!canBeConvertedToDouble(command[2])) { validcheck_ = false; print("COMMAND SWEEPTEMP: TEMPERATURESTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[3])) { validcheck_ = false; print("COMMAND SWEEPTEMP: APPROACHSTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[4])) { validcheck_ = false; print("COMMAND SWEEPTEMP: EQUILIBSTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[5])) { validcheck_ = false; print("COMMAND SWEEPTEMP: MEASUREMENTSTEPS IS NOT VALID "); }
                break;
            case(ScriptOption::SWEEPMAG):
                if (command.size() != 6) { validcheck_ = false; print("COMMAND SWEEPMAG INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SWEEPMAG: FIELDSTRENGTH IS NOT VALID "); }
                if (!canBeConvertedToDouble(command[2])) { validcheck_ = false; print("COMMAND SWEEPMAG: FIELDSTRENGTHSTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[3])) { validcheck_ = false; print("COMMAND SWEEPMAG: APPROACHSTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[4])) { validcheck_ = false; print("COMMAND SWEEPMAG: EQUILIBSTEPS IS NOT VALID "); }
                if (!canBeConvertedToInt(command[5])) { validcheck_ = false; print("COMMAND SWEEPMAG: MEASUREMENTSTEPS IS NOT VALID "); }
                break;
            default:
                validcheck_ = false;
        }
    }
    if (validcheck_)
        print("VALIDATED SCRIPT!");
    else
        print(">>>THERE WERE ERRORS DURING THE VALIDATION OF THE SCRIPT<<<");
}

void ScriptHandler::runScript()
{
    if (!validcheck_)
    {
        print("COULDN'T START DUE TO INVALID SCRIPT!");
        return;
    }
}

void ScriptHandler::print(const std::string toPrint)
{
    auto time = getCurrentTime().time_string;
    auto msg = toPrint;
    tlogRunner_ << (time + " [ScriptHandler] " + msg + "\n");
}