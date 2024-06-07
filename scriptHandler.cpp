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
        {"SAVECHECKPOINT", ScriptOption::SAVECHECKPOINT},
        {"SETMAGMOMAUTOSAVE", ScriptOption::SETMAGMOMAUTOSAVE},
        {"SETTEMP", ScriptOption::SETTEMP},
        {"SETMAGFIELD", ScriptOption::SETMAGFIELD},
        {"APPROACHTEMP", ScriptOption::APPROACHTEMP},
        {"APPROACHMAG", ScriptOption::APPROACHMAG},
        {"EQUILIB", ScriptOption::EQUILIB},
        {"MEASUREMENT", ScriptOption::MEASUREMENT},
        {"SWEEPTEMP", ScriptOption::SWEEPTEMP},
        {"SWEEPMAG", ScriptOption::SWEEPMAG},
        {"SETEASYPLANE", ScriptOption::SETEASYPLANE},
        {"SETCOMPASSANISO", ScriptOption::SETCOMPASSANISO},
        {"SETMAGPATTERN", ScriptOption::SETMAGPATTERN},
        {"SETWORKERCOUNT", ScriptOption::SETWORKERCOUNT}
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

ScriptHandler::ScriptHandler()
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
            case(ScriptOption::SAVECHECKPOINT):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SAVECHECKPOINT INVALID"); }
                break; 
            case(ScriptOption::SETMAGMOMAUTOSAVE):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SETMAGMOMAUTOSAVE INVALID"); }
                break;
            case(ScriptOption::SETTEMP):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SETTEMP INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETTEMP: TEMPERATURE IS NOT VALID "); }
                break;
            case(ScriptOption::SETMAGFIELD):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND SETMAGFIELD INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETMAGFIELD: FIELD STRENGTH IS NOT VALID "); }
                if (!(command[2] == "X")) { print("HERE"); }
                if (command[2] != "X" && command[2] != "Y" && command[2] != "Z") { validcheck_ = false; print("COMMAND SETMAGFIELD: FIELD DIRECTION IS NOT VALID "); }
                break;
            case(ScriptOption::SETEASYPLANE):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND SETEASYPLANE INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETEASYPLANE: INTERACTION STRENGTH IS NOT VALID "); }
                if (command[2] != "X" && command[2] != "Y" && command[2] != "Z") { validcheck_ = false; print("COMMAND SETEASYPLANE: INTERACTION DIRECTION IS NOT VALID "); }
                break;
            case(ScriptOption::SETCOMPASSANISO):
                if (command.size() != 3) { validcheck_ = false; print("COMMAND SETCOMPASSANISO INVALID"); }
                if (!canBeConvertedToDouble(command[1])) { validcheck_ = false; print("COMMAND SETCOMPASSANISO: INTERACTION STRENGTH IS NOT VALID "); }
                if (command[2] != "X" && command[2] != "Y" && command[2] != "Z") { validcheck_ = false; print("COMMAND SETCOMPASSANISO: INTERACTION DIRECTION IS NOT VALID "); }
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
            case(ScriptOption::SETWORKERCOUNT):
                if (command.size() != 2) { validcheck_ = false; print("COMMAND SETWORKERCOUNT INVALID"); }
                if (!canBeConvertedToInt(command[1])) { validcheck_ = false; print("COMMAND SETWORKERCOUNT: WORKER COUNT IS NOT VALID "); }
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
            case(ScriptOption::SETMAGPATTERN):
            {
                bool first = true;
                for (const auto& command_part : command)
                {
                    if (first)
                    {
                        first = false;
                        continue;
                    }
                    if (!canBeConvertedToDouble(command_part))
                    {
                        validcheck_ = false;
                        print("COMMAND SETMAGPATTERN: INVALID VALUE IN LIST");
                    }
                }
                break;
            }
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

    for (const auto& command : script_)
    {
        showScript();
        auto commandType = stringToOptionEnum(command[0]);
        switch (commandType)
        {
        case(ScriptOption::LOAD):
            simRunner_.LOAD(command[1]);
            break;
        case(ScriptOption::SETOUT):
            simRunner_.SETOUT(command[1]);
            break;
        case(ScriptOption::LOADCHECKPOINT):
            simRunner_.LOADCHECKPOINT(command[1]);
            break;
        case(ScriptOption::SETMAGMOMAUTOSAVE):
            simRunner_.SETMAGMOMAUTOSAVE(command[1]);
            break;
        case(ScriptOption::SAVECHECKPOINT):
            simRunner_.SAVECHECKPOINT(command[1]);
            break;
        case(ScriptOption::SETTEMP):
            simRunner_.SETTEMP(std::stod(command[1]));
            break;
        case(ScriptOption::SETMAGFIELD):
            if (command[2] == "X")
                simRunner_.SETMAGFIELD(std::stod(command[1]), XYZ::X);
            if (command[2] == "Y")
                simRunner_.SETMAGFIELD(std::stod(command[1]), XYZ::Y);
            if (command[2] == "Z")
                simRunner_.SETMAGFIELD(std::stod(command[1]), XYZ::Z);
            break;
        case(ScriptOption::SETEASYPLANE):
            if (command[2] == "X")
                simRunner_.SETEASYPLANE(std::stod(command[1]), XYZ::X);
            if (command[2] == "Y")
                simRunner_.SETEASYPLANE(std::stod(command[1]), XYZ::Y);
            if (command[2] == "Z")
                simRunner_.SETEASYPLANE(std::stod(command[1]), XYZ::Z);
            break;
        case(ScriptOption::SETCOMPASSANISO):
            if (command[2] == "X")
                simRunner_.SETCOMPASSANISO(std::stod(command[1]), XYZ::X);
            if (command[2] == "Y")
                simRunner_.SETCOMPASSANISO(std::stod(command[1]), XYZ::Y);
            if (command[2] == "Z")
                simRunner_.SETCOMPASSANISO(std::stod(command[1]), XYZ::Z);
            break;
        case(ScriptOption::APPROACHTEMP):
            simRunner_.APPROACHTEMP(std::stod(command[1]), std::stoi(command[2]));
            break;
        case(ScriptOption::APPROACHMAG):
            simRunner_.APPROACHMAG(std::stod(command[1]), std::stoi(command[2]));
            break;
        case(ScriptOption::EQUILIB):
            simRunner_.EQUILIB(std::stoi(command[1]));
            break;
        case(ScriptOption::MEASUREMENT):
            simRunner_.MEASUREMENT(std::stoi(command[1]));
            break;
        case(ScriptOption::SETWORKERCOUNT):
            simRunner_.SETWORKERCOUNT(std::stoi(command[1]));
            break;
        case(ScriptOption::SWEEPTEMP):
            simRunner_.SWEEPTEMP(std::stod(command[1]), std::stod(command[2]), std::stoi(command[3]), std::stoi(command[4]), std::stoi(command[5]));
            break;
        case(ScriptOption::SWEEPMAG):
            simRunner_.SWEEPMAG(std::stod(command[1]), std::stod(command[2]), std::stoi(command[3]), std::stoi(command[4]), std::stoi(command[5]));
            break;
        case(ScriptOption::SETMAGPATTERN):
        {
            std::vector<double> mag_pattern;
            for (int i = 1; i < command.size(); i++)
            {
                mag_pattern.push_back(std::stod(command[i]));
            }
            simRunner_.SETMAGPATTERN(mag_pattern);
            break;
        }
        default:
            throw std::runtime_error("CRITICAL FAILURE IN THE SCRIPTRUNNER");
        }
        currentPos_ += 1;

    }

}

void ScriptHandler::print(const std::string toPrint)
{
    auto time = getCurrentTime().time_string;
    auto msg = toPrint;
    tlogRunner_ << (time + " [ScriptHandler] " + msg + "\n");
}
