#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>

class ThreadLogger {
public:
    ThreadLogger(bool writeToConsole = false) : writeToConsole(writeToConsole) {
        if (!writeToConsole) {
            std::ostringstream filename;
            filename << "SimLog_" << std::this_thread::get_id() << ".txt";
            file.open(filename.str());
        }
    }

    template <typename T>
    ThreadLogger& operator<<(const T& value) {
        if (writeToConsole) {
            std::cout << value;
        }
        else {
            file << value;
        }
        return *this;
    }

    // Handling manipulators like std::endl
    ThreadLogger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (writeToConsole) {
            std::cout << manip;
        }
        else {
            file << manip;
        }
        return *this;
    }

private:
    bool writeToConsole;
    std::ofstream file;
};