#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <utility>

struct DateTimeInfo {
    std::tm date_object;
    std::time_t unix_time;
    std::string time_string;
};

DateTimeInfo getCurrentTime();