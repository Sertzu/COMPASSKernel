#include "timeUtility.h"

DateTimeInfo getCurrentTime() {
    // Get the current time using the chrono library
    auto now = std::chrono::system_clock::now();
    // Convert the time to a time_t object representing Unix time
    std::time_t unix_time = std::chrono::system_clock::to_time_t(now);
    // Convert the time_t object to a tm struct representing a date object
    std::tm date_object;

#ifdef _WIN32
    localtime_s(&date_object, &unix_time);
#else
    localtime_r(&date_object, &unix_time);
#endif

    // Create a DateTimeInfo object and set the date_object and unix_time fields
    DateTimeInfo dt_info;
    dt_info.date_object = date_object;
    dt_info.unix_time = unix_time;
    char date_buffer[100];
    std::strftime(date_buffer, sizeof(date_buffer), "%Y-%m-%d %H:%M:%S", &date_object);
    std::string formatted_date(date_buffer);
    dt_info.time_string = formatted_date;
    return dt_info;
}