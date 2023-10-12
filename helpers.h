#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <map>


#ifdef _WIN32
#include <windows.h>
#define PATH_SEPARATOR '\\'
#else
#include <unistd.h>
#define PATH_SEPARATOR '/'
#endif

std::vector<double> linspace(double start, double end, int num);
std::vector<double> logspace(double start, double end, int num, double base = 10.0);
std::vector<std::pair<int, int>> splitArray(int x, int y);

void dotProduct(double &sum, const std::vector<double>& vec1, const std::vector<double>& vec2);
int read_args(int argc, char* argv[], std::string &input, std::string& output, double &temperature, int& equilibsteps, int& measuresteps);
void appendToFile(const std::string& filename, const std::string& data);
std::string formatDouble(double value, int width, int precision);
std::vector<double> generate_evenly_spaced_numbers(double a, double b, int n);
void sort_file_by_first_entry(const std::string& filename);
std::string joinPaths(const std::string& path1, const std::string& path2);

std::vector<int> getStartIndices(int size, int threadcount);
