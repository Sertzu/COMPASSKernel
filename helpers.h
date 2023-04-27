#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>

double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2);
int read_args(int argc, char* argv[], std::string &input, std::string& output, double &temperature, int& equilibsteps, int& measuresteps);
void appendToFile(const std::string& filename, const std::string& data);
std::string formatDouble(double value, int width, int precision);
std::vector<double> generate_evenly_spaced_numbers(double a, double b, int n);