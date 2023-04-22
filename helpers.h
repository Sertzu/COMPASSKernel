#pragma once
#include <string>
#include <iostream>
#include <vector>

double dotProduct(const std::vector<double>& vec1, const std::vector<double>& vec2);
int read_args(int argc, char* argv[], std::string &input, std::string& output, double &temperature, int& equilibsteps, int& measuresteps);