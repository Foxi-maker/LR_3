#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "Matrix.h"
#include "Function.h"

const std::vector<std::string> fileNames = { "Test1.txt","EIGEN10.txt" };
//const std::vector<std::string> fileNames = { "Test1.txt" };

const double eps = 1.e-2;

const double ZERO_DOUBLE = 1.e-15;