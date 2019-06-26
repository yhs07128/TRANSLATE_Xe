#include <math.h>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream>
#include <array>
#include <assert.h>
#include <thread>
#include <sys/stat.h>
#include <string>
#include <chrono>
#include <thread>

#include "Vec.h"
#include "Electron.h"

int main()
{
    int dimensions;
    std::cout << "1D or 3D: ";
    std::cin >> dimensions;
    while (std::cin.fail() || (dimensions != 1 && dimensions !=3)) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> dimensions;
    }

    std::string line;
    double volts;
    std::vector<double> volts_list;

    std::cout << "V/cm (multiple values separated by spaces): ";
    std::cin.ignore();
    std::getline(std::cin, line);
    std::istringstream stream(line);
    while (stream >> volts)
        volts_list.push_back(volts);

    double cutoff;
    std::cout << "Stop After (s): ";
    std::cin >> cutoff;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> cutoff;
    }

    int write_every;
    std::cout << "Record Every n Measurements: ";
    std::cin >> write_every;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> write_every;
    }

    int cores = std::thread::hardware_concurrency();
    if (cores < 1) cores = 1;

    int batches;
    std::cout << "Batches (in groups of " + std::to_string(cores) + "): ";
    std::cin >> batches;
    while (std::cin.fail()) {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> batches;
    }

    for (int j = 0; j < volts_list.size(); j++) {
        for (int i = 0; i < batches; i++) {
            std::cout << "Writing " << int(volts_list[j]) << "V file group " << i + 1 << " - 0%" << std::flush;
            std::thread branches[cores];

            for (int k = 0; k < cores; k++) {
                if (dimensions == 1)
                    branches[k] = std::thread(generate_plot<double>, int(volts_list[j]), cutoff, 4 * i + (k + 1), write_every, i + 1, "1d");
                else
                    branches[k] = std::thread(generate_plot<Vec>, int(volts_list[j]), cutoff, 4 * i + (k + 1), write_every, i + 1, "3d");   
            }
            
            for (int k = 0; k < cores; k++)
                    branches[k].join();

            std::cout << "\n";
        }
    }

    std::cout << "Complete." << std::endl;

    return 0;
}