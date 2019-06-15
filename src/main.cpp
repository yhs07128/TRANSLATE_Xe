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

#include "Vec.h"
#include "Electron.h"

int main()
{
    std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> std_gauss(0.0, 1.0);
    std::normal_distribution<double> argon_dist(0.0, sqrt((k_b * T) / (2 * M_a)));
    std::normal_distribution<double> elec_dist(0.0, sqrt((k_b * T) / (2 * m_e)));

    int dimensions;
    std::cout << "1D or 3D: ";
    std::cin >> dimensions;
    while (std::cin.fail() || (dimensions != 1 && dimensions !=3))
    {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> dimensions;
    }

    std::string line;
    int volts;
    std::vector<int> volts_list;

    std::cout << "V/cm (multiple values separated by spaces): ";
    std::cin.ignore();
    std::getline(std::cin, line);
    std::istringstream stream(line);
    while (stream >> volts)
        volts_list.push_back(volts);

    double cutoff;
    std::cout << "Stop After (s): ";
    std::cin >> cutoff;
    while (std::cin.fail())
    {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> cutoff;
    }

    int write_every;
    std::cout << "Record Every n Measurements: ";
    std::cin >> write_every;
    while (std::cin.fail())
    {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> write_every;
    }

    int batches;
    std::cout << "Repeat: ";
    std::cin >> batches;
    while (std::cin.fail())
    {
        std::cin.clear();
        std::cin.ignore();
        std::cin >> batches;
    }

    int ions = 0;

    for (int j = 0; j < volts_list.size(); j++)
    {
        for (int i = 0; i < batches; i++)
        {
            std::cout << "Writing " << int(volts_list[j]) << "V file " << i + 1 << " " << std::flush;

            std::string name = "sim_data/";
            
            if (dimensions == 1)
                name += "1D - " + std::to_string(int(volts_list[j])) + "V - " + std::to_string(i + 1) + ".txt";
            else if (dimensions == 3)
                name += "3D - " + std::to_string(int(volts_list[j])) + "V - " + std::to_string(i + 1) + ".txt";

            std::ofstream file(name);

            if (dimensions == 1)
            {
                Electron<double> elec(ions, volts_list[j], generator, std_gauss, argon_dist, elec_dist);
                generate_plot<double>(elec, cutoff, file, write_every);
            }
            else if (dimensions == 3)
            {
                Electron<Vec> elec(ions, volts_list[j], generator, std_gauss, argon_dist, elec_dist);
                generate_plot<Vec>(elec, cutoff, file, write_every);
            }

            file.close();

            std::cout << "\n";
        }
    }

    std::cout << "Complete." << std::endl;

    return 0;
}