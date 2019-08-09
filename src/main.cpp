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

#include "Electron.h"
#include "ProgressBar.h"

void generate_plot(int volts, double cutoff, int cores, int write_every, int k, int batches, ProgressBar& bar)
{
    for (int i = 0; i < batches; i++) {
        bar.new_file(cores * i + (k + 1), k);
        int total_ions = 0;
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        Electron elec(total_ions, volts, starting_pos(generator), random_velocity(generator, true), generator);

        double starting_x = 0;

        if (!uniform_field)
            starting_x = elec.position().x;

        std::ofstream file("../py/studies/simulation-runs/" + std::to_string(volts) + "V - " + std::to_string(cores * i + (k + 1)) + ".txt");
        assert(file.is_open());
        
        int simulation_step = 1;
        const int total_progress_steps = 100;
        std::array<bool, total_progress_steps> progress;
        progress.fill(false);
        int progress_count = 1;

        while (elec.elapsed_time() < cutoff) {
            elec.update();

            if (!uniform_field && hit_check(elec.position()))
                break;

            if (simulation_step++ % write_every != 0)
                continue;

            if (uniform_field)
                bar.update(elec.elapsed_time() / cutoff, k);
            else
                bar.update(1 - (elec.position().x - z_min) / ((z_max - z_min) * spawn_height_scale), k);

            if (bar.min_prog(k))
                bar.display();
            
            file << elec.elapsed_time() * 1e9 << "," << elec.position().x * 1e6 << "," << elec.position().y * 1e6 << "," << elec.position().z * 1e6 << "," << elec.ke() << "," << ((elec.position().x - starting_x) / elec.elapsed_time()) << "," << total_ions << "\n";
        }

        elec.update();

        bar.update(1, k);

        file << elec.elapsed_time() * 1e9 << "," << elec.position().x * 1e6 << "," << elec.position().y * 1e6 << "," << elec.position().z * 1e6 << "," << elec.ke() << "," << ((elec.position().x - starting_x) / elec.elapsed_time()) << "," << total_ions << "\n";

        file.close();
        assert(!file.is_open());
    }

    return;
}

int main()
{
    std::cout << "V/cm (multiple values separated by spaces): ";

    std::string line;
    double volts;
    std::vector<double> volts_list;
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
        std::thread branches[cores];
    
        ProgressBar bar(cores);

        for (int k = 0; k < cores; k++) {
            branches[k] = std::thread(generate_plot, int(volts_list[j]), cutoff, cores, write_every, k, batches, std::ref(bar));   
        }
    
        for (int k = 0; k < cores; k++) {
            branches[k].join();
        }

        std::cout << "\n";
    }

    std::cout << "Complete." << std::endl;

    return 0;
}