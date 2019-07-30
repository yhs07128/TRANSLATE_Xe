#include <iostream>
#include <chrono>
#include <thread>
#include <string>

#include "ProgressBar.h"

void ProgressBar::update(double prog, int thread)
{
    if (prog > progress[thread])
        progress[thread] = prog;
}

void ProgressBar::new_file(int file, int thread)
{
    file_number[thread] = file;
    progress[thread] = 0;
}

void ProgressBar::display() const
{
    std::cout << "\r";
    for (int k = 0; k < num_threads; k++)
        std::cout << "File " << file_number[k] << ": " << progress[k] * 100 << "%   ";
    std::cout << std::flush;
}

bool ProgressBar::min_prog(int thread) const
{
    for (int k = 0; k < num_threads; k++)
        if (progress[thread] > progress[k])
            return false;

    return true;
}