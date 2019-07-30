#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
#include <iomanip>

class ProgressBar
{
public:
    int num_threads;
    double progress[4];
    int file_number[4];

    ProgressBar(int threads): num_threads(threads), progress{0, 0, 0, 0}, file_number{1, 2, 3, 4}
    {
        std::cout << std::fixed << std::setprecision(3);
    }

    void update(double prog, int thread);
    void new_file(int file, int thread);
    void display() const;
    bool min_prog(int thread) const;
};

#endif