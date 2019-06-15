#ifndef ELECTRON_H
#define ELECTRON_H

#include <random>
#include <array>
#include <vector>
#include <assert.h>

#include "Constants.h"
#include "Utils.h"
#include "Vec.h"

template<typename T>
class Electron
{
public:
    T x, v, accel;

    double volts_per_cm;
    double time_to_collision;
    double total_time;

    double energy;
    int& ionized;
    std::vector<Electron<T>*> ionization_electrons;

    std::mt19937 &generator;
    std::normal_distribution<double> &std_gauss;
    std::normal_distribution<double> &argon_dist;
    std::normal_distribution<double> &elec_dist;

    Electron(int &ions, double volts, std::mt19937 &gen, std::normal_distribution<double> &std, std::normal_distribution<double> &arg, std::normal_distribution<double> &elec):
        ionized(ions), x(0), v(random_velocity<T>(arg, gen)), accel((e / m_e) * volts * 1e2), volts_per_cm(volts), total_time(0),
        generator(gen), std_gauss(std), argon_dist(arg), elec_dist(elec)
    {
        energy = 0.5 * m_e * v * v * 6.242e18;
    }

    int index(double v) const
    {
        double E = 0.5 * m_e * v * v * 6.242e18;
        return ceil(E * 10) - 1;
    }

    inline double cross_section_e(int i) const { return sigma_e[i]; }
    inline double cross_section_p(int i) const { return sigma_p[i]; }
    inline double cross_section_i(int i) const { return sigma_i[i]; }

    double collision_probability(double u, int i) const
    {
        double prob = (u * 1e2 * cross_section_e(i)) / K_max[i];
        assert(prob < 1);
        return prob;
    }

    double ionization_probability(double u, int i) const
    {
        double prob = (u * 1e2 * cross_section_i(i)) / K_max[i];
        assert(prob < 1);
        return prob;
    }

    double next_collision(int i)
    {
        std::exponential_distribution<double> exponential(n * K_max[i]);
        return exponential(generator);
    }

    void update()
    {
        for (auto elec : ionization_electrons)
            elec->update();

        time_to_collision = next_collision(index(abs(v)));

        x += v * time_to_collision;
        v += accel * time_to_collision;

        T vm = random_velocity<T>(argon_dist, generator);

        double u = abs(v - vm);

        int xsec_index = index(u);
        
        double ion_prob = ionization_probability(u, xsec_index);
        double col_prob = collision_probability(u, xsec_index);

        double rand = drand48();

        if (rand < ion_prob)
        {
            ionization_electrons.push_back(new Electron<T>(ionized, volts_per_cm, generator, std_gauss, argon_dist, elec_dist));
            double v_length = abs(v);
            v = ((v_length - 2355743) / v_length) * v; // Only considers the first ionization energy value of Ar, 15.76 eV
            ionized++;
        }
        else if (rand < (ion_prob + col_prob))
        {
            T n = random_unit_vector<T>(std_gauss, generator);

            T a = (M_a * u) * n;
            T b = (m_e) * v;
            T c = (M_a) * vm;

            T v1 = (a + b + c) / (m_e + M_a);

            double r = cross_section_p(xsec_index) / cross_section_e(xsec_index);

            if (drand48() < r)
                v = v1;
            else
                v = abs(v1) / abs(v) * v;
        }
        
        energy = 0.5 * m_e * v * v * 6.242e18;
        total_time += time_to_collision;
    }
};

template<typename T>
void print_data_line(Electron<T> &elec, std::ofstream &file)
{
    file << elec.total_time * 1e9 << "," << elec.x * 1e6 << "," << elec.energy << "," << (elec.x / elec.total_time) << "," << elec.ionized << "\n";
    return;
}

template<>
void print_data_line<Vec>(Electron<Vec> &elec, std::ofstream &file)
{
    file << elec.total_time * 1e9 << "," << elec.x.x * 1e6 << "," << elec.x.y * 1e6 << "," << elec.x.z * 1e6 << "," << elec.energy << "," << (elec.x.x / elec.total_time) << "," << elec.ionized << "\n";
    return;
}

template<typename T>
void generate_plot(Electron<T> &elec, double cutoff, std::ofstream &file, int write_every)
{
    int simulation_step = 1;
    
    const int total_progress_steps = 10;
    std::array<bool, total_progress_steps> progress;
    progress.fill(false);
    int progress_count = 1;

    while (elec.total_time < cutoff)
    {
        elec.update();

        if (((elec.total_time / cutoff) > (float(progress_count) / total_progress_steps)) && progress[(progress_count++ - 1)] == false)
        {
            progress[progress_count - 2] = true;
            std::cout << "." << std::flush;
        }

        if (simulation_step++ % write_every != 0)
            continue;
        
        print_data_line<T>(elec, file);
    }

    elec.update();
    print_data_line<T>(elec, file);

    return;
}

#endif