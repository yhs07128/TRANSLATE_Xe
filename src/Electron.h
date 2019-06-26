#ifndef ELECTRON_H
#define ELECTRON_H

#include <random>
#include <algorithm>
#include <array>
#include <vector>
#include <chrono>
#include <assert.h>

#include "Constants.h"
#include "Utils.h"
#include "Vec.h"

template<typename S>
class Electron
{
public:
    S x, v, accel;

    double volts_per_cm;
    double time_to_collision;
    double total_time;

    double energy;
    int& ionized;
    int child_ions;

    std::vector<Electron<S>*> ionization_electrons;
    std::mt19937 &generator;

    Electron(int &ions, double volts, S velocity, std::mt19937 &gen):
        ionized(ions), child_ions(0), x(0), v(velocity), accel((e / m_e) * volts * 1e2), volts_per_cm(volts), total_time(0), generator(gen)
    {
        energy = J_to_eV(0.5 * m_e * v * v);
    }

    int index(double v) const
    {
        double E = J_to_eV(0.5 * m_e * v * v);

        if (E < 1e-3)
            E = 1e-3;

        int ind = round(167.12 * log10(E) + 501.36);
        assert(ind < 1000);

        return ind;
    }

    inline double x_sec_e(double v) const
    {
        return (is_gas) ? momentum_xsec_gas[index(v)] : effective_xsec_liquid[index(v)];
    }

    inline double x_sec_p(double v) const
    {
        return (is_gas) ? momentum_xsec_gas[index(v)] : momentum_xsec_liquid[index(v)];
    }

    inline double x_sec_i(double v) const
    {
        return ionization_xsec[index(v)];
    }

    inline double x_sec_ex(double v, int level) const
    {
        switch (level) {
            case 11:
                return excite_xsec_11[index(v)];
            case 13:
                return excite_xsec_13[index(v)];
            case 14:
                return excite_xsec_14[index(v)];
            case 15:
                return excite_xsec_15[index(v)];
            default:
                assert(false);
        }
    }

    void remove_energy(double eV)
    {
        energy -= eV;
        std::uniform_real_distribution<double> dist(8.0, 12.0);

        if (energy > 12)
            energy = dist(generator);

        if (energy < 0)
            energy = 0;
        
        v = sqrt((2 * eV_to_J(energy)) / m_e) * (v / abs(v));
    }

    double probability(double u, double x_sec)
    {
        double prob = (u * 1e2 * x_sec) / K_max;
        assert (prob < 1);
        return prob;
    }

    double next_collision(double u)
    {
        double lambda = n * K_max;
        double beta = 1 / lambda;

        std::exponential_distribution<double> exponential(lambda);
        double time_step = exponential(generator);

        return time_step;
    }

    void update()
    {
        for (auto elec : ionization_electrons)
            elec->update();

        time_to_collision = next_collision(abs(v));

        x += v * time_to_collision + 0.5 * accel * time_to_collision * time_to_collision;
        v += accel * time_to_collision;
        
        energy = J_to_eV(0.5 * m_e * v * v);

        S vm = random_velocity<S>(generator);

        double u = abs(v - vm);
        
        double col_prob = probability(u, x_sec_e(u));
        double ion_prob = probability(u, x_sec_i(u));
        
        double ex_11_prob = probability(u, x_sec_ex(u, 11));
        double ex_13_prob = probability(u, x_sec_ex(u, 13));
        double ex_14_prob = probability(u, x_sec_ex(u, 14));
        double ex_15_prob = probability(u, x_sec_ex(u, 15));

        std::uniform_real_distribution<double> rand_event(0.0, 1.0);
        double prob = rand_event(generator);

        if (prob < ion_prob) {
            if (recursive) {
                std::uniform_real_distribution<double> rand_starting_energy(2.0, 4.0);
                ionization_electrons.push_back(new Electron<S>(ionized, volts_per_cm, rand_starting_energy(generator), generator));
            }

            remove_energy(15.76);

            ionized++;
            child_ions++;
        } else if (prob < (ion_prob + ex_11_prob)) {
            remove_energy(11.68);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob)) {
            remove_energy(13.21);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob)) {
            remove_energy(14.10);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob)) {
            remove_energy(15.23);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob + col_prob)) {
            S unit_vec = random_unit_vector<S>(generator);

            S a = (M_a * u) * unit_vec;
            S b = (m_e) * v;
            S c = (M_a) * vm;

            S v1 = (a + b + c) / (m_e + M_a);

            double r = x_sec_p(u) / std::max(x_sec_e(u), x_sec_p(u));

            std::uniform_real_distribution<double> rand_momentum_change(0.0, 1.0);

            if (rand_momentum_change(generator) < r)
                v = v1;
            else
                v = abs(v1) / abs(v) * v;
        }
        
        energy = J_to_eV(0.5 * m_e * v * v);
        total_time += time_to_collision;
    }
};

template<typename T>
void print_data_line(Electron<T> &elec, std::ofstream &file)
{
    file << elec.total_time * 1e9 << "," << elec.x * 1e6 << "," << 0 << "," << 0 << "," << elec.energy << "," << (elec.x / elec.total_time) << "," << elec.ionized << "\n";
    return;
}

template<>
void print_data_line(Electron<Vec> &elec, std::ofstream &file)
{
    file << elec.total_time * 1e9 << "," << elec.x.x * 1e6 << "," << elec.x.y * 1e6 << "," << elec.x.z * 1e6 << "," << elec.energy << "," << (elec.x.x / elec.total_time) << "," << elec.ionized << "\n";
    return;
}

template<typename T>
void generate_plot(int volts, double cutoff, int number, int write_every, int group, std::string folder)
{
    int ions = 0;
    std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    Electron<T> elec(ions, volts, random_velocity<T>(generator, true), generator);
    std::ofstream file;

    file = std::ofstream("../py/simulation-runs/" + folder + "/" + std::to_string(volts) + "V - " + std::to_string(number) + ".txt");
    
    int simulation_step = 1;
    const int total_progress_steps = 100;
    std::array<bool, total_progress_steps> progress;
    progress.fill(false);
    int progress_count = 1;

    while (elec.total_time < cutoff) {
        elec.update();

        if (((elec.total_time / cutoff) > (float(progress_count) / total_progress_steps)) && progress[(progress_count++ - 1)] == false) {
            progress[progress_count - 2] = true;
            std::cout << "\rWriting " << volts << "V file group " << group << " - " << progress_count - 1 << "%" << std::flush;
        }

        if (simulation_step++ % write_every != 0)
            continue;
        
        print_data_line<T>(elec, file);
    }

    elec.update();
    print_data_line<T>(elec, file);

    file.close();

    return;
}

#endif