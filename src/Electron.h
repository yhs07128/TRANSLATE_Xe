#ifndef ELECTRON_H
#define ELECTRON_H

#include <random>
#include <algorithm>
#include <array>
#include <vector>
#include <chrono>
#include <memory>
#include <assert.h>

#include "Constants.h"
#include "LUTs/CrossSections.h"
#include "Vec.h"
#include "Utils.h"
#include "ProgressBar.h"

class Electron
{
private:
    Vec x, v, accel;

    double volts_per_cm;
    double time_to_collision;
    double total_time;

    double energy;
    int child_ions;

    std::vector<Electron*> ionization_electrons;
    std::mt19937& generator;

    int index(double v) const
    {
        double E = J_to_eV(0.5 * m_e * v * v);

        if (E < eV_min)
            E = eV_min;

        int ind = round(((eV_steps - 1) / log10(eV_max / eV_min)) * (log10(E) - log10(eV_min)));
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
        
        v = sqrt((2 * eV_to_J(energy)) / m_e) * (v / norm(v));
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

        if (time_step > 3 * beta)
            time_step = 3 * beta;

        return time_step;
    }

    void update_pos_vel()
    {
        Vec x_old = x;
        
        if (!uniform_field) {
            x += v * time_to_collision + 0.5 * accel_from_E(x, volts_per_cm) * time_to_collision * time_to_collision;
            v += 0.5 * (accel_from_E(x_old, volts_per_cm) + accel_from_E(x, volts_per_cm)) * time_to_collision;
        } else {
            x += v * time_to_collision + 0.5 * accel * time_to_collision * time_to_collision;
            v += accel * time_to_collision;
        }

        return;
    }

    void elastic_collision(double u, Vec vm)
    {
        Vec unit_vec = random_unit_vector(generator);

        Vec a = (M_a * u) * unit_vec;
        Vec b = (m_e) * v;
        Vec c = (M_a) * vm;

        Vec v1 = (a + b + c) / (m_e + M_a);

        double r = x_sec_p(u) / std::max(x_sec_e(u), x_sec_p(u));

        std::uniform_real_distribution<double> rand_momentum_change(0.0, 1.0);

        if (rand_momentum_change(generator) < r)
            v = v1;
        else
            v = norm(v1) / norm(v) * v;

        return;
    }

    void ionization()
    {
        if (use_recursive_ionization) {
            std::uniform_real_distribution<double> rand_eV(1.0, 5.0);
            Vec near_therm_vel = random_unit_vector(generator) * eV_to_v(rand_eV(generator));
            ionization_electrons.push_back(new Electron(ionized, volts_per_cm, x, near_therm_vel, generator));
        }

        remove_energy(15.76);

        ionized++;
        child_ions++;

        return;
    }

public:
    int& ionized;

    Electron(int& ions, double volts, Vec position, Vec velocity, std::mt19937& gen):
        ionized(ions), child_ions(0), x(position), v(velocity), accel((e / m_e) * volts * 1e2, 0, 0), volts_per_cm(volts), total_time(0), generator(gen)
    {
        if (!uniform_field)
            accel = accel_from_E(x, volts_per_cm);
    
        energy = J_to_eV(0.5 * m_e * v * v);

    }

    Vec position() const { return x; }
    Vec velocity() const { return v; }
    double elapsed_time() const { return total_time; }
    double ke() const { return energy; }
    double ions() const { return ionized; }
    
    void update()
    {
        for (auto elec : ionization_electrons) {
            elec->update();
            if (!uniform_field)
                ionization_electrons.erase(std::remove_if(ionization_electrons.begin(), ionization_electrons.end(), [](auto const& i){ return hit_check(i->x); }), ionization_electrons.end());
        }
        
        time_to_collision = next_collision(norm(v));
        update_pos_vel();

        if (!interactions) {
            if (norm(v) > 1000) {
                v /= norm(v);
                v *= 1000;
            }
        }
        
        energy = J_to_eV(0.5 * m_e * v * v);

        Vec vm = random_velocity(generator);
        double u = norm(v - vm);
        
        double col_prob = probability(u, x_sec_e(u));
        double ion_prob = probability(u, x_sec_i(u));

        double ex_11_prob = probability(u, x_sec_ex(u, 11));
        double ex_13_prob = probability(u, x_sec_ex(u, 13));
        double ex_14_prob = probability(u, x_sec_ex(u, 14));
        double ex_15_prob = probability(u, x_sec_ex(u, 15));

        assert(col_prob + ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob <= 1);

        if (!interactions) {
            col_prob = 0;
            ion_prob = 0;
            ex_11_prob = 0;
            ex_13_prob = 0;
            ex_14_prob = 0;
            ex_15_prob = 0;
        }
        
        std::uniform_real_distribution<double> rand_event(0.0, 1.0);
        double prob = rand_event(generator);
        
        if (prob < ion_prob) {
            ionization();
        } else if (prob < (ion_prob + ex_11_prob)) {
            remove_energy(11.68);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob)) {
            remove_energy(13.21);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob)) {
            remove_energy(14.10);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob)) {
            remove_energy(15.23);
        } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob + col_prob)) {
            elastic_collision(u, vm);
        }
        
        energy = J_to_eV(0.5 * m_e * v * v);
        total_time += time_to_collision;
    }
};

#endif