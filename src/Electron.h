#ifndef ELECTRON_H
#define ELECTRON_H

#include <algorithm>
#include <assert.h>
#include <random>
#include <vector>

#include "Constants.h"
#include "LUTs/CrossSections.h"
#include "LUTs/Tip.h"
#include "ProgressBar.h"
#include "Vec.h"

class Electron
{
private:
    Vec x, v, accel;

    double volts_per_cm;
    double time_to_collision;
    double total_time;

    double energy;
    int child_ions;

    std::mt19937& generator;

    int index(double v) const;

    /*
     * Given a velocity, returns the effective scattering cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @return The effective scattering cross-section (in cm^2)
     */
    inline double x_sec_e(double v) const { return (is_gas) ? momentum_xsec_gas[index(v)] : effective_xsec_liquid[index(v)]; }

    /*
     * Given a velocity, returns the momentum-transfer cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @return The momentum-transfer cross-section (in cm^2)
     */
    inline double x_sec_p(double v) const { return (is_gas) ? momentum_xsec_gas[index(v)] : momentum_xsec_liquid[index(v)]; }

    /*
     * Given a velocity, returns the ionization cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @return The ionization cross-section (in cm^2)
     */
    inline double x_sec_i(double v) const { return ionization_xsec[index(v)]; }

    /*
     * Given a velocity, returns the excitation cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @param level The excitation energy 'band'
     * @return The excitation cross-section (in cm^2)
     */
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

    double probability(double u, double x_sec);
    double next_collision(double u);

    void update_pos_vel();
    void remove_energy(double eV);
    void ionization(std::vector<Electron*>& electron_list, int& total_ionizations);
    void elastic_collision(double u, Vec vm);

public:
    Electron(double initial_time, double volts, Vec position, Vec velocity, std::mt19937& gen);

    inline Vec position() const { return x; }
    inline Vec velocity() const { return v; }
    inline double elapsed_time() const { return total_time; }
    inline double ke() const { return energy; }

    void update(std::vector<Electron*>& electron_list, int& total_ionizations);
};

void generate_plot(int volts, double cutoff, int cores, int write_every, int k, int batches, ProgressBar& bar);

#endif