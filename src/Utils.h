#ifndef UTILS_H
#define UTILS_H

#include <random>

#include "Constants.h"
#include "Vec.h"
#include "LUTs/Tip.h"

inline double eV_to_J(double eV)
{
    return eV * 1.60218e-19;
}

inline double J_to_eV(double J)
{
    return J * 6.242e18;
}

inline double eV_to_v(double eV)
{
    double J = eV_to_J(eV);
    return sqrt(2 * J / m_e);
}

Vec random_unit_vector(std::mt19937& gen)
{
    std::normal_distribution<double> dist(0.0, 1.0);

    Vec v(dist(gen), dist(gen), dist(gen));

    if (norm(v) != 0)
        v /= norm(v);

    return v;
}

Vec random_velocity(std::mt19937& gen, bool electron=false)
{
    if (!interactions)
        return Vec(0, 0, 0);

    double a;

    if (!electron)
        a = sqrt(k_b * T / M_a);
    else
        a = sqrt(k_b * T / m_e);

    std::gamma_distribution<double> maxwell(1.5, 1 / (2 * a * a));
    Vec v = random_unit_vector(gen) * sqrt(maxwell(gen));

    return v;
}

Vec starting_pos(std::mt19937& gen)
{
    if (uniform_field)
        return Vec(0, 0, 0);

    std::normal_distribution<double> dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    Vec random_disc(0, dist(gen), dist(gen));

    if (norm(random_disc) != 0)
        random_disc /= norm(random_disc);

    random_disc *= uniform(gen);
    random_disc *= (r_max - r_min) * spawn_width_scale;

    return Vec((z_max - z_min) * spawn_height_scale + z_min, random_disc.y, random_disc.z);
}

double shift_coord(double x)
{
    bool negative = signbit(x);
    double width = r_max - r_min;
    
    if (negative) {
        while ((x + 2 * width) < r_max)
            x += 2 * width;
    } else {
        while ((x - 2 * width) > -r_max)
            x -= 2 * width;
    }

    return x;
}

Vec local_coords(Vec pos)
{
    return Vec(pos.x, shift_coord(pos.y), shift_coord(pos.z));
}

double r(Vec local_pos)
{
    double r_val = sqrt(local_pos.y * local_pos.y + local_pos.z * local_pos.z);

    if (r_val > r_max)
        r_val = r_max;

    return r_val;
}

double z(Vec local_pos)
{
    double z_val = local_pos.x;

    if (z_val > z_max)
        z_val = z_max;

    if (z_val < 0)
        z_val = 0;

    return z_val;
}

bool hit_check(Vec pos)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    if (z_val <= z_min) 
        return true;

    if (rounded_tip) {
        if (top_radius == bot_radius) {
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && r_val <= top_radius) )
                return true;
        } else {
            if (z_val < bot_z) {
                if (r_val < bot_radius)
                    return true;
                else
                    return false;
            }
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z)) )
                return true;
        }
    } else {
        if (top_radius == bot_radius) {
            if (z_val <= top_z && r_val <= top_radius)
                return true;
        } else {
            if (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z))
                return true;
        }
    }

    return false;
}

int r_index(double r_val)
{
    return floor(r_conversion_factor * (r_val - r_min));
}

int z_index(double z_val)
{
    return floor(z_conversion_factor * (z_val - z_min));
}

Vec accel_from_E(Vec pos, double volts)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    int i = r_index(r_val);
    int j = z_index(z_val);
    
    double E_r_val = E_r[num_z_steps * i + j];
    double E_z_val = E_z[num_z_steps * i + j];

    Vec inward_vec;
    
    if (!single_tip)
        inward_vec = Vec(0, -local.y, -local.z);
    else
        inward_vec = Vec(0, -pos.y, -pos.z);

    if (norm(inward_vec) != 0)
        inward_vec /= norm(inward_vec);

    inward_vec *= E_r_val;

    Vec e_field(-E_z_val, inward_vec.y, inward_vec.z);

    return (volts / bulk_field) * (e / m_e) * e_field * 1e2;
}

#endif