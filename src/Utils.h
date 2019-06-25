#include <random>

#include "Constants.h"
#include "Vec.h"

inline double eV_to_J(double eV)
{
    return eV * 1.60218e-19;
}

inline double J_to_eV(double J)
{
    return J * 6.242e18;
}

double maxwell(std::mt19937 &gen, double electron)
{
    std::normal_distribution<double> dist;

    if (electron)
        dist = std::normal_distribution<double>(0.0, sqrt((k_b * T) / (2 * m_e)));
    else
        dist = std::normal_distribution<double>(0.0, sqrt((k_b * T) / (2 * M_a)));

    double s1 = dist(gen);
    double s2 = dist(gen);
    double s3 = dist(gen);

    double speed = sqrt(s1 * s1 + s2 * s2 + s3 * s3);

    if (drand48() < 0.5)
        return -speed;
    return speed;
}

template<typename T>
T random_velocity(std::mt19937 &gen, double electron=false)
{
    T v = maxwell(gen, electron);
    return v;
}

template<>
Vec random_velocity<Vec>(std::mt19937 &gen, double electron)
{
    Vec v(maxwell(gen, electron), maxwell(gen, electron), maxwell(gen, electron));
    return v;
}

template<typename T>
T random_unit_vector(std::mt19937 &gen)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(gen) < 0.5)
        return -1;
    return 1;
}

template<>
Vec random_unit_vector<Vec>(std::mt19937 &gen)
{
    std::normal_distribution<double> dist(0.0, 1.0);
    Vec v(dist(gen), dist(gen), dist(gen));
    v /= abs(v);
    return v;
}