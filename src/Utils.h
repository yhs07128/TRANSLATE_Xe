#include <random>

#include "Constants.h"
#include "Vec.h"

double maxwell(std::normal_distribution<double> &dist, std::mt19937 &gen)
{
    double s1 = dist(gen);
    double s2 = dist(gen);
    double s3 = dist(gen);

    double speed = sqrt(s1 * s1 + s2 * s2 + s3 * s3);

    if (drand48() < 0.5)
        return -speed;
    return speed;
}

template<typename T>
T random_velocity(std::normal_distribution<double> &dist, std::mt19937 &gen)
{
    T v = maxwell(dist, gen);
    return v;
}

template<>
Vec random_velocity<Vec>(std::normal_distribution<double> &dist, std::mt19937 &gen)
{
    Vec v(maxwell(dist, gen), maxwell(dist, gen), maxwell(dist, gen));
    return v;
}

template<typename T>
T random_unit_vector(std::normal_distribution<double> &dist, std::mt19937 &gen)
{
    if (drand48() < 0.5)
        return -1;
    return 1;
}

template<>
Vec random_unit_vector<Vec>(std::normal_distribution<double> &dist, std::mt19937 &gen)
{
    Vec v(dist(gen), dist(gen), dist(gen));
    v /= abs(v);
    return v;
}