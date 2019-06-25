#ifndef Vec_H
#define Vec_H

#include <iostream>

class Vec
{
public:
    double x, y, z;

    Vec() { }
    Vec(double a): x(a), y(0), z(0) { } // This odd constructor is specfically to deal with the x-axis acceleration on 3D simulations
    Vec(double _x, double _y, double _z): x(_x), y(_y), z(_z) { }

    inline Vec& operator+=(const Vec &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    inline Vec& operator/=(const double s)
    {
        double k = 1 / s;
        x *= k;
        y *= k;
        z *= k;
        return *this;
    }
};

inline double abs(const Vec &v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline Vec operator+(const Vec &v, const Vec &w)
{
    return Vec(v.x + w.x, v.y + w.y, v.z + w.z);
}

inline Vec operator-(const Vec &v, const Vec &w)
{
    return Vec(v.x - w.x, v.y - w.y, v.z - w.z);
}

inline double operator*(const Vec &v, const Vec &w)
{
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

inline Vec operator*(float s, const Vec &v)
{
    return Vec(s * v.x, s * v.y, s * v.z);
}

inline Vec operator*(const Vec &v, float s)
{
    return Vec(s * v.x, s * v.y, s * v.z);
}

inline Vec operator/(const Vec v, float s)
{
    return Vec(v.x / s, v.y / s, v.z / s);
}

#endif