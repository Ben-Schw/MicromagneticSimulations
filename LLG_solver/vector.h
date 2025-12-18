#ifndef VEC3_H
#define VEC3_H

#include <cmath>

/* Struct that creates a three dimensional vector.
   The struct also defines simple mathmatics connected to vectors.*/
struct Vec3 {
    double x, y, z;

    constexpr Vec3() : x(0), y(0), z(0) {}
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Vector addition, subtraction, scalar multiplication and division
    Vec3 operator+(const Vec3& b) const {
        return Vec3(x + b.x, y + b.y, z + b.z);
    }
    Vec3 operator-(const Vec3& b) const {
        return Vec3(x - b.x, y - b.y, z - b.z);
    }
    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }
    Vec3 operator/(double scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }

    // Dot product, cross product, norm, and normalization
    double dot(const Vec3& b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec3 cross(const Vec3& b) const {
        return Vec3(
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x
        );
    }

    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 normalized() const {
        double n = norm();
        return (n > 0.0) ? (*this / n) : Vec3();
    }

};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

struct Vec3i {
    int x, y, z;

    constexpr Vec3i() : x(0), y(0), z(0) {}
    constexpr Vec3i(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

    constexpr bool periodic_x() const { return x != 0; }
    constexpr bool periodic_y() const { return y != 0; }
    constexpr bool periodic_z() const { return z != 0; }

};

#endif