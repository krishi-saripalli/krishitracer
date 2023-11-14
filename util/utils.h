#ifndef UTILS_H
#define UTILS_H




#include <vector>
#include <string>
#include <iostream>


#include <Eigen/Dense>
#include <random>
using namespace Eigen;

//initializes RNG
std::random_device randdev;
std::mt19937 generator(randdev());
std::uniform_real_distribution<> distrib(0.0, 1.0);


// Returns a random floating point number between 0.0 and 1.0 uniformly
inline float randomFloat() {
    return distrib(generator);
}


// prints the values of a Vector3f in one line
inline void printVec3(const Vector3f &vec) {
    std::cout << vec(0) << " " << vec(1) << " " << vec(2) << std::endl;
}

// Reflects a vec3
inline Vector3f reflectVec3(const Vector3f& w_i, const Vector3f& normal) {
    return w_i - (2.0f * w_i.dot(normal) * normal);
}


// TODO: Add annotation after understanding
inline Matrix3f rotateHemisphereMatrix(Vector3f normal) {
    Vector3f y = Vector3f(0.0f,1.0f,0.0f);

    Matrix3f R;
    R = Quaternionf().setFromTwoVectors(y,normal);

    return R;

}

// Schlick's approximation for the probability of reflectance
inline float schlick(float cosine, float omega_i, float omega_o) {
    float r0 = pow((omega_i - omega_o) / (omega_i + omega_o),2);
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}

#endif // UTILS_H
