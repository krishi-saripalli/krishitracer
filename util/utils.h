#ifndef UTILS_H
#define UTILS_H




#include "scene/shape/triangle.h"
#include <vector>
#include <string>
#include <iostream>


#include <Eigen/Dense>
#include <random>
using namespace Eigen;

//initializes RNG
thread_local std::random_device randdev;
thread_local std::mt19937 generator(randdev());
thread_local std::uniform_real_distribution<> distrib(0.0, 1.0);


// Returns a random floating point number between 0.0 and 1.0 uniformly
static inline float randomFloat() {
    return distrib(generator);
}


// prints the values of a Vector3f in one line
static inline void printVec3(const Vector3f &vec) {
    std::cout << vec(0) << " " << vec(1) << " " << vec(2) << std::endl;
}

// Reflects a vec3
static inline Vector3f reflectVec3(const Vector3f& w_i, const Vector3f& normal) {
    return w_i - (2.0f * w_i.dot(normal) * normal);
}


// TODO: Add annotation after understanding
static inline Matrix3f rotateHemisphereMatrix(const Vector3f& normal) {
    Vector3f y = Vector3f(0.0f,1.0f,0.0f);

    Matrix3f R;
    R = Quaternionf().setFromTwoVectors(y,normal);

    return R;

}

// Schlick's approximation for the probability of reflectance
static inline float schlick(float cosine, float omega_i, float omega_o) {
    float r0 = pow((omega_i - omega_o) / (omega_i + omega_o),2);
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}


// Samples a random point on a triangle
static inline Vector3f sampleTriangle(Triangle* triangle) {
    Vector3<Vector3f> verts = triangle->getVertices();
    Vector3f ba = verts[1]-verts[0];
    Vector3f ca = verts[2]-verts[0];
    float random1 = randomFloat();
    float random2 = randomFloat();
    if (random1 + random2 > 1.0f) {
        random1 = 1.0f - random1;
        random2 = 1.0f - random2;
     }
    Vector3f sample = verts[0] + random1*ba + random2*ca;
    return sample;
}

// All code for ACES tone mapping
// Define the ACES input and output matrices as Matrix3f
static Matrix3f createAcesInputMatrix() {
    Matrix3f mat;
    mat << 0.59719f, 0.35458f, 0.04823f,
           0.07600f, 0.90834f, 0.01566f,
           0.02840f, 0.13383f, 0.83777f;
    return mat;
}

static Matrix3f createAcesOutputMatrix() {
    Matrix3f mat;
    mat << 1.60475f, -0.53108f, -0.07367f,
           -0.10208f,  1.10813f, -0.00605f,
           -0.00327f, -0.07276f,  1.07602f;
    return mat;
}

static Matrix3f aces_input_matrix = createAcesInputMatrix();
static Matrix3f aces_output_matrix = createAcesOutputMatrix();

static inline Vector3f rtt_and_odt_fit(Vector3f v)
{
    Vector3f a = v.cwiseProduct(v + Vector3f(0.0245786f, 0.0245786f, 0.0245786f)) - Vector3f(0.000090537f, 0.000090537f, 0.000090537f);
    Vector3f b = v.cwiseProduct(0.983729f * v + Vector3f(0.4329510f, 0.4329510f, 0.4329510f)) + Vector3f(0.238081f, 0.238081f, 0.238081f);
    return a.cwiseQuotient(b);
}

// Returns a tone-mapped RGB tuple given an intensity vector
static inline Vector3f acesFitted(Vector3f v)
{
    v = aces_input_matrix * v;
    v = rtt_and_odt_fit(v);
    v = aces_output_matrix * v;

    // Scale and clamp the values to [0, 255] range
    v *= 255.0f;
    v = v.cwiseMin(Vector3f(255.0f, 255.0f, 255.0f)).cwiseMax(Vector3f(0.0f, 0.0f, 0.0f));

    return v;
}
#endif // UTILS_H
