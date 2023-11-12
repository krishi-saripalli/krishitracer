#ifndef SOLVER_H
#define SOLVER_H




#include <vector>
#include <string>
#include <iostream>


#include <Eigen/Dense>
#include <random>
#include "util/CS123SceneData.h"
#include "util/tiny_obj_loader.h"
#include "util/utils.h"

using namespace Eigen;

enum BRDF_TYPE {BRDF_DIFFUSE = 1, BRDF_PHONG = 2,BRDF_MIRROR = 3,BRDF_REFRACT = 4};


// Return the BRDF type given the material information at the point of intersection
BRDF_TYPE BRDFType( const tinyobj::material_t& mat) {


    const tinyobj::real_t *s = mat.specular;
    const tinyobj::real_t *d = mat.diffuse;
    Vector3f vSpecular = Vector3f(s[0],s[1],s[2]);
    Vector3f vDiffuse = Vector3f(d[0],d[1],d[2]);

    BRDF_TYPE type = BRDF_DIFFUSE; //default diffuse

    if (vDiffuse.norm() > 0.0f) {
        type = BRDF_DIFFUSE;
    }
    if ( vSpecular.norm() > 0.0f) {
        type = BRDF_PHONG;
    }
    if (mat.illum == 5) {
        type = BRDF_MIRROR;
    }
    //TODO: Find condition for Refractive Surface
    if (1==2) {
        type = BRDF_REFRACT;
    }

    return type;

}


// Returns a Diffuse BRDF and cosine weighted sampling
std::tuple<Vector3f,float> sampleDiffuseDir(const Vector3f w_i, const Vector3f normal) {
    Vector3f direction;

    float prob = (1.0f/M_PI) * abs(w_i.dot(normal));

    float random_f = randomFloat();
    float random_f_2 = randomFloat();

    //sample angle according to cosine weighting
    float theta = acos(sqrt(random_f)), phi = 2.0f*M_PI*random_f_2;
    //convert to rectangular coords
    float x = sin(theta)*cos(phi);
    float y = cos(theta);
    float z = sin(theta)*sin(phi);


    direction = Vector3f(x,y,z).normalized();
    return std::tuple<Vector3f,float>(direction,prob);


}

// Returns w_o and probability biased towards the normal (dull mirror/glossy effect)
// shininess is the specular exponent from Phong
std::tuple<Vector3f,float> samplePhongDir(const Vector3f w_i, const Vector3f normal, const float ks, const float kd, const float shininess) {

    Vector3f direction;
    float prob;

    Vector3f reflection = reflectVec3(w_i,normal);
    float random_f = randomFloat();
    float random_f_2 = randomFloat();
    float random_f_3 = randomFloat();
    float theta, phi;


    //diffuse contribution
    if (random_f > kd) {


        //cosine weighted pdf
        prob = (1.0f/M_PI)* abs(w_i.dot(normal));

        //sample angle based on pdf
        theta = acos(sqrt(random_f_2)), phi = 2.0f*M_PI*random_f_3;


    }

    //specular contribution
    else if (kd <= random_f && random_f < kd+ks) {


        //sample angle based on pdf
        float power = powf(random_f_2,(1.0f/(shininess+1.0f)));
        theta = acos(power), phi = 2.0f*M_PI*random_f_3;

        //convert to rectangular coords
        float x = sin(theta)*cos(phi);
        float y = cos(theta);
        float z = sin(theta)*sin(phi);
        Vector3f w_o = Vector3f(x,y,z).normalized();

        //specular pdf
        //find angle alpha between the (transformed!) w_o and perfect reflection
        float cosine = std::clamp((rotateHemisphereMatrix(normal)*w_o).dot(reflection),0.0f,1.0f);
        prob = (shininess+2)/(2.0f*M_PI) * powf(cosine,shininess);
        return std::tuple<Vector3f,float>(w_o,prob);

    }

    //no contribution
    else {
        return std::tuple<Vector3f,float>(Vector3f{0.0f,0.0f,0.0f},1.0f);
    }

    float x = sin(theta)*cos(phi);
    float y = cos(theta);
    float z = sin(theta)*sin(phi);
    direction = Vector3f(x,y,z).normalized();


    return std::tuple<Vector3f,float>(direction,prob);


}

// Returns a Mirror BRDF with prob 1.0
std::tuple<Vector3f,float> sampleMirrorDir(Vector3f w_i, Vector3f normal) {
    Vector3f direction = reflectVec3(w_i,normal).normalized();
    float prob = 1.0f;


    return std::tuple<Vector3f,float>(direction,prob);
}


// Returns the BRDF Vector quantity given a BRDF type
Vector3f BRDF( const BRDF_TYPE brdfType , const Vector3f w_i, const Vector3f w_o, const Vector3f normal, const CS123SceneGlobalData& sceneData, const tinyobj::material_t& mat) {
    const tinyobj::real_t *d = mat.diffuse;
    const tinyobj::real_t *s = mat.specular;
    Vector3f reflection = reflectVec3(w_i,normal);
    float cosine = std::clamp(w_o.dot(reflection),0.0f,1.0f);
    float t = sceneData.ks/(sceneData.ks + sceneData.kd);
    Vector3f spec = t*Vector3f(s[0],s[1],s[2])* (mat.shininess + 2/(2.0f*M_PI))*powf(cosine,mat.shininess);

    switch(brdfType) {
      case BRDF_DIFFUSE:
        return Vector3f(d[0],d[1],d[2])/M_PI;
        break;
      case BRDF_PHONG:
        return   (1.0f-t)*Vector3f(d[0],d[1],d[2])/M_PI + t*Vector3f(s[0],s[1],s[2])* (mat.shininess + 2/(2.0f*M_PI))*powf(cosine,mat.shininess);
        break;
      case BRDF_MIRROR:
        //TODO: Look at the 1/cos(theta) constant to multiply here and understand if we need to add it
        return (1.0f/reflection.dot(normal))*Vector3f(1.0f,1.0f,1.0f);

      //TODO: Add BRDF_REFRACT case
      default:
        return Vector3f(d[0],d[1],d[2])/M_PI;


    }

}



// Returns the next sampled outgoing direction and the prob of choosing it based on the given BRDF type
std::tuple<Vector3f,float> sampleNextDir(BRDF_TYPE brdfType, const CS123SceneGlobalData& sceneData,  const tinyobj::material_t& material,  const Vector3f& w_i, const Vector3f& normal) {

    switch(brdfType) {
      case BRDF_DIFFUSE:
        return sampleDiffuseDir(w_i,normal);
        break;
      case BRDF_PHONG:
        return samplePhongDir(w_i, normal, sceneData.ks,sceneData.kd, material.shininess);
        break;
      case BRDF_MIRROR:
        return sampleMirrorDir(w_i, normal);

      //TODO: Add BRDF_REFRACT case
      default:
        return sampleDiffuseDir(w_i,normal);


    }


}






#endif // SOLVER_H
