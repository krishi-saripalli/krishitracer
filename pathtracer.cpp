#include "pathtracer.h"
#include "solver/solver.h"

#include <iostream>

#include <Eigen/Dense>
#include <omp.h>

#include <util/CS123Common.h>
#include <util/utils.h>

using namespace Eigen;

PathTracer::PathTracer(int width, int height)
    : m_width(width), m_height(height)
{
}

//Top-level function for rendering scene
void PathTracer::traceScene(QRgb *imageData, const Scene& scene)
{
    std::vector<Vector3f> intensityValues(m_width * m_height);
    Matrix4f invViewMat = (scene.getCamera().getScaleMatrix() * scene.getCamera().getViewMatrix()).inverse();

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int y = 0; y < m_height; ++y) {

        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            intensityValues[offset] = tracePixel(x, y, scene, invViewMat);

        }
    }

    toneMap(imageData, intensityValues);
}


//Traces K rays through a pixel, where K is specified by the user
Vector3f PathTracer::tracePixel(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix)
{
    int K = sampleCount;
    Vector3f total_radiance = Vector3f(0,0,0);


   for (int i=0; i < K; i++) {
       float rand_float = randomFloat(); //between 0.0 and 1.0
       Vector3f p(0, 0, 0);
       Vector3f d(((2.f*float(x) + rand_float) / m_width) - 1.0f, 1.0f - ( (2.f*float(y) + rand_float) / m_height), -1);
       d.normalize();

       Ray r(p, d);
       r = r.transform(invViewMatrix);

       //trace ray
       total_radiance += radiance(r,scene,0);
   }
  return total_radiance/K;
}

// Traces Ray and return whether there was a intersection and the corresponding IntersectionInfo
// Computes the outgoing radiance at x along -w (reflected incoming ray)
Vector3f PathTracer::radiance(const Ray& r, const Scene& scene, int depth ) {

    IntersectionInfo intersectionInfo;
    Ray ray(r);
    bool intersected = scene.getIntersection(ray,&intersectionInfo);

    Vector3f L_total = Vector3f(0.0f,0.0f,0.0f);
    if (intersected) {
        const Triangle *tri = static_cast<const Triangle *>(intersectionInfo.data);
        const tinyobj::material_t& material = tri->getMaterial();
        const tinyobj::real_t *e = material.emission;
        BRDF_TYPE brdfType = BRDFType(material);

        if (brdfType == BRDF_MIRROR || brdfType == BRDF_DIELECTRIC || tri->isEmissive()) {
             L_total += Vector3f(e[0],e[1],e[2]);
        }

        float probRussian = 0.90f;

        // russian roulette (but we bounce at least 1 time)
        if (probRussian > randomFloat() || depth == 0) {


            //get BRDF, w_i, object normal and sample next direction
            Vector3f w_i = r.d;
            Vector3f normal = intersectionInfo.object->getNormal(intersectionInfo);
            std::tuple<Vector3f,float> newDirection = sampleNextDir(brdfType,scene.m_globalData, material,w_i,normal);
            Vector3f w_o = std::get<0>(newDirection);
            float pdf = std::get<1>(newDirection);

             // transform w_o if it is not a reflection/refraction
            if (brdfType != BRDF_MIRROR && brdfType != BRDF_DIELECTRIC) {
                w_o = rotateHemisphereMatrix(normal)*w_o;
            }


            Vector3f brdf = BRDF(brdfType,w_i,w_o,normal,pdf,scene.m_globalData,material);

            // Direct lighting
              if (brdfType != BRDF_MIRROR && brdfType != BRDF_DIELECTRIC) {
              L_total +=  directLighting(scene, w_i, intersectionInfo.hit, normal, brdf);

              }

             //trace new ray recursivley
            Ray newRay = Ray(intersectionInfo.hit,w_o);
            Vector3f outgoingRadiance = radiance(newRay,scene,depth + 1);

            // Indirect lighting
            L_total += (outgoingRadiance.cwiseProduct(brdf)) * (-w_i).dot(normal)*(1.0f/(pdf*probRussian));
        }

    }

    return L_total;
}


// Tone map using the Reinhard operator
void PathTracer::toneMap(QRgb *imageData, std::vector<Vector3f> &intensityValues) {
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            float red =  intensityValues[offset][0]/(1.0f+intensityValues[offset][0])*255.0f;
            float green = intensityValues[offset][1]/(1.0f+intensityValues[offset][1])*255.0f;
            float blue =  intensityValues[offset][2]/(1.0f+intensityValues[offset][2]) *255.0f;
            imageData[offset] = qRgb(red,green,blue);

        }
    }

}
