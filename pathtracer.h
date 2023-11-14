#ifndef PATHTRACER_H
#define PATHTRACER_H

#include <QImage>

#include "scene/scene.h"

using namespace Eigen;


class PathTracer
{
public:
    PathTracer(int width, int height);
    int sampleCount;

    void traceScene(QRgb *imageData, const Scene &scene);

private:
    int m_width, m_height;


    void toneMap(QRgb *imageData, std::vector<Vector3f> &intensityValues);
    Vector3f radiance(const Ray& r, const Scene& scene, int depth);
    Vector3f tracePixel(int x, int y, const Scene &scene, const Matrix4f &invViewMatrix);
    std::tuple<bool,IntersectionInfo> traceRay(const Ray& r, const Scene &scene);
};

#endif // PATHTRACER_H
