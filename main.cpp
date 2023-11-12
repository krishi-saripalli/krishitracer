#include <QCoreApplication>
#include <QCommandLineParser>

#include <iostream>

#include "pathtracer.h"
#include "scene/scene.h"

#include <QImage>

#include "util/CS123Common.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("scene", "Scene file to be rendered");
    parser.addPositionalArgument("output", "Image file to write the rendered image to");
    parser.addPositionalArgument("sampleCount", "Number of rays to be traced through each pixel");


    parser.process(a);

    const QStringList args = parser.positionalArguments();
    if(args.size() < 2 || args.size() > 3) {
        std::cerr << "Error: Wrong number of arguments" << std::endl;
        a.exit(1);
        return 1;
    }
    QString scenefile = args[0];
    QString output = args[1];
    QString sampleCount = args.size() == 3 ? args[2] : QString("100"); //default 100 samples if not specified

    QImage image(IMAGE_WIDTH, IMAGE_HEIGHT, QImage::Format_RGB32);

    Scene *scene;
    if(!Scene::load(scenefile, &scene)) {
        std::cerr << "Error parsing scene file " << scenefile.toStdString() << std::endl;
        a.exit(1);
        return 1;
    }

    PathTracer tracer(IMAGE_WIDTH, IMAGE_HEIGHT);

    QRgb *data = reinterpret_cast<QRgb *>(image.bits());

    bool ok;
    // Try to convert the sampleCount QString to an int
    int K = sampleCount.toInt(&ok);

    if (ok && K >= 1 ) {
        std::cout << "Tracing " << sampleCount.toStdString() << " rays per pixel" << std::endl;

    } else {
        // Conversion failed
        std::cerr << "Invalid sampleCount " << scenefile.toStdString() << std::endl;
        a.exit(1);
        return 1;
    }

    tracer.sampleCount = K;
    tracer.traceScene(data, *scene);

    delete scene;

    bool success = image.save(output);
    if(!success) {
        success = image.save(output, "PNG");
    }
    if(success) {
        std::cout << "Wrote rendered image to " << output.toStdString() << std::endl;
    } else {
        std::cerr << "Error: failed to write image to " << output.toStdString() << std::endl;
    }
    a.exit();
}
