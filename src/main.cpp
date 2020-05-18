// Taken from github.com/NCCA/SimpleNGL
#include "NGLScene.h"
#include <QtGui/QGuiApplication>
#include <iostream>
#include <gflags/gflags.h>

DECLARE_bool(colour);
DECLARE_int32(resolution);
DECLARE_int32(num_particles);
DECLARE_double(viscosity);

int main(int argc, char** argv)
{
    gflags::SetUsageMessage("Simulation program for MAC Grid fluids.");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    QGuiApplication app(argc, argv);
    QSurfaceFormat format;
    format.setSamples(4);
    // #if defined(__APPLE__)
    //   format.setMajorVersion(4);
    //   format.setMinorVersion(1);
    // #else
    format.setMajorVersion(4);
    format.setMinorVersion(5);
    // #endif
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setDepthBufferSize(24);
    QSurfaceFormat::setDefaultFormat(format);
    NGLScene window;
    window.resize(1024, 720);
    window.show();
    return app.exec();
}
