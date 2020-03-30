TARGET=sim
OBJECTS_DIR=obj
# as I want to support 4.8 and 5 this will set a flag for some of the mac stuff
# mainly in the types.h file for the setMacVisual which is native in Qt5
isEqual(QT_MAJOR_VERSION, 5) {
cache()
}

MOC_DIR=moc
CONFIG-=app_bundle
CONFIG+=c++11
QT+= opengl core
include ($(HOME)/NGL/UseNGL.pri)

SRC_DIR=../

SOURCES+=$$SRC_DIR/src/main.cpp
SOURCES+=$$SRC_DIR/src/NGLScene.cpp
SOURCES+=$$SRC_DIR/src/NGLSceneMouseControls.cpp
SOURCES+=$$SRC_DIR/src/Point.cpp

HEADERS+=../include/NGLScene.h
HEADERS+=../include/WindowParams.h
HEADERS+=../include/Point.h

INCLUDEPATH+=../include/

OTHER_FILES+=../shaders/*

DESTDIR=../

CONFIG += console
CONFIG -= app_bundle
LIBS+=-lGLEW

