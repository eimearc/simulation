TARGET=tests

include ($(HOME)/NGL/UseNGL.pri)

TEST_DIR=./
SRC_DIR=../src/
SOURCES+=$$TEST_DIR/PointTests.cpp $$SRC_DIR/Point.cpp
SOURCES+=$$TEST_DIR/VectorFieldTests.cpp $$SRC_DIR/VectorField.cpp
SOURCES+=$$TEST_DIR/GridTests.cpp $$SRC_DIR/Grid.cpp

HEADER_DIR=../include/
HEADERS+=$$HEADER_DIR/Point.h
HEADERS+=$$HEADER_DIR/Grid.h
HEADERS+=$$HEADER_DIR/VectorField.h
INCLUDEPATH+=$$HEADER_DIR

#LIBS+=-lGLEW
LIBS+=$$system(pkg-config --libs glfw3)
LIBS+=-L/usr/local/lib -lgtest -lgtest_main -pthread
