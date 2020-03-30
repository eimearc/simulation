TARGET=tests

include ($(HOME)/NGL/UseNGL.pri)

SRC_DIR=./
SOURCES+=$$SRC_DIR/PointTests.cpp

HEADER_DIR=../include/
HEADERS+=$$HEADER_DIR/Point.h
INCLUDEPATH+=$$HEADER_DIR

LIBS+=-L/usr/local/lib -lgtest -lgtest_main -pthread
