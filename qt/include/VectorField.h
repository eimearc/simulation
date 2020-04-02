#pragma once

#include <ngl/Vec3.h>
#include <Point.h>
#include <vector>
#include <ngl/AbstractVAO.h>
#include <gtest/gtest.h>
#include "Grid.h"

class VectorField
{
public:
    VectorField()=default;
    VectorField(size_t _width, size_t _height, size_t _depth, float _size);
    ~VectorField()=default;
    VectorField& operator=(VectorField &&_other);
    void startCoords(ngl::Vec3 &_coords);
    void update();
    void draw();

private:
    std::vector<Point> m_points;
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec3> m_vbo;

    Grid m_grid;

    void updateVBO();

    FRIEND_TEST(VectorField, ctor);
    FRIEND_TEST(VectorField, startCoords);
    FRIEND_TEST(VectorField, startCoords);
};
