#pragma once

#include <ngl/Vec3.h>
#include <Point.h>
#include <vector>
#include <ngl/AbstractVAO.h>

class VectorField
{
public:
    VectorField()=default;
    VectorField(size_t _width, size_t _height, size_t _depth, float _size, size_t _steps);
    ~VectorField()=default;
    VectorField& operator=(VectorField &&_other);
    void startCoords(ngl::Vec3 &_coords);
    void update();
    void draw();

private:
    std::vector<Point> m_points;
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec3> m_vbo;
    size_t m_width;
    size_t m_height;
    size_t m_depth;
    float m_size;
    size_t m_steps;
    float m_stepSize;

    void updateVBO();
};
