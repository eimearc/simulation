#pragma once

#include <ngl/AbstractVAO.h>
#include <ngl/Vec3.h>
#include <vector>

class Grid
{
public:
    Grid()=default;
    Grid(size_t m_width, size_t m_height, size_t m_depth, float m_size);
    Grid(Grid &&_grid);
    Grid& operator=(Grid &&_other);
    ~Grid()=default;
   void startCoords(ngl::Vec3 &_coords);
   void draw();

private:
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec3> m_vbo;

    size_t m_width;
    size_t m_height;
    size_t m_depth;
    float m_size;
    float m_stepSize;

    void makeVBO();
};
