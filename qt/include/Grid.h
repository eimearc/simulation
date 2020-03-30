#pragma once

#include <ngl/AbstractVAO.h>
#include <ngl/Vec3.h>
#include <vector>

class Grid
{
public:
    Grid()=default;
    Grid(float m_size, int m_numSteps);
    Grid(Grid &&_grid);
    Grid& operator=(Grid &&_other);
    ~Grid()=default;
   void startCoords(ngl::Vec3 &_coords);
   void draw();

private:
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec3> m_vbo;

    float m_size=1.5f;
    size_t m_numSteps=10;
    float m_stepSize;
    bool m_2d=false;

    void makeVBO();
    void makeVBOXY(ngl::Real _u, ngl::Real _v, ngl::Real _z);
    void makeVBOXZ(ngl::Real _u, ngl::Real _y, ngl::Real _v);
};
