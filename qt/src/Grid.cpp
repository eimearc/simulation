#include "Grid.h"

#include <ngl/SimpleVAO.h>
#include <ngl/VAOFactory.h>

Grid::Grid(float _size, int _numSteps)
{
    m_size = _size;
    m_numSteps = _numSteps;
    m_stepSize = _size/static_cast<float>(_numSteps);

    makeVBO();

    size_t size = m_vbo.size();
    m_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
    m_vao->bind();

    m_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_vbo[0].m_x));
    m_vao->setNumIndices(size);
    m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);

    glPointSize(5);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(50);

    m_vao->setMode(GL_LINES);

    m_vao->unbind();
}

Grid::Grid(Grid &&_grid)
{
    m_vao = std::move(_grid.m_vao);
    m_vbo = _grid.m_vbo;
    m_size = _grid.m_size;
    m_numSteps = _grid.m_numSteps;
    m_stepSize = _grid.m_stepSize;
    m_2d = _grid.m_2d;
}

Grid& Grid::operator=(Grid &&_grid)
{
    m_vao = std::move(_grid.m_vao);
    m_vbo = _grid.m_vbo;
    m_size = _grid.m_size;
    m_numSteps = _grid.m_numSteps;
    m_stepSize = _grid.m_stepSize;
    m_2d = _grid.m_2d;
    return *this;
}

void Grid::startCoords(ngl::Vec3 &_coords)
{
    _coords.m_x = m_size/2.0f;
    _coords.m_y = -(_coords.m_x);
    _coords.m_z = -(_coords.m_x);
}

void Grid::makeVBO()
{
    ngl::Vec3 pos;
    startCoords(pos);

    float u = pos.m_x;
    float v = pos.m_y;

    if (m_2d)
    {
        float d = 0;
        for(size_t i=0; i<=m_numSteps; ++i)
        {
            makeVBOXY(u,v,d);
            v+=m_stepSize;
        }
    }
    else
    {
        for(size_t i=0; i<=m_numSteps; ++i)
        {
            float d = -u;
            for (size_t j=0; j<=m_numSteps; ++j)
            {
                makeVBOXY(u,v,d);

                d+=m_stepSize;
            }
            v+=m_stepSize;
        }

        v = -u;
        for(size_t i=0; i<=m_numSteps; ++i)
        {
            float d = -u;
            for (size_t j=0; j<=m_numSteps; ++j)
            {
                makeVBOXZ(u,d,v);

                d+=m_stepSize;
            }
            v+=m_stepSize;
        }
    }
}

void Grid::makeVBOXY(ngl::Real _u, ngl::Real _v, ngl::Real _z)
{
    m_vbo.push_back({-_u, _v, _z}); // Left vert
    m_vbo.push_back({_u, _v, _z});  // Right vert
    m_vbo.push_back({_v, _u, _z});  // Top vert
    m_vbo.push_back({_v, -_u, _z}); // Bottom vert
}

void Grid::makeVBOXZ(ngl::Real _u, ngl::Real _y, ngl::Real _v)
{
    m_vbo.push_back({-_u, _y, _v}); // Left vert
    m_vbo.push_back({_u, _y, _v});  // Right vert
    m_vbo.push_back({_v, _y, _u});  // Top vert
    m_vbo.push_back({_v, _y, -_u}); // Bottom vert
}

void Grid::draw()
{
    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
}
