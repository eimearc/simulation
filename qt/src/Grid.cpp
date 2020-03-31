#include "Grid.h"

#include <ngl/SimpleVAO.h>
#include <ngl/VAOFactory.h>
#include <iostream>

Grid::Grid(size_t _width, size_t _height, size_t _depth, float _size)
{
    m_size = _size;
    m_width = _width;
    m_height = _height;
    m_depth = _depth;
    m_stepSize = m_size/m_width;

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

Grid& Grid::operator=(Grid &&_grid)
{
    m_vao = std::move(_grid.m_vao);
    m_vbo = _grid.m_vbo;
    m_size = _grid.m_size;
    m_width = _grid.m_width;
    m_height = _grid.m_height;
    m_depth = _grid.m_depth;
    m_stepSize = _grid.m_stepSize;
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
    float w = pos.m_z;

    for (size_t j=0; j<=m_depth; ++j)
    {
        // Vertical lines - along x
        for (size_t i=0; i<=m_width; ++i)
        {
            m_vbo.push_back({u,v,w});
            m_vbo.push_back({u,v+m_height*m_stepSize,w});

            u -= m_stepSize;
        }
        w += m_stepSize;
        u = pos.m_x;
    }

    u = pos.m_x;
    v = pos.m_y;
    w = pos.m_z;

    for (size_t j=0; j<=m_depth; ++j)
    {
        // Horizontal lines - varying y
        for (size_t i=0; i<=m_height; ++i)
        {
            m_vbo.push_back({u,v,w});
            m_vbo.push_back({u-m_width*m_stepSize,v,w});
            v += m_stepSize;
        }
        w += m_stepSize;
        v = pos.m_y;
    }

    u = pos.m_x;
    v = pos.m_y;
    w = pos.m_z;

    for (size_t j=0; j<=m_width; ++j)
    {
        // Horizontal lines - varying z
        for (size_t i=0; i<=m_height; ++i)
        {
            m_vbo.push_back({u,v,w});
            m_vbo.push_back({u,v,w+m_depth*m_stepSize});

            v += m_stepSize;
        }
        u -= m_stepSize;
        v = pos.m_y;
    }

    for(size_t i=0; i<=m_height; ++i)
    {
//        for (size_t j=0; j<=m_depth; ++j)
//        {
//            makeVBOXY(u,v,w);

//            w+=m_stepSize;
//        }
//        v+=m_stepSize;
    }

//    v = -u;
//    for(size_t i=0; i<=m_width; ++i)
//    {
//        float d = -u;
//        for (size_t j=0; j<=m_height; ++j)
//        {
//            makeVBOXZ(u,d,v);

//            d+=m_stepSize;
//        }
//        v+=m_stepSize;
//    }
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
