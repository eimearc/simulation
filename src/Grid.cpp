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

    makeInnerVBO();
    makeOuterVBO();

    size_t size = m_inner_vbo.size();
    m_inner_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
    m_inner_vao->bind();
    m_inner_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_inner_vbo[0].m_x));
    m_inner_vao->setNumIndices(size);
    m_inner_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
    glPointSize(5);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1);
    m_inner_vao->setMode(GL_LINES);
    m_inner_vao->unbind();

    size = m_outer_vbo.size();
    m_outer_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
    m_outer_vao->bind();
    m_outer_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_outer_vbo[0].m_x));
    m_outer_vao->setNumIndices(size);
    m_outer_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
    glPointSize(5);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(50);
    m_outer_vao->setMode(GL_LINES);
    m_outer_vao->unbind();
}

Grid& Grid::operator=(Grid &&_grid)
{
    m_inner_vao = std::move(_grid.m_inner_vao);
    m_inner_vbo = _grid.m_inner_vbo;
    m_outer_vao = std::move(_grid.m_outer_vao);
    m_outer_vbo = _grid.m_outer_vbo;
    m_size = _grid.m_size;
    m_width = _grid.m_width;
    m_height = _grid.m_height;
    m_depth = _grid.m_depth;
    m_stepSize = _grid.m_stepSize;
    return *this;
}

void Grid::startCoords(ngl::Vec3 &_coords) const
{
    _coords.m_x = m_size/2.0f;
    _coords.m_y = -(_coords.m_x);
    _coords.m_z = -(_coords.m_x);
}

float Grid::gridSize() const
{
    return m_size;
}

float Grid::stepSize() const
{
    return m_stepSize;
}

size_t Grid::width() const
{
    return m_width;
}

size_t Grid::height() const
{
    return m_height;
}

size_t Grid::depth() const
{
    return m_depth;
}

void Grid::makeInnerVBO()
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
            m_inner_vbo.push_back({u,v,w});
            m_inner_vbo.push_back({u,v+m_height*m_stepSize,w});

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
            m_inner_vbo.push_back({u,v,w});
            m_inner_vbo.push_back({u-m_width*m_stepSize,v,w});
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
            m_inner_vbo.push_back({u,v,w});
            m_inner_vbo.push_back({u,v,w+m_depth*m_stepSize});

            v += m_stepSize;
        }
        u -= m_stepSize;
        v = pos.m_y;
    }
}

void Grid::makeOuterVBO()
{
    ngl::Vec3 pos;
    startCoords(pos);

    float u = pos.m_x;
    float v = pos.m_y;

    // Outer
    m_outer_vbo.push_back({u,v,0});
    m_outer_vbo.push_back({u,-v,0});
    m_outer_vbo.push_back({-u,v,0});
    m_outer_vbo.push_back({-u,-v,0});

    m_outer_vbo.push_back({-u,v,0});
    m_outer_vbo.push_back({u,v,0});
    m_outer_vbo.push_back({-u,-v,0});
    m_outer_vbo.push_back({u,-v,0});

    // Inner
    u-=m_stepSize;
    v+=m_stepSize;
    m_outer_vbo.push_back({u,v,0});
    m_outer_vbo.push_back({u,-v,0});
    m_outer_vbo.push_back({-u,v,0});
    m_outer_vbo.push_back({-u,-v,0});

    m_outer_vbo.push_back({-u,v,0});
    m_outer_vbo.push_back({u,v,0});
    m_outer_vbo.push_back({-u,-v,0});
    m_outer_vbo.push_back({u,-v,0});
}

void Grid::drawOuter() const
{
    m_outer_vao->bind();
    m_outer_vao->draw();
    m_outer_vao->unbind();
}

void Grid::drawInner() const
{
    m_inner_vao->bind();
    m_inner_vao->draw();
    m_inner_vao->unbind();
}

bool Grid::operator==(const Grid &_other) const
{
    bool result = true;

    result &= (m_inner_vbo == _other.m_inner_vbo);
    result &= (m_outer_vbo == _other.m_outer_vbo);
    result &= (m_width == _other.m_width);
    result &= (m_height == _other.m_height);
    result &= (m_depth == _other.m_depth);
    result &= (m_size == _other.m_size);
    result &= (m_stepSize == _other.m_stepSize);

    return result;
}
