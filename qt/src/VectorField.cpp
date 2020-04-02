#include "VectorField.h"

#include <ngl/VAOFactory.h>
#include <ngl/SimpleVAO.h>
#include <ngl/ShaderLib.h>

VectorField::VectorField(size_t _width, size_t _height, size_t _depth, float _size) :
    m_grid(_width, _height, _depth, _size)
{
    m_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);

    ngl::Vec3 position;
    ngl::Vec3 direction;
    ngl::Vec3 velocity(0.0f,0.1f,0.0f);

    ngl::Vec3 coords;
    startCoords(coords);
    float u = coords.m_x;
    float v = coords.m_y;
    float w = coords.m_z;

    const size_t &width = m_grid.width();
    const size_t &height = m_grid.height();
    const size_t &depth = m_grid.depth();
    const float &step = m_grid.stepSize();

    for(size_t i = 0; i < width; ++i)
    {
        for(size_t j = 0; j < height; ++j)
        {
            for(size_t k = 0; k < depth; ++k)
            {
                position = ngl::Vec3(u+step*i,v+step*j,w+step*k);

                float x = (rand()%10)/10.0f;
                float y = (rand()%10)/10.0f;
                float z = (rand()%10)/10.0f;
                direction = ngl::Vec3(x, y, z);

                m_points.push_back({position, direction, velocity});
            }
        }
    }

    updateVBO();
}

VectorField& VectorField::operator=(VectorField &&_other)
{
    m_points = _other.m_points;
    m_grid = std::move(_other.m_grid);
    m_vao = std::move(_other.m_vao);
    m_vbo = _other.m_vbo;
    return *this;
}

void VectorField::updateVBO()
{
    m_vbo.clear();

    for (const auto& point : m_points)
    {
        m_vbo.push_back(point.position());
        m_vbo.push_back(point.direction());
    }

    size_t size = m_vbo.size();

    m_vao->bind();
        m_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_vbo[0].m_x));
        m_vao->setNumIndices(size);
        m_vao->setVertexAttributePointer(0,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),0); // Position.
        m_vao->setVertexAttributePointer(1,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),3); // Direction.
        m_vao->setMode(GL_POINTS);
    m_vao->unbind();
}

void VectorField::update()
{
    for (auto &point : m_points)
    {
        point.update();
    }
    updateVBO();
}

void VectorField::draw()
{
    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
}

void VectorField::startCoords(ngl::Vec3 &_coords)
{
    _coords.m_x = m_grid.gridSize()/2.0f;
    _coords.m_y = -(_coords.m_x);
    _coords.m_z = -(_coords.m_x);

    _coords.m_x *= -1;
    _coords += m_grid.stepSize()/2.0f;
}
