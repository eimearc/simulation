#include "MAC.h"

#include <gflags/gflags.h>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdlib>
#include <ngl/SimpleVAO.h>
#include <ngl/NGLInit.h>
#include <algorithm>
#include <functional>
#include <ngl/ShaderLib.h>

DEFINE_bool(colour, false, "Render particles in different areas with different colours.");
DEFINE_int32(num_particles, 2500, "Number of particles to render.");
DEFINE_double(viscosity, 1.308, "Viscosity of the fluid.");
DEFINE_int32(obstacles, 0, "The configuration of obstacles to use.");
DEFINE_double(time_step, 0.005, "Time step. Descrease this to slow down the simulation.");
DEFINE_int32(threshold, 1, "Threshold.");

MAC::MAC(size_t _resolution) : m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 0.0f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 0.0f));
    m_pressure = std::vector<std::vector<double>>(m_resolution, std::vector<double>(m_resolution, 0.0));
    m_density = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution, AIR_DENSITY));

    m_type = std::vector<std::vector<Type>>(m_resolution, std::vector<Type>(m_resolution, FLUID));
    m_particles = std::vector<Position>(FLAGS_num_particles, Position(0.0f, 0.0f));
    m_numParticles = std::vector<std::vector<size_t>>(m_resolution, std::vector<size_t>(m_resolution, 0));
    for (int i = 0; i < m_resolution; ++i)
    {
        for (int j = 0; j < m_resolution; ++j)
        {
            if (j==0 || j==m_resolution-1 || i==0 || i==m_resolution-1)
            {
                m_type[i][j] = SOLID;
            }
        }
    }

    setupObstacles();

    int j = m_resolution;
    for (int i = 1; i < m_resolution-1; ++i)
    {
        m_x[i][j] = 0;
    }

    size_t i = m_resolution;
    for (int j = 1; j < m_resolution-1; ++j)
    {
        m_y[i][j] = 0;
    }

    cellWidth = gridWidth/m_resolution;

    float x_width = gridWidth-4.0f*cellWidth*0.9f;
    float y_height = 0.2f;
    float y_offset = 0.25f;

    for (Position &p: m_particles)
    {
        Index index;
        do
        {
        p.m_x = ((rand() % 1000)/1000.0f - 0.5f) * x_width;
        p.m_y = (((rand() % 1000)/1000.0f - 0.5f) * y_height) + y_offset;
        positionToCellIndex(p,index);
        } while(isSolidCell(index));
        if (FLAGS_colour)
        {
            if (p.m_x < 0)
            {
                m_particleColours.push_back({1.0,0.4,0.0});
            }
            else
            {
                m_particleColours.push_back({0,0.4,1.0});
            }
        }
        else
        {
             m_particleColours.push_back({0,0.4,1.0});
        }
    }

    setupVAO();
    setupVBO();
    setupGridVAO();
    setupGridVBO();
}

void MAC::setupObstacles()
{
    struct obstacle
    {
        Index center;
        size_t size;
    };

    std::vector<obstacle> obstacles;
    size_t size = m_resolution/10;

    if (FLAGS_obstacles >= 1)
    {
        obstacles.push_back({{m_resolution/3,m_resolution/2},size*2});
    }
    if (FLAGS_obstacles >= 2)
    {
        obstacles.push_back({{m_resolution/4,m_resolution/4},size});
    }
    if (FLAGS_obstacles >= 3)
    {
        obstacles.push_back({{m_resolution/4,m_resolution-m_resolution/4},size});
    }

    for (const auto& obstacle: obstacles)
    {
        Index center=obstacle.center;
        size_t size=obstacle.size;
        size_t startRow=center.row-((size)/2);
        size_t endRow=startRow+(size-1);
        size_t startCol=center.col-((size)/2);
        size_t endCol=startCol+(size-1);
        for (size_t i=startRow; i<=endRow; ++i)
        {
            for (size_t j=startCol; j<=endCol; ++j)
            {
                m_type[i][j] = SOLID;
            }
        }
    }
}

void MAC::setupGridVAO()
{
    m_grid_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
}

void MAC::setupGridVBO()
{
    // Set up grid VBO.
    m_grid_vbo.clear();
    // Outer
    float u=0.5;
    float v=-0.5;
    m_grid_vbo.push_back({u,v,0});
    m_grid_vbo.push_back({u,-v,0});
    m_grid_vbo.push_back({-u,v,0});
    m_grid_vbo.push_back({-u,-v,0});

    m_grid_vbo.push_back({-u,v,0});
    m_grid_vbo.push_back({u,v,0});
    m_grid_vbo.push_back({-u,-v,0});
    m_grid_vbo.push_back({u,-v,0});

    // Inner
    u-=cellWidth;
    v+=cellWidth;
    m_grid_vbo.push_back({u,v,0});
    m_grid_vbo.push_back({u,-v,0});
    m_grid_vbo.push_back({-u,v,0});
    m_grid_vbo.push_back({-u,-v,0});

    m_grid_vbo.push_back({-u,v,0});
    m_grid_vbo.push_back({u,v,0});
    m_grid_vbo.push_back({-u,-v,0});
    m_grid_vbo.push_back({u,-v,0});

    // Set up boxes.
    Index index;
    Position p;
    for (index.row=0;index.row<m_resolution;index.row++)
    {
        for(index.col=0;index.col<m_resolution;index.col++)
        {
            if(isSolidCell(index))
            {
                cellIndexToPosition(index,p);
                const float x_pos=p.m_x+0.5f*cellWidth;
                const float y_pos=p.m_y+0.5f*cellWidth;
                const float x_neg=p.m_x-0.5f*cellWidth;
                const float y_neg=p.m_y-0.5f*cellWidth;
                // Vertical
                m_grid_vbo.push_back({x_pos,y_pos,0.0f});
                m_grid_vbo.push_back({x_pos,y_neg,0.0f});
                m_grid_vbo.push_back({x_neg,y_pos,0.0f});
                m_grid_vbo.push_back({x_neg,y_neg,0.0f});
                // Horizontal
                m_grid_vbo.push_back({x_neg,y_pos,0.0f});
                m_grid_vbo.push_back({x_pos,y_pos,0.0f});
                m_grid_vbo.push_back({x_neg,y_neg,0.0f});
                m_grid_vbo.push_back({x_pos,y_neg,0.0f});
            }
        }
    }

    const size_t &size = m_grid_vbo.size();
    m_grid_vao->bind();
    m_grid_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_grid_vbo[0].m_x));
        m_grid_vao->setNumIndices(size);
        m_grid_vao->setVertexAttributePointer(0,3,GL_FLOAT,1*(GLsizei)sizeof(ngl::Vec3),0); // Position.
        m_grid_vao->setMode(GL_LINES);
    m_grid_vao->unbind();
}

MAC& MAC::operator=(MAC&& other)
{
    m_x = other.m_x;
    m_y = other.m_y;
    m_pressure = other.m_pressure;
    m_density = other.m_density;
    m_type = other.m_type;
    m_numParticles = other.m_numParticles;
    m_particles = other.m_particles;
    m_resolution = other.m_resolution;
    gridWidth = other.gridWidth;
    cellWidth = other.cellWidth;
    m_vao = std::move(other.m_vao);
    m_vbo = other.m_vbo;
    m_grid_vao = std::move(other.m_grid_vao);
    m_grid_vbo = other.m_grid_vbo;
    m_particleColours = other.m_particleColours;
    return *this;
}

MAC::MAC(MAC&& other)
{
    m_x = other.m_x;
    m_y = other.m_y;
    m_pressure = other.m_pressure;
    m_density = other.m_density;
    m_type = other.m_type;
    m_numParticles = other.m_numParticles;
    m_particles = other.m_particles;
    m_resolution = other.m_resolution;
    gridWidth = other.gridWidth;
    cellWidth = other.cellWidth;
    m_vao = std::move(other.m_vao);
    m_vbo = other.m_vbo;
    m_grid_vao = std::move(other.m_grid_vao);
    m_grid_vbo = other.m_grid_vbo;
    m_particleColours = other.m_particleColours;
}

void MAC::setupVAO()
{
    m_vao = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
}

void MAC::setupVBO()
{
    updateVBO();
}

void MAC::updateVBO()
{
    ngl::Vec3 backgroundColour = {0.2f, 0.2f, 0.2f};
    m_vbo.clear();
    for (size_t i=0; i<m_particles.size(); ++i)
    {
        const Position &v = m_particles[i];
        m_vbo.push_back({v.m_x, v.m_y, 0.0f});
        Index index;
        positionToCellIndex(v, index);
        if (isInSolidCell(v)) m_vbo.push_back(backgroundColour);
        else m_vbo.push_back(m_particleColours[i]);
    }
    const size_t &size = m_particles.size();
    m_vao->bind();
        m_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3)*2, m_vbo[0].m_x));
        m_vao->setNumIndices(size);
        m_vao->setVertexAttributePointer(0,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),0); // Position.
        m_vao->setVertexAttributePointer(1,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),3); // Colour.
        m_vao->setMode(GL_POINTS);
    m_vao->unbind();
}

void MAC::update()
{
    static size_t time_elapsed = 0;
    const size_t step = 25;
    if (time_elapsed%step == 0)
    {
    }
    updateVectorField();

    // Only update if we hit a frame.
    if (m_frame)
    {
        updateVBO();
        m_frame=false;
    }
    time_elapsed++;
}

void MAC::draw()
{
    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    shader->use("Grid");
    shader->setUniform("colour", ngl::Vec3(0.0,0.0,0.0));
    m_grid_vao->bind();
    m_grid_vao->draw();
    m_grid_vao->unbind();

    shader->use("Particle");
    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
}
