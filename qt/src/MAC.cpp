#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdlib>
#include <ngl/SimpleVAO.h>
#include <ngl/NGLInit.h>

const std::string FLUID = "FLUID";
const std::string SOLID = "SOLID";
const std::string AIR = "AIR";
constexpr float ATMOSPHERIC_PRESSURE = 101325.0f;
//constexpr float WATER_DENSITY = 1000.0f;
constexpr float WATER_DENSITY = 1.0f; // According to notes, water density is always 1.
constexpr float AIR_DENSITY = 1.0f;
constexpr size_t MAX_PARTICLES_PER_CELL = 100;

MAC::MAC(size_t _resolution) : m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 0.0f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 0.0f));
    m_pressure = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution, 0.0f));
    m_type = std::vector<std::vector<std::string>>(m_resolution, std::vector<std::string>(m_resolution, FLUID));
    m_particles = std::vector<ngl::Vec2>(2000, ngl::Vec2(0.0f, 0.0f));
    m_numParticles = std::vector<std::vector<size_t>>(m_resolution, std::vector<size_t>(m_resolution, 0));
    for (size_t i = 0; i < m_resolution; ++i)
    {
        for (size_t j = 0; j < m_resolution; ++j)
        {
            if (j==0 || j==m_resolution-1 || i==0 || i==m_resolution-1)
            {
                m_type[i][j] = SOLID;
            }
        }
    }

    size_t j = m_resolution;
    for (size_t i = 1; i < m_resolution-1; ++i)
    {
        m_x[i][j] = 0;
    }

    size_t i = m_resolution;
    for (size_t j = 1; j < m_resolution-1; ++j)
    {
        m_y[i][j] = 0;
    }

    cellWidth = gridWidth/m_resolution;
    for (ngl::Vec2 &p: m_particles)
    {
        float ratio = (m_resolution-2)/float(m_resolution);
        p.m_x = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f) * ratio * 0.5;
        p.m_y = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f) * ratio * 0.5;
    }

    setupVAO();
    setupVBO();
}

MAC& MAC::operator=(MAC&& other)
{
    m_x = other.m_x;
    m_y = other.m_y;
    m_pressure = other.m_pressure;
    m_type = other.m_type;
    m_numParticles = other.m_numParticles;
    m_particles = other.m_particles;
    m_resolution = other.m_resolution;
    gridWidth = other.gridWidth;
    cellWidth = other.cellWidth;
    m_vao = std::move(other.m_vao);
    m_vbo = other.m_vbo;
    return *this;
}

MAC::MAC(MAC&& other)
{
    m_x = other.m_x;
    m_y = other.m_y;
    m_pressure = other.m_pressure;
    m_type = other.m_type;
    m_numParticles = other.m_numParticles;
    m_particles = other.m_particles;
    m_resolution = other.m_resolution;
    gridWidth = other.gridWidth;
    cellWidth = other.cellWidth;
    m_vao = std::move(other.m_vao);
    m_vbo = other.m_vbo;
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
    m_vbo = m_particles;
    const size_t &size = m_vbo.size();
    m_vao->bind();
        m_vao->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec2), m_vbo[0].m_x));
        m_vao->setNumIndices(size);
        m_vao->setVertexAttributePointer(0,2,GL_FLOAT,1*(GLsizei)sizeof(ngl::Vec2),0); // Position.
        m_vao->setMode(GL_POINTS);
    m_vao->unbind();
}

void MAC::draw(float _time)
{
    static size_t time_elapsed = 0;
    const size_t step = 10;
    if (time_elapsed%step == 0)
    {
        updateVectorField(_time);
        updateVBO();
        std::cout<<"\n*****\nNewRound\n\n" << *this<<std::endl;
    }

    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
    time_elapsed++;
}

void MAC::updateVectorField(float _time)
{
    _time = calculateTimeStep();
    _time *= 0.5;
    updateGrid();
    std::cout << "Grid before convection:\n" << *this << std::endl;
    applyConvection(_time);
    std::cout << "Grid after convection:\n" << *this << std::endl;
    applyExternalForces(_time);
    std::cout << "Grid after external forces:\n" << *this << std::endl;
//    applyViscosity(_time);
    calculatePressure(_time);
    applyPressure(_time);
    fixBorderVelocities();
    moveParticles(_time);
}

float MAC::calculateTimeStep()
{
    float maxU = 0.0f;
    for (size_t row = 0; row <= m_resolution; ++row)
    {
        for (size_t col = 0; col <= m_resolution; ++col)
        {
            if (row < m_resolution)
            {
                if ((m_x[row][col])*(m_x[row][col]) > maxU) maxU = (m_x[row][col])*(m_x[row][col]);
            }
            if (col < m_resolution)
            {
                if ((m_y[row][col])*(m_y[row][col]) > maxU) maxU = (m_y[row][col])*(m_y[row][col]);
            }
        }
    }
    std::cout << "maxU " << sqrt(maxU) << std::endl;
    maxU += 0.000001f;
    return 2*(cellWidth/sqrt(maxU));
}

void MAC::updateGrid()
{
    for (size_t row = 1; row < m_resolution-1 ; ++row)
    {
        for (size_t col = 1; col < m_resolution-1 ; ++col)
        {
            if (!isSolidCell(row, col)) m_type[row][col] = AIR;
        }
    }

    for (const auto &p: m_particles)
    {
        Index index;
        positionToCellIndex(Position(p.m_x, p.m_y), index);
        if (!outOfBounds(index.row, index.col) && !isSolidCell(index.row, index.col)) m_type[index.row][index.col] = FLUID; // Ensure not boundary.
    }
}

void MAC::cellIndexToPositionX(size_t row, size_t col, float &x, float &y)
{
    cellIndexToPosition(row, col, x, y);
    x -= 0.5*cellWidth;
}

void MAC::cellIndexToPositionY(size_t row, size_t col, float &x, float &y)
{
    cellIndexToPosition(row, col, x, y);
    y -= 0.5*cellWidth;
}

void MAC::cellIndexToPosition(size_t row, size_t col, float &x, float &y)
{
    float l = -gridWidth/2.0f;
    x = l + col*(cellWidth) + 0.5*(cellWidth);
    y = l + row*(cellWidth) + 0.5*(cellWidth);
}

void MAC::applyConvection(float _time)
{
    float x = 0.0f, y = 0.0f;
    ngl::Vec2 updated;
    MAC tmp(m_resolution);
    for (size_t row = 0; row <= m_resolution; ++row)
    {
        for (size_t col = 0; col <= m_resolution; ++col)
        {
            if (bordersFluidCellX(row,col) && !bordersSolidCellX(row,col))
            {
                cellIndexToPositionX(row, col, x, y);
                std::cout << " X  row: " << row << ", col: " << col;
                updated = traceParticle(x, y, _time);
                tmp.m_x[row][col] = updated.m_x;
            }
            if(bordersFluidCellY(row,col) && !bordersSolidCellY(row,col))
            {
                cellIndexToPositionY(row, col, x, y);
                std::cout << " Y  row: " << row << ", col: " << col << " current velocity:" << m_y[row][col];
                updated = traceParticle(x, y, _time);
                tmp.m_y[row][col] = updated.m_y;
            }
        }
    }

    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

void MAC::applyExternalForces(float _time)
{
    ngl::Vec2 gravityVector = {0, -9.80665};
    gravityVector*=_time;
    for (size_t col = 0; col <= m_resolution; ++col)
    {
        for (size_t row = 0; row <= m_resolution; ++row)
        {
            if (!outOfBounds(row,col) && (!bordersSolidCellY(row,col)) && (bordersFluidCellY(row,col)))
            {
                m_y[row][col] += gravityVector.m_y;
            }
            if (!outOfBounds(row,col) && (!bordersSolidCellX(row,col)) && (bordersFluidCellY(row,col)))
            {
                m_x[row][col] += gravityVector.m_x;
            }
        }
    }
}

void MAC::calculatePressure(float _time)
{
    m_indices = std::vector<std::vector<int>>(m_resolution, std::vector<int>(m_resolution, -1));
    size_t index = 0;
    for (size_t row = 0; row < m_resolution ; ++row)
    {
        for (size_t col = 0; col < m_resolution; ++col)
        {
            if (isFluidCell(row, col))
            {
                m_indices[row][col] = index;
                ++index;
            }
        }
    }

    auto A = constructCoefficientMatrix();
    auto b = constructDivergenceVector(_time);

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd p = solver.solve(b);

    m_pressure = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution, 0.0f));
    for (size_t row = 0; row < m_resolution; ++row)
    {
        for (size_t col = 0; col < m_resolution; ++col)
        {
            if (isFluidCell(row, col))
            {
                size_t j = m_indices[row][col];
                m_pressure[row][col] = p[j];
            }
        }
    }
}

void MAC::applyPressure(float _time)
{
    MAC tmp(m_resolution);
    tmp.m_x = m_x; // Need this, otherwise particles float.
    tmp.m_y = m_y;

    // Fix atmospheric pressure for air cells.
    for (size_t row = 0; row < m_resolution; ++row)
    {
        for (size_t col = 0; col < m_resolution; ++col)
        {
            if (m_type[row][col] == AIR)
            {
                m_pressure[row][col] = ATMOSPHERIC_PRESSURE;
            }
        }
    }

    // Fix atmospheric pressure for solid cells.
    for (size_t row = 0; row < m_resolution; ++row)
    {
        m_pressure[row][0] = m_pressure[row][1];
        m_pressure[row][m_resolution-1] = m_pressure[row][m_resolution-2];
    }
    for (size_t col = 0; col < m_resolution; ++col)
    {
        m_pressure[0][col] = m_pressure[1][col];
        m_pressure[m_resolution-1][col] = m_pressure[m_resolution-2][col];
    }

    std::cout << "m_pressure\n" << m_pressure << std::endl;

    float x = 0.0f, y = 0.0f;

    for (size_t row = 0; row <= m_resolution; ++row)
    {
        for (size_t col = 0; col <= m_resolution; ++col)
        {
            if (bordersFluidCellX(row, col) && !(bordersSolidCellX(row, col)))
            {
                cellIndexToPositionX(row, col, x, y);
                std::cout << "\tX  Applying to " << row << "," << col << " --> " << x  << "," << y << std::endl;
//                ngl::Vec2 v = applyPressureToPoint(x,y,_time);
                ngl::Vec2 v = applyPressureToPointX(row,col,_time);
                tmp.m_x[row][col] = v.m_x;
                std::cout << "\t\t" << m_x[row][col] << " ---> " << v.m_x << std::endl;
            }

            if (bordersFluidCellY(row, col) && !(bordersSolidCellY(row, col)))
            {
                cellIndexToPositionY(row, col, x, y);
                std::cout << "\tY  Applying to " << row << "," << col << " --> " << x  << "," << y << std::endl;
//                ngl::Vec2 v = applyPressureToPoint(x,y,_time);
                ngl::Vec2 v = applyPressureToPointY(row,col,_time);
                tmp.m_y[row][col] = v.m_y;
                std::cout << "\t\t" << m_y[row][col] << " ---> " << v.m_y << std::endl;
            }
        }
    }

    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

ngl::Vec2 MAC::calculatePressureGradient(size_t row, size_t col)
{
    // Air has a pressure of atmospheric pressure.
    // Air has a density of 1.
    // {p(x,y) - p(x-1,y), p(x,y) - p(x,y-1)}
    float x1=0.f, x2=0.f, y1=0.f, y2=0.f;
    x1 = m_pressure[row][col];
    y1 = m_pressure[row][col];
    if (col > 0) x2 = m_pressure[row][col-1];
    if (row > 0) y2 = m_pressure[row-1][col];
    std::cout << "\t\t x1:" << x1 << " y1:" << y1 << " x2:" << x2 << " y2:" << y2 << std::endl;
    return ngl::Vec2(x1-x2,y1-y2);
}

ngl::Vec2 MAC::applyPressureToPoint(float x, float y, float _time)
{
    Position p(x,y);
    Index index;
    ngl::Vec2 v = velocityAtPosition(p);
    positionToCellIndex(p,index);
    ngl::Vec2 gradient = calculatePressureGradient(index.row, index.col);
    float density = WATER_DENSITY;
//    if (isAirCell(row,col)) density = AIR_DENSITY;

    auto rhs = (_time/(density*cellWidth))*gradient;
    ngl::Vec2 result = v-rhs;
    std::cout << "\t\t" << index.row << "," << index.col;
    std::cout << " Velocity: " << v << " Gradient: " << gradient << " Density:" <<density;
    std::cout << " rhs:" << rhs << std::endl;
    return result;
}

ngl::Vec2 MAC::applyPressureToPointY(const int row, const int col, float _time)
{
    ngl::Vec2 v = {0,m_y[row][col]};
    float p1 = (m_pressure[row][col] - m_pressure[row-1][col])/2.0f;
    float p2 = (m_pressure[row-1][col] - m_pressure[row-2][col])/2.0f;
    ngl::Vec2 gradient = {0,p1-p2};
    float density = WATER_DENSITY;
    if (isAirCell(row,col)) density = AIR_DENSITY;

    auto rhs = (_time/(density*cellWidth))*gradient;
    ngl::Vec2 result = v-rhs;
    std::cout << "\t\t" << row << "," << col;
    std::cout << " Velocity: " << v << " Gradient: " << gradient << " Density:" <<density;
    std::cout << " rhs:" << rhs << std::endl;
    return result;
}

ngl::Vec2 MAC::applyPressureToPointX(const int row, const int col, float _time)
{
    ngl::Vec2 v = {m_x[row][col],0};
    float p1 = (m_pressure[row][col] - m_pressure[row][col-1])/2.0f;
    float p2 = (m_pressure[row][col-1] - m_pressure[row][col-2])/2.0f;
    ngl::Vec2 gradient = {p1-p2,0};
    float density = WATER_DENSITY;
    if (isAirCell(row,col)) density = AIR_DENSITY;

    auto rhs = (_time/(density*cellWidth))*gradient;
    ngl::Vec2 result = v-rhs;
    std::cout << "\t\t" << row << "," << col;
    std::cout << " Velocity: " << v << " Gradient: " << gradient << " Density:" <<density;
    std::cout << " rhs:" << rhs << std::endl;
    return result;
}

bool MAC::isOutsideGrid(const Position &p)
{
    Index index;
    positionToCellIndex(p, index);
    const size_t &row = size_t(index.row);
    const size_t &col = size_t(index.col);
    if (row > m_resolution-1 || row < 1 || col > m_resolution-1 || col < 1)
    {
        return true;
    }
    return false;
}

bool MAC::isOutsideFluid(const Position &p)
{
    if (isOutsideGrid(p)) return true;
    Index index;
    positionToCellIndex(p, index);
    if (isSolidCell(index.row, index.col)) return true;
    return false;
}

void MAC::moveParticles(float _time)
{
    for (ngl::Vec2 &p : m_particles)
    {
        ngl::Vec2 velocity = velocityAtPosition(Position(p.m_x, p.m_y));
        ngl::Vec2 newPos = p + _time*velocity;
        if (isOutsideFluid(newPos))
        {
            newPos = p; // TODO: Change.
        }
        p = newPos;
    }
}


float distance(float x, float y)
{
    return sqrt(x*x)-sqrt(y*y);
}

float MAC::interpolate(const std::vector<std::vector<float>> &m, const float x, const float y, const ngl::Vec2 center, std::string type)
{
    float result = 0.0f;
    Index index;
    Position p(x,y);
    positionToCellIndex(p,index);
    size_t row=index.row;
    size_t col=index.col;

    float x1, x2, y1, y2;
    int tmpRow = row, tmpCol = col;

    float q1=0.0f, q2 = 0.0f, q3 = 0.0f, q4 = 0.0f;

    if (y<center.m_y)
    {
        tmpRow--;
    }
    if (x<center.m_x)
    {
        tmpCol--;
    }

    if(tmpRow<0&&tmpCol<0)
    {
        return m[0][0];
    }
    else if(tmpRow<0)
    {
         return m[0][tmpCol];
    }
    else if(tmpCol<0)
    {
         return m[tmpRow][0];
    }

    if(tmpRow>int(m.size()-1)&&tmpCol>int(m[0].size()-1))
    {
        return m[m.size()-1][m[0].size()-1];
    }
    else if(tmpRow>int(m.size()-1))
    {
        return m[m.size()-1][tmpCol];
    }
    else if(tmpCol>int(m[0].size()-1))
    {
        return m[tmpRow][m.size()-1];
    }

    if (type=="x")
    {
        cellIndexToPositionX(tmpRow, tmpCol, x1, y1);
        cellIndexToPositionX(tmpRow+1, tmpCol+1, x2, y2);
    }
    else
    {
        cellIndexToPositionY(tmpRow, tmpCol, x1, y1);
        cellIndexToPositionY(tmpRow+1, tmpCol+1, x2, y2);
    }
    q1 = m[tmpRow][tmpCol];
    if (tmpCol < int(m[0].size()-1)) q2 = m[tmpRow][tmpCol+1];
    if (tmpRow < int(m.size()-1)) q3 = m[tmpRow+1][tmpCol];
    if (tmpCol < int(m[0].size()-1) && tmpRow < int(m.size()-1)) q4 = m[tmpRow+1][tmpCol+1];

    // X direction.
    auto fx1 = (x2-x)/(x2-x1)*q1 + (x-x1)/(x2-x1)*q2;
    auto fx2 = (x2-x)/(x2-x1)*q3 + (x-x1)/(x2-x1)*q4;
    // Y direction.
    result = (y2-y)/(y2-y1)*fx1 + (y-y1)/(y2-y1)*fx2;
    return result;
}

// =====================
// Velocity Methods
// =====================
ngl::Vec2 MAC::velocityAtPosition(const Position p)
{
    // Separately bilinearly interpolate x and y.
    ngl::Vec2 v;

    ngl::Vec2 center;
    float centerX, centerY;
    Index index;
    positionToCellIndex(p,index);

    cellIndexToPositionX(index.row,index.col,centerX,centerY);
    center.m_x = centerX;
    center.m_y = centerY;
    v.m_x = interpolate(m_x,p.m_x,p.m_y,center, "x");

    cellIndexToPositionY(index.row,index.col,centerX,centerY);
    center.m_x = centerX;
    center.m_y = centerY;
    v.m_y = interpolate(m_y,p.m_x,p.m_y,center, "y");

    return v;
}

ngl::Vec2 MAC::traceParticle(float _x, float _y, float _time)
{
    // Trace particle from point (_x, _y) using RK2.
    Position p(_x,_y);
    ngl::Vec2 v = velocityAtPosition(p);
//    size_t row, col;
    Index index;
    positionToCellIndex(p, index);
    std::cout << "\n\tcurrent pos:" << ngl::Vec2(_x,_y) << " current velocity:" << v << " current row,col:"<< index.row << "," << index.col;
    ngl::Vec2 halfV = v*_time*0.5;
    ngl::Vec2 half_prev_pos = {_x-halfV.m_x, _y-halfV.m_y};
    std::cout << "\n\thalf prev pos:" << half_prev_pos << std::endl;
    ngl::Vec2 half_prev_v = velocityAtPosition(half_prev_pos);
    std::cout << "\thalf prev v: " << half_prev_v;
    ngl::Vec2 full_prev_pos = ngl::Vec2(_x, _y) - _time*half_prev_v;
    std::cout << "\n\tfull prev pos:" << full_prev_pos;
    positionToCellIndex(full_prev_pos, index);
    std::cout << "\n\t\tprev row:" << index.row << ", col:" << index.col;
    ngl::Vec2 new_velocity = velocityAtPosition(full_prev_pos);
    std::cout << "\n\tprev velocity: " << new_velocity << std::endl;

    return new_velocity;
}

void MAC::fixBorderVelocities()
{
    // Top row y.
    for (float &v : m_y[m_resolution])
    {
        v = 0.0f;
    }
    // Bottom row y.
    for (float &v : m_y[0])
    {
        v = 0.0f;
    }
    // Right and left col y.
    for (auto &col : m_y)
    {
        col[0] = 0.0f;
        col[m_resolution-1] = 0.0f;
    }
    // Top row x.
    for (auto &v : m_x[m_resolution-1])
    {
        v = 0.0f;
    }
    // Bottom row x.
    for (auto &v : m_x[0])
    {
        v = 0.0f;
    }
    // Left and right solids x.
    for (auto &col : m_x)
    {
        col[m_resolution] = 0.0f;
        col[0] = 0.0f;
    }
}


// =====================
// Pressure Methods
// =====================
Eigen::SparseMatrix<double> MAC::constructCoefficientMatrix()
{
    size_t n = numFluidCells();
    Eigen::SparseMatrix<double> m(n,n);
    auto tripletList = constructNeighbourTriplets();
    m.setFromTriplets(tripletList.begin(), tripletList.end());
    return m;
}

size_t typeToIndex(std::string type)
{
    if (type==FLUID)
    {
        return 1;
    }
    return 0;
}

std::vector<Eigen::Triplet<double>> MAC::constructNeighbourTriplets()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::Triplet<double> t;
    size_t i;
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
        {
            if(isFluidCell(row, col))
            {
                i = vectorIndex(row, col);
                i = m_indices[row][col];
                t = Eigen::Triplet<double>(i,i,-1*int(getNumNonSolidNeighbours(row, col)));
                tripletList.push_back(t);

                auto neighbours = getNeighbourIndices(row, col);
                for ( const auto &neighbour : neighbours)
                {
                    const size_t &row = neighbour.first;
                    const size_t &col = neighbour.second;
                    if(isFluidCell(row, col))
                    {
                        size_t j = m_indices[row][col];
                        t = Eigen::Triplet<double>(i,j,1);
                        tripletList.push_back(t);
                    }
                }
            }
        }
    }
    return tripletList;
}

Eigen::VectorXd MAC::constructDivergenceVector(float _time)
{
    size_t n = numFluidCells();

    Eigen::VectorXd v(n);
    v.setZero();
    std::vector<std::vector<size_t>> numParticles(m_resolution, std::vector<size_t>(m_resolution, 0));
    Index index;
    for (const auto &p : m_particles)
    {
        positionToCellIndex(Position(p.m_x, p.m_y), index);
        if (!outOfBounds(index.row, index.col))
        {
            numParticles[index.row][index.col]++;
        }
    }

    m_numParticles = numParticles;
    std::cout << "Number of particles:\n" << m_numParticles << std::endl;

    for (size_t col = 0; col <= m_resolution; ++col)
    {
        for (size_t row = 0; row <= m_resolution; ++row)
        {
            if (isFluidCell(row, col))
            {
                size_t i = m_indices[row][col];
                float density = WATER_DENSITY;
                float h = cellWidth;
                float divergence = calculateModifiedDivergence(row,col);

                size_t numNeighbourAirCells = 0;
                if (isAirCell(row,col-1)) numNeighbourAirCells++;
                if (isAirCell(row,col+1)) numNeighbourAirCells++;
                if (isAirCell(row-1,col)) numNeighbourAirCells++;
                if (isAirCell(row+1,col)) numNeighbourAirCells++;

                auto result = ((density*h)/_time)*divergence - numNeighbourAirCells*ATMOSPHERIC_PRESSURE;

                v[i] = result;
            }
        }
    }

    return v;
}

ngl::Vec2 MAC::velocityAtIndex(const Index index)
{
    float x, y;
    size_t row = index.row;
    size_t col = index.col;
    cellIndexToPosition(row, col, x, y);
    Position p(x,y);
    return velocityAtPosition(p);
}

float MAC::calculateModifiedDivergence(size_t row, size_t col)
{
    // Div(u)
    // Velocity components between fluid cells and solid cells
    // are considered to be zero.
    Index index{int(row),int(col)};
    const ngl::Vec2 v = velocityAtIndex(index);

    float x1 = v.m_x;
    float x2 = 0;
    if (!isSolidCell(row,col+1)) // U between fluid and solid considered to be zero.
    {
        index = Index{int(row),int(index.col+1)};
        ngl::Vec2 v2  = velocityAtIndex(index);
        x2 = v2.m_x;
    }
    float xDiv = x2-x1;

    float y1 = v.m_x;
    float y2 = 0;
    if (!isSolidCell(row+1,col))
    {
        index = Index{int(row+1),int(col)};
        ngl::Vec2 v2  = velocityAtIndex(index);
        y2 = v2.m_y;
    }
    float yDiv = y2-y1;

    return xDiv + yDiv;
}


// =====================
// Helper Methods
// =====================
void MAC::positionToCellIndex(const Position &position, Index &index)
{
    float x = position.m_x;
    float y = position.m_y;
    x += gridWidth/2.0f;
    y += gridWidth/2.0f;
    index.col = x/cellWidth;
    index.row = y/cellWidth;
}

std::string MAC::getType(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return AIR; // Return air (non-solid)
    }
    return m_type[row][col];
}

size_t MAC::getNumNonLiquidNeighbours(size_t row, size_t col)
{
    auto neighbours = getNeighbourIndices(row, col);
    size_t num = 0;
    for (const auto& pair : neighbours)
    {
        if (!isFluidCell(pair.first, pair.second))
        {
            num++;
        }
    }
    return num;
}

size_t MAC::getNumNonSolidNeighbours(size_t row, size_t col)
{
    auto neighbours = getNeighbourIndices(row, col);
    size_t num = 0;
    for (const auto& pair : neighbours)
    {
        if (!isSolidCell(pair.first, pair.second))
        {
            num++;
        }
    }
    return num;
}

std::vector<std::pair<size_t, size_t>> MAC::getNeighbourIndices(size_t row, size_t col)
{
    std::vector<std::pair<size_t, size_t>> indices;

    if (row < m_resolution-1) indices.push_back({row+1, col});
    if (row > 0) indices.push_back({row-1, col});
    if (col < m_resolution-1) indices.push_back({row, col+1});
    if (col > 0) indices.push_back({row, col-1});

    return indices;
}

std::map<size_t, std::string> MAC::getNeighbourType(size_t row, size_t col)
{
    std::map<size_t, std::string> m;
    size_t i = 0;
    std::string type;
    // Get upper neighbour.
    if (row < m_resolution-1)
    {
        type = getType(row+1,col);
        i = vectorIndex(row+1, col);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get lower neighbour.
    if (row > 0)
    {
        type = getType(row-1,col);
        i = vectorIndex(row-1, col);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get right neighbour.
    if (col < m_resolution-1)
    {
        type = getType(row,col+1);
        i = vectorIndex(row, col+1);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get left neighbour.
    if (col > 0)
    {
        type = getType(row,col-1);
        i = vectorIndex(row, col-1);
        m.insert(std::pair<size_t, std::string>(i,type));
    }

    return m;
}

size_t MAC::vectorIndex(size_t row, size_t col)
{
    return row*m_resolution + col;
}

void MAC::coordinate(size_t index, size_t &row, size_t &col)
{
    col = index%m_resolution;
    row = index/m_resolution;
}

bool MAC::outOfBounds(size_t row, size_t col)
{
    return ((row >= m_resolution) || (col >= m_resolution) || (row < 0) || (col < 0));
}

bool MAC::isFluidCell(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return false;
    }
    return m_type[row][col] == FLUID;
}

bool MAC::isSolidCell(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return false;
    }
    return m_type[row][col] == SOLID;
}

bool MAC::isAirCell(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return false;
    }
    return m_type[row][col] == AIR;
}

size_t MAC::numFluidCells()
{
    size_t num = 0;
    for (size_t row = 0; row <m_resolution; ++row)
    {
        for (size_t col = 0; col < m_resolution; ++col)
        {
            if (isFluidCell(row, col))
            {
                num++;
            }
        }
    }
    return num;
}

bool MAC::bordersSolidCellX(size_t row, size_t col)
{
    if (outOfBounds(row, col)) return false;
    if (isSolidCell(row, col)) return true;
    if(outOfBounds(row,col-1)) return false;
    if (isSolidCell(row, col-1)) return true;
    return false;
}

bool MAC::bordersSolidCellY(size_t row, size_t col)
{
    if (outOfBounds(row, col)) return false;
    if (isSolidCell(row, col)) return true;
    if(outOfBounds(row-1,col)) return false;
    if (isSolidCell(row-1, col)) return true;
    return false;
}

bool MAC::bordersFluidCellX(size_t row, size_t col)
{
    if (outOfBounds(row, col)) return false;
    if (isFluidCell(row, col)) return true;
    if(outOfBounds(row,col-1)) return false;
    if (isFluidCell(row, col-1)) return true;
    return false;
}

bool MAC::bordersFluidCellY(size_t row, size_t col)
{
    if (outOfBounds(row, col)) return false;
    if (isFluidCell(row, col)) return true;
    if(outOfBounds(row-1,col)) return false;
    if (isFluidCell(row-1, col)) return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, MAC& mac)
{
    for (int i = mac.m_type.size()-1; i >= 0; --i)
    {
        os << i << ":\t";
        for (const auto &y: mac.m_type[i])
        {
            os << y << ' ';
        }
        os << '\n';
    }
    os << '\n';

    os << std::fixed << std::setprecision(4) << std::setfill('0') << std::showpos;
    for (int i = mac.m_y.size()-1; i >= 0 ; --i)
    {
        if (i < int(mac.m_x.size()))
        {
            os << "X" << i << "  ";
            for (const auto &x : mac.m_x[i])
            {
                os << x;
                os << "     ";
            }
        }
        os << "\n\n";
        os << "Y" << i << "  ";
        for (const auto &y : mac.m_y[i])
        {
            os <<"  *  ";
            os << y;
        }
        os <<"  *  ";

        os << "\n\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<float>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<size_t>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<int>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}
