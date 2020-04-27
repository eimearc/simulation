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
constexpr int numWaterParticlesPerPoint = 100;
constexpr int ATMOSPHERIC_PRESSURE = 101325;
constexpr int WATER_DENSITY = 1000;
constexpr float AIR_DENSITY = 1.3f;

MAC::MAC(size_t _resolution) : m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 0.5f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 0.5f));
    m_pressure = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution, 0.0f));
    m_type = std::vector<std::vector<std::string>>(m_resolution, std::vector<std::string>(m_resolution, FLUID));
    m_particles = std::vector<ngl::Vec2>(1000, ngl::Vec2(0.0f, 0.0f));
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
        p.m_x = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f) * ratio * 0.8;
        p.m_y = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f) * ratio * 0.8;
    }

    fixBorderVelocities();

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
    _time = 5; // Needs to be this high.
    static size_t time_elapsed = 0;
    const size_t step = 10;
    if (time_elapsed%step == 0)
    {
        updateVectorField(_time);
        updateVBO();
        std::cout<<*this<<std::endl;
    }

    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
    time_elapsed++;
}

void MAC::updateVectorField(float _time)
{
    _time = 0.5f;
    updateGrid();
    applyConvection(_time);
    applyExternalForces(_time);
//    applyViscosity(_time);
    calculatePressure(_time);
    applyPressure(_time);
    moveParticles(_time);
    fixBorderVelocities();
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
                std::cout << col << std::endl;
                if (m_x[row][col] > maxU) maxU = m_x[row][col];
            }
            if (col < m_resolution)
            {
                if (m_y[row][col] > maxU) maxU = m_y[row][col];
            }
        }
    }
    return 100*(cellWidth/sqrt(maxU*maxU));
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
        size_t row, col;
        positionToCellIndex(p.m_x, p.m_y, row, col);
        if (!outOfBounds(row, col) && !isSolidCell(row, col)) m_type[row][col] = FLUID; // Ensure not boundary.
    }
}

void MAC::cellIndexToPositionX(size_t row, size_t col, float &x, float &y)
{
    cellIndexToPosition(row, col, x, y);
    y += 0.5*cellWidth;
}

void MAC::cellIndexToPositionY(size_t row, size_t col, float &x, float &y)
{
    cellIndexToPosition(row, col, x, y);
    x += 0.5*cellWidth;
}

void MAC::cellIndexToPosition(size_t row, size_t col, float &x, float &y)
{
    float l = -gridWidth/2.0f;
    x = l + col*(cellWidth);
    y = l + row*(cellWidth);
}

void MAC::applyConvection(float _time)
{
    float x = 0.0f, y = 0.0f;
    MAC tmp(m_resolution);
    for (size_t row = 1; row < m_resolution-1; ++row)
    {
        for (size_t col = 1; col < m_resolution; ++col)
        {
            cellIndexToPosition(row, col, x, y);
            ngl::Vec2 updated = traceParticle(x, y, _time);
            tmp.m_x[row][col] = updated.m_x;
        }
    }

    for (size_t row = 1; row < m_resolution; ++row)
    {
        for (size_t col = 1; col < m_resolution-1; ++col)
        {
            cellIndexToPosition(row, col, x, y);
            ngl::Vec2 updated = traceParticle(x, y, _time);
            tmp.m_y[row][col] = updated.m_y;
        }
    }

    tmp.fixBorderVelocities();
    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

void MAC::applyExternalForces(float _time)
{
    ngl::Vec2 gravityVector = {0, -9.80665};
    gravityVector*=_time;
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution+1; ++row)
        {
            m_y[row][col] += gravityVector.m_y;
        }
    }
}

//void applyViscosity(float _time) {}

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
//    Eigen::IncompleteCholesky<Eigen::SparseMatrix<double>> solver;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    std::cout << "A:\n" << A << std::endl;
    std::cout << "b:\n" << b << std::endl;
    Eigen::VectorXd p = solver.solve(b);
    std::cout << "p:\n" << p << std::endl;

    std::cout << "Answer: " << (A*p-b).norm()/b.norm() << "\n";

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

    float x = 0.0f, y = 0.0f;

    for (size_t row = 0; row <= m_resolution; ++row)
    {
        for (size_t col = 0; col <= m_resolution; ++col)
        {
            if (bordersFluidCellX(row, col) && !(bordersSolidCellX(row, col)))
            {
                cellIndexToPositionX(row, col, x, y);
                ngl::Vec2 v = applyPressureToPoint(x,y,_time);
                tmp.m_x[row][col] = v.m_x;
            }

            if (bordersFluidCellY(row, col) && !(bordersSolidCellY(row, col)))
            {
                cellIndexToPositionY(row, col, x, y);
                ngl::Vec2 v = applyPressureToPoint(x,y,_time);
                tmp.m_y[row][col] = v.m_y;
            }
        }
    }

    tmp.fixBorderVelocities();
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
    return ngl::Vec2(x1-x2,y1-y2);
}

ngl::Vec2 MAC::applyPressureToPoint(float x, float y, float _time)
{
    ngl::Vec2 v = velocityAt(x,y);
    size_t row, col;
    positionToCellIndex(x,y,row,col);
    ngl::Vec2 gradient = calculatePressureGradient(row, col);
    float density = 1000.0f; // Water density is 1000kg/m^3.
    if (y>=0)
    {
        density = 1.3f; // Air density is 1.3kg/m^3.
    }

    std::cout << col << "," << row << " Gradient: " << gradient << std::endl;

    auto rhs = (_time/(density*cellWidth))*gradient;
    return (v - rhs);
}

bool MAC::isOutsideGrid(ngl::Vec2 p)
{
    size_t row, col;
    positionToCellIndex(p.m_x, p.m_y, row, col);
    if (row > m_resolution-1 || row < 1 || col > m_resolution-1 || col < 1)
    {
        return true;
    }
    return false;
}

void MAC::moveParticles(float _time)
{
    for (ngl::Vec2 &p : m_particles)
    {
        ngl::Vec2 velocity = velocityAt(p.m_x, p.m_y);
        ngl::Vec2 newPos = p + _time*velocity*0.1;
//        p += velocity;
        if (isOutsideGrid(newPos))
        {
//            newPos = p; // TODO: Change.
        }
        p = newPos;
    }
}


// =====================
// Velocity Methods
// =====================
ngl::Vec2 MAC::velocityAt(const float x, const float y)
{
    ngl::Vec2 v;
    size_t row, col;
    positionToCellIndex(x,y,row,col);

    float x1Pos, x2Pos, y1Pos, y2Pos;
    cellIndexToPosition(row, col, x1Pos, y1Pos);
    cellIndexToPosition(row+1, col+1, x2Pos, y2Pos);

    if (outOfBounds(row, col))
    {
        return ngl::Vec2();
    }

    float x1 = m_x[row][col], x2 = 0.0f, x3 = 0.0f, x4 = 0.0f;
    float y1 = m_y[row][col], y2 = 0.0f, y3 = 0.0f, y4 = 0.0f;
    if (row < m_resolution-1)
    {
        // Top row of the grid.
        x3 = m_x[row+1][col];
        y3 = m_y[row+1][col];
        if (col < m_resolution-1)
        {
            x2 = m_x[row][col+1];
            y2 = m_y[row][col+1];
            x4 = m_x[row+1][col+1];
            y4 = m_y[row+1][col+1];
        }
    }

    v.m_x = (
        (x2Pos-x) * (y2Pos-y) * x1 +
        (x-x1Pos) * (y2Pos-y) * x2 +
        (x2Pos-x) * (y-y1Pos) * x3 +
        (x-x1Pos) * (y-y1Pos) * x4
    );

    v.m_y = (
        (x2Pos-x) * (y2Pos-y) * y1 +
        (x-x1Pos) * (y2Pos-y) * y2 +
        (x2Pos-x) * (y-y1Pos) * y3 +
        (x-x1Pos) * (y-y1Pos) * y4
    );

    return v;
}

ngl::Vec2 MAC::traceParticle(float _x, float _y, float _time)
{
    // Trace particle from point (_x, _y) using simple forward Euler.
    // TODO: update to use RK2.
    ngl::Vec2 v = velocityAt(_x, _y);
    ngl::Vec2 prev_pos = ngl::Vec2(_x, _y) - _time*v;
    ngl::Vec2 new_velocity = velocityAt(prev_pos.m_x, prev_pos.m_y);
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
    // One above bottom row y.
    for (float &v : m_y[1])
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
//        col[m_resolution-1] = 0.0f;
        col[0] = 0.0f;
//        col[1] = 0.0f;
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
                i = index(row, col);
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
    size_t row, col;
    for (const auto &p : m_particles)
    {
        positionToCellIndex(p.m_x, p.m_y, row, col);
        if (!outOfBounds(row, col))
        {
            numParticles[row][col]++;
        }
    }

    m_numParticles = numParticles;

    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
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

ngl::Vec2 MAC::velocityAt(size_t row, size_t col)
{
    float x, y;
    cellIndexToPosition(row, col, x, y);
    return velocityAt(x,y);
}

float MAC::calculateModifiedDivergence(size_t row, size_t col)
{
    // Div(u)
    // Velocity components between fluid cells and solid cells
    // are considered to be zero.
    // TODO:ensure above holds.

    if (outOfBounds(row, col)) return 0.0f;
    if ((row>=m_resolution-1) || (col>=m_resolution-1)) return 0.0f; // At edges.
    if ((row==0) || (col==0)) return 0.0f; // At edges.

    const ngl::Vec2 v = velocityAt(row, col);

    float x1 = v.m_x;
    float x2 = 0;
    if (m_type[row][col+1] == FLUID)
    {
        ngl::Vec2 v2  = velocityAt(row, col+1);
        x2 = v2.m_x;
    }
    float xDiv = x2-x1;
    float y1 = v.m_x;
    float y2 = 0;
    if (m_type[row+1][col] == FLUID)
    {
        ngl::Vec2 v2  = velocityAt(row+1, col);
        y2 = v2.m_y;
    }
    float yDiv = y2-y1;

    return xDiv + yDiv;
}


// =====================
// Helper Methods
// =====================
void MAC::positionToCellIndex(float x, float y, size_t &row, size_t &col)
{
    x += gridWidth/2.0f;
    y += gridWidth/2.0f;
    col = x/cellWidth;
    row = y/cellWidth;
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
        i = index(row+1, col);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get lower neighbour.
    if (row > 0)
    {
        type = getType(row-1,col);
        i = index(row-1, col);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get right neighbour.
    if (col < m_resolution-1)
    {
        type = getType(row,col+1);
        i = index(row, col+1);
        m.insert(std::pair<size_t, std::string>(i,type));
    }
    // Get left neighbour.
    if (col > 0)
    {
        type = getType(row,col-1);
        i = index(row, col-1);
        m.insert(std::pair<size_t, std::string>(i,type));
    }

    return m;
}

size_t MAC::index(size_t row, size_t col)
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
