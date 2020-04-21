#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdlib>
#include <ngl/SimpleVAO.h>
#include <ngl/NGLInit.h>

const std::string FLUID = "FLUID";
const std::string SOLID = "SOLID";
constexpr int numWaterParticlesPerPoint = 10;

MAC::MAC(size_t _resolution) : m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 1.0f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 0.5f));
    m_type = std::vector<std::vector<std::string>>(m_resolution, std::vector<std::string>(m_resolution, FLUID));
    m_particles = std::vector<ngl::Vec2>(m_resolution*m_resolution, ngl::Vec2(0.0f, 0.0f));
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
//        p.m_x = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5*gridWidth)*0.9;
//        p.m_y = (((rand() % (int(gridWidth*100))) / 100.0f) - 0.5*gridWidth)*0.9;
        p.m_x = ((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f;
        p.m_y = ((rand() % (int(gridWidth*100))) / 100.0f) - 0.5f;
    }

    fixBorderVelocities();

    setupVAO();
    setupVBO();
}

MAC& MAC::operator=(MAC&& other)
{
    m_x = other.m_x;
    m_y = other.m_y;
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
    time_elapsed++;
    const size_t step = 20;
    if (time_elapsed%step == 0)
    {
        std::cout << time_elapsed/step << '\n' << *this << std::endl;
        updateVectorField(_time);
        updateVBO();
    }
    m_vao->bind();
    m_vao->draw();
    m_vao->unbind();
}

void MAC::updateVectorField(float _time)
{
    applyConvection(_time);
    applyExternalForces(_time);
//    applyViscosity(_time);
    calculatePressure(_time);
//    applyPressure(_time);
    moveParticles(_time);
    fixBorderVelocities();
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
    const float startOffset = 1.5f;
    MAC tmp(m_resolution);
    for (float row = startOffset; row < m_resolution; row+=1.0f)
    {
        for (size_t col = 1; col < m_resolution-1; ++col)
        {
            ngl::Vec2 updated = traceParticle(row, col, _time);
            tmp.m_x[floor(row)][col] = updated.m_x;
        }
    }

    for (size_t row = 1; row < m_resolution-1; ++row)
    {
        for (float col = startOffset; col < m_resolution; col+=1.0f)
        {
            ngl::Vec2 updated = traceParticle(row, col, _time);
            tmp.m_y[row][floor(col)] = updated.m_y;
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

void applyViscosity(float _time) {}

void MAC::calculatePressure(float _time)
{
    auto A = constructCoefficientMatrix();
    auto b = constructDivergenceVector(_time);
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
    Eigen::VectorXd p = chol.solve(b);
}

void applyPressure(float _time) {}

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
//        ngl::Vec2 newPos = p + _time*velocity;
        p += velocity;
        if (isOutsideGrid(p))
        {
            p -= velocity;
//            newPos = p + _time*velocity*0.001;
        }
//        p = newPos;
    }
}


// =====================
// Velocity Methods
// =====================
ngl::Vec2 MAC::velocityAt(const float x, const float y)
{
    ngl::Vec2 v;
//    const int row = floor(x);
//    const int col = floor(y);
    size_t row, col;
    positionToCellIndex(x,y,row,col);

    if (outOfBounds(row, col))
    {
        std::cout << row << ", " << col << " Point is OOB\n";
        return ngl::Vec2();
    }

    float x1 = m_x[row][col], x2 = 0.0f, x3 = 0.0f, x4 = 0.0f;
    float y1 = m_y[row][col], y2 = 0.0f, y3 = 0.0f, y4 = 0.0f;
    if (row < m_resolution-1)
    {
        // Top row of the grid.
        x3 = m_x[row+1][col];
        y3 = m_y[row+1][col];
        if (col < int(m_resolution-1))
        {
            x2 = m_x[row][col+1];
            y2 = m_y[row][col+1];
            x4 = m_x[row+1][col+1];
            y4 = m_y[row+1][col+1];
        }
    }

    v.m_x = (
        (row+1-x) * (col+1-y) * x1 +
        (x-row) * (col+1-y) * x2 +
        (row+1-x) * (y-col) * x3+
        (x-row) * (y-col) * x4
    );

    v.m_y = (
        (row+1-x) * (col+1-y) * y1 +
        (x-row) * (col+1-y) * y2 +
        (row+1-x) * (y-col) * y3+
        (x-row) * (y-col) * y4
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
    // Right col y.
    for (auto &col : m_y)
    {
        col[m_resolution-1] = 0.0f;
    }
    // Left col y.
    for (auto &col : m_y)
    {
        col[0] = 0.0f;
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
    // Right column x.
    for (auto &col : m_x)
    {
        col[m_resolution] = 0.0f;
    }
    // Left column x.
    for (auto &col : m_x)
    {
        col[0] = 0.0f;
    }
}


// =====================
// Pressure Methods
// =====================
Eigen::SparseMatrix<double> MAC::constructCoefficientMatrix()
{
    size_t n = m_resolution;
    Eigen::SparseMatrix<double> m(n*n,n*n);

    auto tripletList = constructNeighbourTriplets();

    m.setFromTriplets(tripletList.begin(), tripletList.end());
    return m;
}

std::vector<Eigen::Triplet<double>> MAC::constructNeighbourTriplets()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::Triplet<double> t;
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
        {
            if (m_type[row][col] == FLUID)
            {
                auto i = index(row, col);
                auto m = getNeighbours(row, col);
                size_t nonSolidNeighbours = 0;
                for ( const auto &e : m)
                {
                    size_t r, c;
                    coordinate(e.first, r, c);
                    auto neighbourIndex = index(r,c);
                    t = Eigen::Triplet<double>(i, neighbourIndex, e.second);
                    tripletList.push_back(t);
                    nonSolidNeighbours+=e.second;
                }
                size_t num = getNumNonLiquidNeighbours(row, col);
                t = Eigen::Triplet<double>(i, i, num);
                tripletList.push_back(t);
            }
        }
    }
    return tripletList;
}

Eigen::VectorXd MAC::constructDivergenceVector(float _time)
{
    Eigen::VectorXd v(m_resolution*m_resolution);

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

    // For each fluid cell...
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
        {
            if (isFluidCell(row, col))
            {
                size_t i = index(row, col);
                double cellArea = cellWidth*cellWidth;
                double density = 0.0f;
                if (numParticles[row][col] > 0)
                {
                    density = numWaterParticlesPerPoint*numParticles[row][col] / cellArea;
                }
                float h = cellWidth;
                float divergence = calculateModifiedDivergence(row,col);
                size_t numNeighbourAirCells = 0;
                int atmosphericPressure = 101325;

                auto result = ((density*h)/_time)*divergence - numNeighbourAirCells*atmosphericPressure;

                v[i] = result;
            }
        }
    }

    return v;
}

float MAC::calculateModifiedDivergence(size_t row, size_t col)
{
    // Div(u)
    // Velocity components between fluid cells and solid cells
    // are considered to be zero.

    if (outOfBounds(row, col)) return 0.0f;
    if ((row==m_resolution-1) || (col==m_resolution-1)) return 0.0f; // At edges.

    float xDiv = 0.0f;
    if (m_type[row+1][col] == FLUID)
    {
        xDiv = m_x[row+1][col] - m_x[row][col];
    }
    float yDiv = 0.0f;
    if (m_type[row][col+1] == FLUID)
    {
        yDiv = m_x[row][col+1] - m_x[row][col];
    }

    return xDiv + yDiv;
}


// =====================
// Helper Methods
// =====================
void MAC::positionToCellIndex(float x, float y, size_t &row, size_t &col)
{
    std::cout << x << "," << y;
    x += gridWidth/2.0f;
    y += gridWidth/2.0f;
    col = x/cellWidth;
    row = y/cellWidth;
    std::cout << " --> " << row << "," << col << " gridWidth: " << gridWidth << " cellWidth: " << cellWidth << std::endl;
}

size_t MAC::getType(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return 0; // Return air (non-solid)
    }
    if (m_type[row][col] == SOLID)
    {
        return 0;
    }
    return 1;
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

std::vector<std::pair<size_t, size_t>> MAC::getNeighbourIndices(size_t row, size_t col)
{
    std::vector<std::pair<size_t, size_t>> indices;

    if (row < m_resolution-1) indices.push_back({row+1, col});
    if (row > 0) indices.push_back({row-1, col});
    if (col < m_resolution-1) indices.push_back({row, col+1});
    if (col > 0) indices.push_back({row, col-1});

    return indices;
}

std::map<size_t, size_t> MAC::getNeighbours(size_t row, size_t col)
{
    std::map<size_t, size_t> m;
    size_t i = 0;
    size_t type = 0;
    // Get upper neighbour.
    if (row < m_resolution-1)
    {
        type = getType(row+1,col);
        i = index(row+1, col);
        m.insert(std::pair<size_t, size_t>(i,type));
    }
    // Get lower neighbour.
    if (row > 0)
    {
        type = getType(row-1,col);
        i = index(row-1, col);
        m.insert(std::pair<size_t, size_t>(i,type));
    }
    // Get right neighbour.
    if (col < m_resolution-1)
    {
        type = getType(row,col+1);
        i = index(row, col+1);
        m.insert(std::pair<size_t, size_t>(i,type));
    }
    // Get left neighbour.
    if (col > 0)
    {
        type = getType(row,col-1);
        i = index(row, col-1);
        m.insert(std::pair<size_t, size_t>(i,type));
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

    std::cout << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int i = mac.m_y.size()-1; i >= 0 ; --i)
    {
        if (i < int(mac.m_x.size()))
        {
            std::cout << "X" << i << "  ";
            for (const auto &x : mac.m_x[i])
            {
                std::cout << x;
                std::cout << "     ";
            }
        }
        std::cout << "\n\n";
        std::cout << "Y" << i << "  ";
        for (const auto &y : mac.m_y[i])
        {
            std::cout <<"     ";
            std::cout << y;
        }

        std::cout << "\n\n";
    }

    return os;
}

