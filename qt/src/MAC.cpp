#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdlib>

MAC::MAC(size_t _resolution) :
    m_pressure(_resolution), m_velocityX(_resolution),
    m_velocityY(_resolution), m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 1.0f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 1.0f));
    m_type = std::vector<std::vector<std::string>>(m_resolution, std::vector<std::string>(m_resolution, "fluid"));
    m_particles = std::vector<ngl::Vec2>(m_resolution, ngl::Vec2(0.0f, 0.0f));
    for (size_t i = 0; i < m_resolution; ++i)
    {
        for (size_t j = 0; j < m_resolution; ++j)
        {
            if (j==0 || j==m_resolution-1 || i==0 || i==m_resolution-1)
            {
                m_type[i][j] = "solid";
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

    for (ngl::Vec2 &p: m_particles)
    {
        p.m_x = (rand() % (m_resolution-2))+1;
        p.m_y = (rand() % (m_resolution-2))+1;
    }

    fixBorderVelocities();
}

void MAC::calculatePressure(float _time)
{
    auto A = constructCoefficientMatrix();
    auto b = constructDivergenceVector(_time);
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(A);
    Eigen::VectorXd p = chol.solve(b);
}

Eigen::SparseMatrix<double> MAC::constructCoefficientMatrix()
{
    size_t n = m_resolution;
    Eigen::SparseMatrix<double> m(n*n,n*n);

    auto tripletList = constructTriplets();

    m.setFromTriplets(tripletList.begin(), tripletList.end());
    return m;
}

std::vector<Eigen::Triplet<double>> MAC::constructTriplets()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::Triplet<double> t;
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
        {
            if (m_type[row][col] == "fluid")
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

    // For each fluid cell...
    for (size_t col = 0; col < m_resolution; ++col)
    {
        for (size_t row = 0; row < m_resolution; ++row)
        {
            if (m_type[row][col] == "fluid")
            {
                size_t i = index(row, col);
                double density = 1.0;
                float h = 1.0f;
                float divergence = 1.0f;
                size_t numNeighbourAirCells = 0;
                int atmosphericPressure = 101325;

                auto result = ((density*h)/_time)*divergence - numNeighbourAirCells*atmosphericPressure;

                v[i] = result;
            }
        }
    }

    return v;
}

size_t MAC::getType(size_t row, size_t col)
{
    if (outOfBounds(row, col))
    {
        return 0; // Return air (non-solid)
    }
    if (m_type[row][col] == "solid")
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
        if (m_type[pair.first][pair.second] != "fluid")
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
    return m_type[row][col] == "fluid";
}

ngl::Vec2 MAC::velocityAt(float _i, float _j)
{
    ngl::Vec2 v;

    const float &x=_i;
    const float &y=_j;
    const int row = floor(x);
    const int col = floor(y);

    if (outOfBounds(row, col))
    {
        return ngl::Vec2();
    }

    float x1 = m_x[row][col], x2 = 0.0f, x3 = 0.0f, x4 = 0.0f;
    float y1 = m_y[row][col], y2 = 0.0f, y3 = 0.0f, y4 = 0.0f;
    if (row < int(m_resolution-1))
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

void MAC::updateVectorField(float _time)
{
    std::cout << "Updating vector field.\n";
    applyConvection(_time);
    //    applyExternalForces(_time);
    //    applyViscosity(_time);
    //    applyPressure(_time);
//    calculatePressure(_time);
    moveParticles(_time);
}

void MAC::applyConvection(float _time)
{
    MAC tmp(m_resolution);
    for (float y = 0.5f; y < m_resolution; y+=1.0f)
    {
        for (size_t x = 0; x <= m_resolution; ++x)
        {
            ngl::Vec2 updated = traceParticle(x, y, _time);
            tmp.m_x[floor(y)][x] = updated.m_x;
        }
    }

    for (size_t y = 0; y <= m_resolution; ++y)
    {
        for (float x = 0.5f; x < m_resolution; x+=1.0f)
        {
            ngl::Vec2 updated = traceParticle(x, y, _time);
            tmp.m_y[y][floor(x)] = updated.m_y;
        }
    }

    tmp.fixBorderVelocities();
    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

void MAC::moveParticles(float _time)
{
    for (ngl::Vec2 &p : m_particles)
    {
        ngl::Vec2 velocity = velocityAt(p.m_x, p.m_y);
        p = p + _time*velocity;
    }
}

ngl::Vec2 MAC::traceParticle(float _x, float _y, float _time)
{
    // Trace particle from point (_x, _y) using simple forward Euler.
    // TODO: update to use RK2.
    ngl::Vec2 v = velocityAt(_x, _y);
    ngl::Vec2 prev_pos = ngl::Vec2(_x, _y) - _time*v;
    ngl::Vec2 prev_velocity = velocityAt(prev_pos.m_x, prev_pos.m_y);
    return prev_velocity;
}

void MAC::fixBorderVelocities()
{
    // Top row y.
    for (float &v : m_y[m_resolution])
    {
        v = 0.0f;
    }

    // Right column x.
    for (auto &row : m_x)
    {
        row[m_resolution] = 0.0f;
    }
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


void applyExternalForces(float _time)
{

}

void applyViscosity(float _time)
{

}

void applyPressure(float _time)
{

}

void moveMarkers(float _time)
{

}

float MAC::pressureDiff(size_t _x, size_t _y)
{
    return m_pressure.diff(_x, _y);
}

// Returns the staggered central difference for velocity
// at grid index x, y.
// Should take a position (ngl::Vec2) and bilinearly interpolate
// each component.
// Write test first.
// interpolate x component of each face of the grid cell containing the point
// interpolate y component of each face of the grid cell containing the point
// Maybe start with only points in the center of each grid.
ngl::Vec2 MAC::velocityDiff(size_t _x, size_t _y)
{
    ngl::Vec2 v;
    v.m_x = m_velocityX.diff(_x, _y);
    v.m_y = m_velocityY.diff(_x, _y);
    return v;
}

float MAC::Grid::at(size_t _i, size_t _j) const
{
    return m_v[index(_i, _j)];
}

void MAC::Grid::set(size_t _x, size_t _y, float _v)
{
    m_v[index(_x, _y)] = _v;
}

size_t MAC::Grid::index(size_t _x, size_t _y) const
{
    return _x*m_y + _y;
}

MAC::Grid::Grid(size_t _x, size_t _y)
{
    m_x = _x;
    m_y = _y;
    m_v = std::vector<float>(_x*_y);
}

bool MAC::Grid::operator==(const Grid &_other) const
{
    bool result = true;

    Grid o = _other;
    Grid t = *(this);

    if (m_x != _other.m_x) return false;
    if (m_y != _other.m_y) return false;

    for (Grid::iterator i = t.begin(), j=o.begin(); i!=t.end() && j!=o.end(); ++i, ++j)
    {
        result &= (*i == *j);
    }

    return result;
}
