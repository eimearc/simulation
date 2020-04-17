#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>

MAC::MAC(size_t _resolution) :
    m_pressure(_resolution), m_velocityX(_resolution),
    m_velocityY(_resolution), m_resolution(_resolution)
{
    m_x = std::vector<std::vector<float>>(m_resolution, std::vector<float>(m_resolution+1, 1.0f));
    m_y = std::vector<std::vector<float>>(m_resolution+1, std::vector<float>(m_resolution, 1.0f));
    m_type = std::vector<std::vector<std::string>>(m_resolution, std::vector<std::string>(m_resolution, "fluid"));
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
}

ngl::Vec2 MAC::velocityAt(float _i, float _j)
{
    ngl::Vec2 v;

    std::cout << "\t x:" << _i << " y:" << _j << std::endl;

    const float &x=_i;
    const float &y=_j;
    int i = floor(x);
    int j = floor(y);
    const int &x_floor = i;
    const int &y_floor = j;

    if (x_floor>=(int(m_resolution)) && y_floor>=(int(m_resolution)))
    {
        // Top right corner of the grid.
        return ngl::Vec2(0.0f,0.0f);
    }

    if (y_floor > int(m_resolution-1))
    {
        // Above the top row of the grid.
        v.m_x = 0.0f;
    }
    else if (x_floor > int(m_resolution-1))
    {
        v.m_y = 0.0f;
    }
    else
    {
        float x1 = m_x[j][i], x2 = 0.0, x3 = 0.0f, x4 = 0.0f;
        if (y_floor < int(m_resolution-1))
        {
            // Top row of the grid.
            x3 = m_x[j+1][i];
            if (x_floor < int(m_resolution))
            {
                x2 = m_x[j][i+1];
                x4 = m_x[j+1][i+1];
            }
        }
        std::cout << "\tX: " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";
        v.m_x = (
            (i+1-x) * (j+1-y) * x1 +
            (x-i) * (j+1-y) * x2 +
            (i+1-x) * (y-j) * x3+
            (x-i) * (y-j) * x4
        );

        float y1 = m_y[j][i], y2 = 0.0, y3 = 0.0f, y4 = 0.0f;
        if (x_floor < int(m_resolution-1))
        {
            // Top row of the grid.
            y3 = m_y[j][i+1];
            if (y_floor < int(m_resolution))
            {
                y2 = m_y[j+1][i];
                y4 = m_y[j+1][i+1];
            }
        }
        std::cout << "\tY: " << y1 << ", " << y2 << ", " << y3 << ", " << y4 << "\n";
        v.m_y = (
            (i+1-x) * (j+1-y) * y1 +
            (x-i) * (j+1-y) * y2 +
            (i+1-x) * (y-j) * y3+
            (x-i) * (y-j) * y4
        );
    }

    std::cout << "\t" << v << '\n';

    return v;
}

void MAC::updateVectorField()
{
    std::cout << "Updating vector field.\n";
    // For every grid point.
    // Move particle back.
    // Find old velocity.
    // Update new velocity to old one.
    MAC tmp(m_resolution);
    for (float y = 0.5f; y < m_resolution; y+=1.0f)
    {
        for (size_t x = 0; x < m_resolution; ++x)
        {
            std::cout << "X Update: x:" << x << " y: " << y << std::endl;
            ngl::Vec2 updated = traceParticle(x, y, 1.0f);
            tmp.m_x[floor(y)][x] = updated.m_x;
        }
    }

    for (size_t y = 0; y < m_resolution; ++y)
    {
        for (float x = 0.5f; x < m_resolution; x+=1.0f)
        {
            std::cout << "Y Update: " << x << " y: " << y << std::endl;
            ngl::Vec2 updated = traceParticle(x, y, 1.0f);
            tmp.m_y[y][floor(x)] = updated.m_y;
        }
    }

    std::cout << '\n';

//    std::cout << "Old grid\n" << *this << std::endl;

//    std::cout << "Udpated tmp\n" << tmp << std::endl;

    m_x = tmp.m_x;
    m_y = tmp.m_y;
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

    for (int i = mac.m_x.size()-1; i >= 0; --i)
    {
        os << i << ":\t";
        for (const auto &y: mac.m_x[i])
        {
            os << y << ' ';
        }
        os << '\n';
    }
    os << '\n';

    for (int i = mac.m_y.size()-1; i >= 0; --i)
    {
        os << i << ":\t";
        for (const auto &y: mac.m_y[i])
        {
            os << y << ' ';
        }
        os << '\n';
    }
    os << '\n';

    return os;
}

void MAC::advance(float _time)
{
//    applyConvection(_time);
//    applyExternalForces(_time);
//    applyViscosity(_time);
//    applyPressure(_time);
}

void MAC::applyConvection(float _time)
{
    GridX tmpX(m_resolution);
    GridY tmpY(m_resolution);
    Grid tmp(m_resolution);
    // For each grid point (each grid cell),
    // for each x, y, traceParticle and then update.
    // Store temp copies and then copy across once whole step is complete.
    for (int i = 0; i < m_resolution; ++i)
    {
        for (int j = 0; j < m_resolution; ++j)
        {

        }
    }
}

ngl::Vec2 MAC::traceParticle(float _x, float _y, float _time)
{
    // Trace particle from point (_x, _y) using RK2.
    ngl::Vec2 v = getVelocity(_x, _y);
    float x_pos = _x+0.5*_time*v.m_x;
    float y_pos = _y+0.5*_time*v.m_y;
    v = getVelocity(x_pos, y_pos);
    return ngl::Vec2(_x, _y) + _time * v;
}

ngl::Vec2 MAC::getVelocity(float _x, float _y)
{
    return velocityAt(_x, _y);
//    ngl::Vec2 result =
//    {
//        getInterpolatedValueX(_x, _y),
//        getInterpolatedValueY(_x, _y)
//    };
//    return result;
}

float MAC::getInterpolatedValueX(float _x, float _y)
{
    int i = floor(_x);
    int j = floor(_y);
    return 0.0f;
}

float MAC::getInterpolatedValueY(float _x, float _y)
{
    return 0.0f;
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
