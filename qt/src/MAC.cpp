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
}

ngl::Vec2 MAC::velocityAt(float _i, float _j)
{
    ngl::Vec2 v;

    const float &x=_i;
    const float &y=_j;
    int i = floor(x);
    int j = floor(y);

    float x1 = m_x[j][i];
    float x2 = m_x[j][i+1];
    float x3 = m_x[j+1][i];
    float x4 = m_x[j+1][i+1];

    v.m_x = (
        (i+1-x) * (j+1-y) * x1 +
        (x-i) * (j+1-y) * x2 +
        (i+1-x) * (y-j) * x3+
        (x-i) * (y-j) * x4
    );

    float y1 = m_y[j][i];
    float y2 = m_y[j][i+1];
    float y3 = m_y[j+1][i];
    float y4 = m_y[j+1][i+1];

    v.m_y = (
        (i+1-y) * (j+1-y) * y1 +
        (y-i) * (j+1-y) * y2 +
        (i+1-y) * (y-j) * y3 +
        (y-i) * (y-j) * y4
    );

    return v;
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
    v = getVelocity(_x+0.5*_time*v.m_x, _y+0.5*_time*v.m_y);
    return ngl::Vec2(_x, _y) + _time * v;
}

ngl::Vec2 MAC::getVelocity(float _x, float _y)
{
    ngl::Vec2 result =
    {
        getInterpolatedValueX(_x, _y),
        getInterpolatedValueY(_x, _y)
    };
    return result;
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
