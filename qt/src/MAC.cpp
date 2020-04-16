#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>

MAC::MAC(size_t _resolution) :
    m_pressure(_resolution), m_velocityX(_resolution),
    m_velocityY(_resolution), m_resolution(_resolution)
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
