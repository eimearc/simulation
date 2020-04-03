#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>

MAC::MAC(size_t _resolution) :
    m_pressure(_resolution), m_velocityX(_resolution),
    m_velocityY(_resolution), m_velocityZ(_resolution),
    m_resolution(_resolution)
{
}

float MAC::Grid::at(size_t _i, size_t _j, size_t _k) const
{
    return m_v[index(_i, _j, _k)];
}

void MAC::Grid::set(size_t _i, size_t _j, size_t _k, float _v)
{
    m_v[index(_i, _j, _k)] = _v;
}

size_t MAC::Grid::index(size_t _x, size_t _y, size_t _z) const
{
    return _x*(m_x*m_y) + _y*(m_y) + _z;
}

MAC::Grid::Grid(size_t _x, size_t _y, size_t _z)
{
    m_x = _x;
    m_y = _y;
    m_z = _z;
    m_v = std::vector<float>(_x*_y*_z);
}

bool MAC::Grid::operator==(const Grid &_other) const
{
    bool result = true;

    Grid o = _other;
    Grid t = *(this);

    if (m_x != _other.m_x) return false;
    if (m_y != _other.m_y) return false;
    if (m_z != _other.m_z) return false;

    for (Grid::iterator i = t.begin(), j=o.begin(); i!=t.end() && j!=o.end(); ++i, ++j)
    {
        result &= (*i == *j);
    }

    return result;
}
