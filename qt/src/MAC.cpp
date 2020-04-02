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
    return m_v[_i*(m_resolution+m_resolution) + _j*(m_resolution) + _k];
}

void MAC::Grid::set(size_t _i, size_t _j, size_t _k, float _v)
{
    m_v[_i*(m_resolution+m_resolution) + _j*(m_resolution) + _k] = _v;
}

MAC::Grid::Grid(size_t _resolution)
{
    m_resolution = _resolution;
    std::vector<float> a(_resolution);
    std::vector<std::vector<float>> b(_resolution);
    m_v = std::vector<float>(_resolution*_resolution*_resolution);
}

bool MAC::Grid::operator==(const Grid &_other) const
{
    bool result = true;

    Grid o = _other;
    Grid t = *(this);

    if (m_resolution != _other.m_resolution) return false;

    for (Grid::iterator i = t.begin(), j=o.begin(); i!=t.end() && j!=o.end(); ++i, ++j)
    {
        result &= (*i == *j);
    }

    return result;
}
