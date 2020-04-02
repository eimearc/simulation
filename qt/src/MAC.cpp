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
    return m_v[_i][_j][_k];
}

void MAC::Grid::set(size_t _i, size_t _j, size_t _k, float _v)
{
    m_v[_i][_j][_k] = _v;
}

MAC::Grid::Grid(size_t _resolution)
{
    m_resolution = _resolution;
    std::vector<float> a(_resolution);
    std::vector<std::vector<float>> b(_resolution);
    m_v = std::vector<std::vector<std::vector<float>>>(_resolution);

    for (std::vector<std::vector<float>> &i : m_v)
    {
        i = std::vector<std::vector<float>>(_resolution);
        for (std::vector<float> &j : i)
        {
            j = std::vector<float>(_resolution);
        }
    }
}

bool MAC::Grid::operator==(const Grid &_other)
{
    bool result = true;

    // TODO;

    return result;
}
