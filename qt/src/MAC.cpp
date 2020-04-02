#include "MAC.h"

#include <stdexcept>
#include <string>
#include <iostream>

float MAC::Grid::at(size_t _i, size_t _j, size_t _k) const
{
    return m_v[_i][_j][_k];
}

MAC::Grid::Grid(size_t _resolution)
{
    m_resolution = _resolution;
//    std::vector<float> a(_resolution);
//    std::vector<std::vector<float>> b(_resolution);
    m_v = std::vector<std::vector<std::vector<float>>>(_resolution);
}
