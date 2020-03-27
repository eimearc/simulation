#include "Point.h"

Point::Point(ngl::Vec3 _position, ngl::Vec3 _direction, ngl::Vec3 _velocity)
{
    m_position = _position;
    m_direction = _direction;
    m_velocity = _velocity;
}

void Point::update()
{
    // Rotate around the z axis.
    float theta = 1.0f;
    ngl::Real x = m_direction.m_x;
    ngl::Real y = m_direction.m_y;
    m_direction.m_x = x*cos(theta) - y*sin(theta);
    m_direction.m_y = x*sin(theta) + y*cos(theta);
}

ngl::Vec3 Point::position() const
{
    return m_position;
}

ngl::Vec3 Point::direction() const
{
    return m_direction;
}

ngl::Vec3 Point::velocity() const
{
    return m_velocity;
}
