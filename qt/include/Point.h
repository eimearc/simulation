#pragma once

#include <ngl/Vec3.h>
#include <ngl/NGLStream.h>

class Point
{
public:
    Point() = default;
    Point(ngl::Vec3 _position, ngl::Vec3 _direction, ngl::Vec3 _velocity);
    ~Point() = default;

    void update();
    ngl::Vec3 direction() const;
    ngl::Vec3 position() const;
    ngl::Vec3 velocity() const;

private:
    ngl::Vec3 m_position;
    ngl::Vec3 m_direction;
    ngl::Vec3 m_velocity;
};
