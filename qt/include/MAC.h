#pragma once

#include <vector>
#include <gtest/gtest.h>

class MAC
{
public:
    MAC();
    MAC(size_t resolution);
    ~MAC() noexcept=default;

    class Grid
    {
    public:
        Grid();
        Grid(size_t _resolution);
        ~Grid() noexcept = default;
        float at(size_t _i, size_t _j, size_t _k) const;

    private:
        std::vector<std::vector<std::vector<float>>> m_v;
        size_t m_resolution;
        FRIEND_TEST(MACGrid, ctor);
    };

private:
    Grid pressure;
    Grid velocityX;
    Grid velocityY;
    Grid velocityZ;

    size_t m_width;
    size_t m_height;
    size_t m_depth;

    FRIEND_TEST(MAC, ctor);
};
