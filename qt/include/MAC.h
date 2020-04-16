#pragma once

#include <vector>
#include <gtest/gtest.h>
#include <ngl/Vec3.h>

typedef ngl::Vec3 Index;

class MAC
{
public:
    MAC();
    MAC(size_t _resolution);
    ~MAC() noexcept=default;

    ngl::Vec3 velocityDiff(size_t _x, size_t _y, size_t _z);
    float pressureDiff(size_t _x, size_t _y, size_t _z);

    class Grid
    {
    public:
        Grid()=default;
        Grid(size_t _x, size_t _y, size_t _z);
        Grid(size_t _resolution) : Grid(_resolution, _resolution, _resolution) {}
        ~Grid() noexcept = default;
        float at(size_t _i, size_t _j, size_t _k) const;
        virtual float diff(size_t _i, size_t _j, size_t _k) const {};
        void set(size_t _i, size_t _j, size_t _k, float _v);
        bool operator==(const Grid &_other) const;

        struct iterator
        {
        public:
            iterator(float* f) : m_ptr(f) {}
            iterator operator++(){
                iterator i=*this;
                m_ptr++;
                return i;
            }
            float& operator*(){return *m_ptr;}
            float* operator&(){return m_ptr;}
            bool operator==(const iterator &_rhs){return m_ptr == _rhs.m_ptr;}
            bool operator!=(const iterator &_rhs){return m_ptr != _rhs.m_ptr;}

        private:
            float* m_ptr;
        };

        iterator begin() {
            return iterator(&m_v[0]);
        }
        iterator end() {
            return iterator(&(m_v[0])+(m_x)*(m_y)*(m_z));
        }

    private:
        size_t index(size_t _x, size_t _y, size_t _z) const;

        std::vector<float> m_v;
        size_t m_x;
        size_t m_y;
        size_t m_z;
        FRIEND_TEST(MACGrid, ctor);
        FRIEND_TEST(MACGrid, set);
        FRIEND_TEST(MACGrid, index);

        // TODO: need a step size for delta(pos).
    };

    class PressureGrid : Grid
    {
    public:
        PressureGrid()=default;
        PressureGrid(size_t _resolution) : Grid(_resolution, _resolution, _resolution) {}
        ~PressureGrid() noexcept = default;
        float diff(size_t _x, size_t _y, size_t _z) const override
        {
            float result = at(_x, _y, _z) - at(_x, _y, _z);
            return result;
        }
    };

    class GridX : Grid
    {
    public:
        GridX()=default;
        GridX(size_t _resolution) : Grid(_resolution+1, _resolution, _resolution) {}
        ~GridX() noexcept = default;
        float diff(size_t _x, size_t _y, size_t _z) const override
        {
            float result = at(_x+1, _y, _z) - at(_x, _y, _z);
            return result;
        }
    };

    class GridY : Grid
    {
    public:
        GridY()=default;
        GridY(size_t _resolution) : Grid(_resolution, _resolution+1, _resolution) {}
        ~GridY() noexcept = default;
        float diff(size_t _x, size_t _y, size_t _z) const override
        {
            float result = at(_x, _y+1, _z) - at(_x, _y, _z);
            return result;
        }
    };

    class GridZ : Grid
    {
    public:
        GridZ()=default;
        GridZ(size_t _resolution) : Grid(_resolution, _resolution, _resolution+1) {}
        ~GridZ() noexcept = default;
        float diff(size_t _x, size_t _y, size_t _z) const override
        {
            float result = at(_x, _y, _z+1) - at(_x, _y, _z);
            return result;
        }
    };

private:
    float velocityX(size_t _x, size_t _y, size_t _z);
    float velocityY(size_t _x, size_t _y, size_t _z);
    float velocityZ(size_t _x, size_t _y, size_t _z);

    PressureGrid m_pressure;
    GridX m_velocityX;
    GridY m_velocityY;
    GridZ m_velocityZ;
    size_t m_resolution;

    FRIEND_TEST(MAC, ctor);
};
