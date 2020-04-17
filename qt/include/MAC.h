#pragma once

#include <vector>
#include <gtest/gtest.h>
#include <ngl/Vec2.h>
#include <ngl/NGLStream.h>

typedef ngl::Vec2 Index;

class MAC
{
public:
    MAC()=default;
    MAC(size_t _resolution);
    ~MAC() noexcept=default;

    ngl::Vec2 velocityAt(float _i, float _j);

    void advance(float _time);

    void applyConvection(float _time);
    ngl::Vec2 traceParticle(float _x, float _y, float _time);
    ngl::Vec2 getVelocity(float _x, float _y);
    float getInterpolatedValueX(float _x, float _y);
    float getInterpolatedValueY(float _x, float _y);

    void updateVectorField();

    void applyExternalForces(float _time);
    void applyViscosity(float _time);
    void applyPressure(float _time);

    void moveMarkers(float _time);

    ngl::Vec2 velocityDiff(size_t _x, size_t _y);
    float pressureDiff(size_t _x, size_t _y);

    std::vector<std::vector<float>> m_x;
    std::vector<std::vector<float>> m_y;
    std::vector<std::vector<std::string>> m_type;

    class Grid
    {
    public:
        Grid()=default;
        Grid(size_t _x, size_t _y);
        Grid(size_t _resolution) : Grid(_resolution, _resolution) {}
        ~Grid() noexcept = default;
        float at(size_t _x, size_t _y) const;
        virtual float diff(size_t _x, size_t _y) const {};
        void set(size_t _x, size_t _y, float _v);
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
            return iterator(&(m_v[0])+(m_x)*(m_y));
        }

    private:
        size_t index(size_t _x, size_t _y) const;

        std::vector<float> m_v;
        size_t m_x;
        size_t m_y;
        FRIEND_TEST(MACGrid, ctor);
        FRIEND_TEST(MACGrid, set);
        FRIEND_TEST(MACGrid, index);
    };

    class PressureGrid : Grid
    {
    public:
        PressureGrid()=default;
        PressureGrid(size_t _resolution) : Grid(_resolution, _resolution) {}
        ~PressureGrid() noexcept = default;
        float diff(size_t _x, size_t _y) const override
        {
            float result = at(_x, _y) - at(_x, _y);
            return result;
        }
    };

    class GridX : Grid
    {
    public:
        GridX()=default;
        GridX(size_t _resolution) : Grid(_resolution+1, _resolution) {}
        ~GridX() noexcept = default;
        float diff(size_t _x, size_t _y) const override
        {
            float result = at(_x+1, _y) - at(_x, _y);
            return result;
        }
    };

    class GridY : Grid
    {
    public:
        GridY()=default;
        GridY(size_t _resolution) : Grid(_resolution, _resolution+1) {}
        ~GridY() noexcept = default;
        float diff(size_t _x, size_t _y) const override
        {
            float result = at(_x, _y+1) - at(_x, _y);
            return result;
        }
    };

private:
    float velocityX(size_t _x, size_t _y, size_t _z);
    float velocityY(size_t _x, size_t _y, size_t _z);

    PressureGrid m_pressure;
//    Grid m_velocity;
    GridX m_velocityX;
    GridY m_velocityY;
    size_t m_resolution;

    FRIEND_TEST(MAC, ctor);
    FRIEND_TEST(MAC, velocityAt);
};

std::ostream& operator<<(std::ostream& os, MAC& mac);
