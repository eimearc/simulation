#pragma once

#include <vector>
#include <gtest/gtest.h>
#include <ngl/Vec3.h>

class MAC
{
public:
    MAC();
    MAC(size_t _resolution);
    ~MAC() noexcept=default;

    ngl::Vec3 velocity(size_t _x, size_t _y, size_t _z);
    ngl::Vec3 pressure(size_t _x, size_t _y, size_t _z);

    class Grid
    {
    public:
        Grid()=default;
        Grid(size_t _resolution);
        ~Grid() noexcept = default;
        float at(size_t _i, size_t _j, size_t _k) const;
        void set(size_t _i, size_t _j, size_t _k, float _v);
        bool operator==(const Grid &_other);

        struct iterator
        {
        public:
            iterator(float* f) : m_ptr(f) {}
            iterator operator++(){iterator i=*this; m_ptr++; return i;}
            float& operator*(){return *m_ptr;}
            float* operator&(){return m_ptr;}
            bool operator==(const iterator &_rhs){return m_ptr == _rhs.m_ptr;}
            bool operator!=(const iterator &_rhs){return m_ptr != _rhs.m_ptr;}

        private:
            float* m_ptr;
        };

        iterator begin() {return iterator(&m_v[0][0][0]);}
        iterator end() {return iterator(&(m_v[0][0][0])+sizeof(float)*(m_resolution-1)*(m_resolution-1)*(m_resolution-1));}

    private:
        std::vector<std::vector<std::vector<float>>> m_v;
        size_t m_resolution;
        FRIEND_TEST(MACGrid, ctor);
        FRIEND_TEST(MACGrid, set);
    };

private:
    float velocityX(size_t _x, size_t _y, size_t _z);
    float velocityY(size_t _x, size_t _y, size_t _z);
    float velocityZ(size_t _x, size_t _y, size_t _z);

    Grid m_pressure;
    Grid m_velocityX;
    Grid m_velocityY;
    Grid m_velocityZ;
    size_t m_resolution;

    FRIEND_TEST(MAC, ctor);
};
