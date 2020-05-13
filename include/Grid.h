#pragma once

#include <ngl/AbstractVAO.h>
#include <ngl/Vec3.h>
#include <vector>
#include <gtest/gtest.h>

class Grid
{
public:
    Grid()=default;
    Grid(size_t m_width, size_t m_height, size_t m_depth, float m_size);
    Grid(Grid &&_grid);
    Grid& operator=(Grid &&_other);
    bool operator==(const Grid &_other) const;
    ~Grid()=default;
   void drawOuter() const;
   void drawInner() const;
   void startCoords(ngl::Vec3 &_coords) const;
   float gridSize() const;
   float stepSize() const;
   size_t width() const;
   size_t height() const;
   size_t depth() const;

private:
    std::unique_ptr<ngl::AbstractVAO> m_outer_vao;
    std::vector<ngl::Vec3> m_outer_vbo;

    std::unique_ptr<ngl::AbstractVAO> m_inner_vao;
    std::vector<ngl::Vec3> m_inner_vbo;

    size_t m_width;
    size_t m_height;
    size_t m_depth;
    float m_size;
    float m_stepSize;

    void makeInnerVBO();
    void makeOuterVBO();

    FRIEND_TEST(Grid, startCoords);
    FRIEND_TEST(Grid, stepSize);
    FRIEND_TEST(Grid, lines);
};
