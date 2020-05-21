#include "MAC.h"

#include <gtest/gtest.h>
#include <set>

TEST(MAC, velocityAt)
{
    MAC grid(4);

    Index index{1,1};
    Position p{0,0};
    grid.cellIndexToPosition(index,p);
    EXPECT_EQ(Position(0,0), grid.velocityAtPosition(p));

    grid.m_x[2][2] = 2.0f;
    grid.m_y[2][2] = 2.0f;
    grid.cellIndexToPositionX({1,1},p);
    EXPECT_EQ(ngl::Vec2(0.0f,0.0f), grid.velocityAtPosition(p));

    grid.m_x[3][3] = 5.0f;
    grid.m_y[3][3] = 5.0f;
    grid.cellIndexToPositionY({3,3},p);

    EXPECT_EQ(ngl::Vec2(0.0f,0.0f), grid.velocityAtIndex({1,1}));
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAtIndex({2,2}));

    grid.updateVectorField();
    grid.updateVectorField();
    grid.updateVectorField();
}

TEST(MAC, getType)
{
    MAC grid(3);
    EXPECT_EQ(SOLID, grid.getType({0,0}));
    EXPECT_EQ(SOLID, grid.getType({2,2}));
    EXPECT_EQ(SOLID, grid.getType({0,2}));
    EXPECT_EQ(FLUID, grid.getType({1,1}));
}

TEST(MAC, getNumNonLiquidNeighbours)
{
    MAC grid(5);
    size_t expected = 0;
    EXPECT_EQ(grid.getNumNonLiquidNeighbours({2,2}), expected);
    expected = 2;
    EXPECT_EQ(grid.getNumNonLiquidNeighbours({0,0}), expected);
}

TEST(MAC, getOwningCellIndex)
{
    MAC grid(4);
    size_t expectRow=0, expectCol=0;
    Index index;
    grid.positionToCellIndex({-0.5,-0.5}, index);
    EXPECT_EQ(expectRow, index.row);
    EXPECT_EQ(expectCol, index.col);
    expectRow=3;
    expectCol=3;
    grid.positionToCellIndex({0.49,0.49}, index);
    EXPECT_EQ(expectRow, index.row);
    EXPECT_EQ(expectCol, index.col);
    expectRow=1;
    expectCol=0;
    grid.positionToCellIndex({-0.26,-0.1}, index);
    EXPECT_EQ(expectRow, index.row);
    EXPECT_EQ(expectCol, index.col);
}

TEST(MAC, cellIndexToPosition)
{
    MAC grid(4);
    grid.gridWidth=1.0f;
    Position p;
    float expectX=-0.5f, expectY=-0.375f;
    grid.cellIndexToPositionX({0,0},p);
    EXPECT_EQ(p.m_x, expectX);
    EXPECT_EQ(p.m_y, expectY);

    expectX=-0.375f, expectY=-0.5f;
    grid.cellIndexToPositionY({0,0},p);
    EXPECT_EQ(p.m_x, expectX);
    EXPECT_EQ(p.m_y, expectY);

    expectX=0.5f, expectY=0.375f;
    grid.cellIndexToPositionX({3,4},p);
    EXPECT_EQ(p.m_x, expectX);
    EXPECT_EQ(p.m_y, expectY);

    expectX=0.375f, expectY=0.5f;
    grid.cellIndexToPositionY({4,3},p);
    EXPECT_EQ(p.m_x, expectX);
    EXPECT_EQ(p.m_y, expectY);

    ngl::Vec2 v;
    grid.m_x[1][1] = 3.0f;
    grid.cellIndexToPositionX({0,0},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_x, 0.0f);
}

TEST(MAC, interpolate)
{
    MAC grid(3);
    Index index = {1,1};
    grid.m_y[1][1]=1.0f;
    grid.m_y[2][1]=3.0f;
    Position p;
    grid.cellIndexToPosition(index,p);
    float y = grid.interpolate(p, Dimension::y);
    EXPECT_FLOAT_EQ(y,2.0f);

    p.m_y-=0.5*grid.cellWidth;
    y = grid.interpolate(p, Dimension::y);
    EXPECT_FLOAT_EQ(y,1.0f);

    index={1,1};
    grid.cellIndexToPosition(index,p);
    p.m_x-=0.5*grid.cellWidth;
    y = grid.interpolate(p, Dimension::y);
    EXPECT_FLOAT_EQ(y,1.0f);

    grid.m_x[1][2]=3.0f;
    grid.m_x[1][1]=1.0f;
    index={1,1};
    grid.cellIndexToPosition(index,p);
    float x = grid.interpolate(p, Dimension::x);
    EXPECT_FLOAT_EQ(x,2.0f);

    p.m_x-=0.5*grid.cellWidth;
    x = grid.interpolate(p, Dimension::x);
    EXPECT_FLOAT_EQ(x,1.0f);

    index={1,1};
    grid.cellIndexToPositionX(index,p);
    p.m_y-=0.5*grid.cellWidth;
    x = grid.interpolate(p, Dimension::x);
    EXPECT_FLOAT_EQ(x,0.5f);
}

TEST(MAC, bordersFluidCell)
{
    MAC grid(4);
    grid.m_type[2][2] = AIR;
    EXPECT_EQ(grid.bordersFluidCellX({0,0}), false);
    EXPECT_EQ(grid.bordersFluidCellX({3,4}), false);
    EXPECT_EQ(grid.bordersFluidCellX({3,3}), false);
    EXPECT_EQ(grid.bordersFluidCellX({2,3}), false);
    EXPECT_EQ(grid.bordersFluidCellX({2,2}), true);
    EXPECT_EQ(grid.bordersFluidCellX({1,1}), true);

    EXPECT_EQ(grid.bordersFluidCellY({0,0}), false);
    EXPECT_EQ(grid.bordersFluidCellY({4,3}), false);
    EXPECT_EQ(grid.bordersFluidCellY({3,3}), false);
    EXPECT_EQ(grid.bordersFluidCellY({2,3}), false);
    EXPECT_EQ(grid.bordersFluidCellY({2,2}), true);
    EXPECT_EQ(grid.bordersFluidCellY({1,1}), true);
}
