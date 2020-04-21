#include "MAC.h"

#include <gtest/gtest.h>

TEST(MAC, velocityAt)
{
    MAC grid(5);

    EXPECT_EQ(ngl::Vec2(1,1), grid.velocityAt(1,1));

    grid.m_x[0][0] = 2.0f;
    grid.m_y[0][0] = 2.0f;
    EXPECT_EQ(ngl::Vec2(2.0f,2.0f), grid.velocityAt(0,0));
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(1,1));
    EXPECT_EQ(ngl::Vec2(1.25f,1.25f), grid.velocityAt(0.5,0.5));

    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(2,2));
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(0,2));

    grid.updateVectorField(0.1f);
    grid.updateVectorField(0.1f);
    grid.updateVectorField(0.1f);
}

TEST(MAC, pressure)
{
    MAC grid(4);
    auto m = grid.constructCoefficientMatrix();
}

TEST(MAC, index)
{
    MAC grid(4);
    size_t expected = 0;
    EXPECT_EQ(expected, grid.index(0,0));
    expected = 15;
    EXPECT_EQ(expected, grid.index(3,3));
    expected = 4;
    EXPECT_EQ(expected, grid.index(1,0));
}

TEST(MAC, getType)
{
    MAC grid(3);
    EXPECT_EQ(0u, grid.getType(0,0));
    EXPECT_EQ(0u, grid.getType(2,2));
    EXPECT_EQ(0u, grid.getType(0,2));
    EXPECT_EQ(1u, grid.getType(1,1));
}

TEST(MAC, getNeighbours)
{
    MAC grid(3);
    std::map<size_t, size_t> m;
    for (size_t row = 0; row < 3; ++row)
    {
        for (size_t col = 0; col < 3; ++col)
        {
            m = grid.getNeighbours(row, col);
            for (const auto& e: m)
            {
                if (e.first == 4)
                {
                    EXPECT_EQ(e.second, 1u);
                }
                else
                {
                    EXPECT_EQ(e.second, 0u);
                }
            }
        }
    }
}

TEST(MAC, getNumNonLiquidNeighbours)
{
    MAC grid(5);
    size_t expected = 0;
    EXPECT_EQ(grid.getNumNonLiquidNeighbours(2,2), expected);
    expected = 2;
    EXPECT_EQ(grid.getNumNonLiquidNeighbours(0,0), expected);
}

TEST(MAC, constructDivergenceVector)
{
    MAC grid(5);
    auto v = grid.constructDivergenceVector(0.01f);
    grid.calculatePressure(0.1f);
}

TEST(MAC, getOwningCellIndex)
{
    MAC grid(4);
    size_t row, col;
    size_t expectRow=0, expectCol=0;
    grid.positionToCellIndex(-0.5,-0.5, row, col);
    EXPECT_EQ(expectRow, row);
    EXPECT_EQ(expectCol, col);
    expectRow=3;
    expectCol=3;
    grid.positionToCellIndex(0.49,0.49, row, col);
    EXPECT_EQ(expectRow, row);
    EXPECT_EQ(expectCol, col);
    expectRow=1;
    expectCol=0;
    grid.positionToCellIndex(-0.26,-0.1, row, col);
    EXPECT_EQ(expectRow, row);
    EXPECT_EQ(expectCol, col);
}

TEST(MAC, cellIndexToPosition)
{
    MAC grid(4);
    float x, y;
    float expectX=-0.5f, expectY=-0.375f;
    grid.cellIndexToPositionX(0,0,x,y);
    EXPECT_EQ(x, expectX);
    EXPECT_EQ(y, expectY);

    expectX=-0.375f, expectY=-0.5f;
    grid.cellIndexToPositionY(0,0,x,y);
    EXPECT_EQ(x, expectX);
    EXPECT_EQ(y, expectY);

    expectX=0.5f, expectY=0.375f;
    grid.cellIndexToPositionX(3,4,x,y);
    EXPECT_EQ(x, expectX);
    EXPECT_EQ(y, expectY);

    expectX=0.375f, expectY=0.5f;
    grid.cellIndexToPositionY(4,3,x,y);
    EXPECT_EQ(x, expectX);
    EXPECT_EQ(y, expectY);
}
