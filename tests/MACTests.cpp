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

    grid.m_x[0][0] = 2.0f;
    grid.m_y[0][0] = 2.0f;
    grid.cellIndexToPosition({0,0},p);
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAtPosition(p));

    grid.cellIndexToPositionX({0,0},p);
    EXPECT_EQ(ngl::Vec2(2.0f,1.5f), grid.velocityAtPosition(p));

    grid.m_x[2][2] = 2.0f;
    grid.m_y[2][2] = 2.0f;
    grid.cellIndexToPositionX({1,1},p);
    EXPECT_EQ(ngl::Vec2(0.0f,0.0f), grid.velocityAtPosition(p));

    grid.m_x[3][3] = 5.0f;
    grid.m_y[3][3] = 5.0f;
    grid.cellIndexToPositionY({3,3},p);
    auto expected = grid.velocityAtPosition(p);

    EXPECT_EQ(5.0f, expected.m_y);
    EXPECT_EQ(ngl::Vec2(0.0f,0.0f), grid.velocityAtIndex({1,1}));
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAtIndex({0,0}));
    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAtIndex({2,2}));
    grid.m_x[2][0]=3.0f;
    EXPECT_EQ(ngl::Vec2(1.5f,0.0f), grid.velocityAtIndex({2,0}));

    grid.updateVectorField();
    grid.updateVectorField();
    grid.updateVectorField();
}

TEST(MAC, moveParticles)
{
    MAC grid(5);
    float time = 0.001f;
    Position p;
    grid.m_x[2][1] = 0.1f;
    grid.m_x[2][2] = 0.1f;
    grid.cellIndexToPosition({2,2}, p);
    Velocity v =  grid.traceParticle(p, time);
    std::cout << v << p-time*v << std::endl;
}

//TEST(MAC, pressure)
//{
//    MAC grid(4);
//    auto m = grid.constructCoefficientMatrix();
//}

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

//TEST(MAC, constructDivergenceVector)
//{
//    MAC grid(5);
//    auto v = grid.constructDivergenceVector(0.01f);
//    grid.calculatePressure(0.1f);
//}

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

    std::cout << grid.m_x << std::endl;

    grid.m_x[0][0] = 3.0f;
    grid.cellIndexToPositionX({0,0},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_x, 3.0f);

    grid.m_x[3][4] = 5.0f;
    grid.cellIndexToPositionX({3,4},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_x, 5.0f);

    grid.m_y[0][0] = 7.0f;
    grid.cellIndexToPositionY({0,0},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_y, 7.0f);

    grid.m_y[4][3] = 12.0f;
    grid.cellIndexToPositionY({4,3},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_y, 12.0f);

    grid.m_y[1][3] = 12.0f;
    grid.cellIndexToPositionY({1,3},p);
    v = grid.velocityAtPosition(p);
    EXPECT_EQ(v.m_y, 12.0f);
}

//TEST(MAC, constructCoefficientMatrix)
//{
//    MAC grid(4);

//    std::set<size_t> all = {};
//    for (size_t i = 0; i < 4*4; ++i)
//    {
//        all.insert(i);
//    }
//    std::set<size_t> corner = {0,3,12,15};
//    std::set<size_t> border = {0,1,2,3,4,7,8,11,12,13,14,15};

//    auto m = grid.constructCoefficientMatrix();

//    for (size_t i = 0; i < 4*4; ++i)
//    {
//        for (size_t k = 0; k < 4*4; ++k)
//        {
//            if(border.find(i)!=border.end())
//            {
//                ASSERT_EQ(m.coeffRef(i,k), 0);
//            }
//            else
//            {
//                if (i==k)
//                {
//                    ASSERT_GT(m.coeffRef(i,k),0);
//                }
//            }
//        }
//    }
//}

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

TEST(MAC, borderVelocities)
{
    MAC grid(10);
    auto &types = grid.m_type;
    int start=5;
    for (int i = 0; i < 3; ++i)
    {
        grid.m_type[i+start][start]=SOLID;
        grid.m_type[i+start][start+1]=SOLID;
        grid.m_type[i+start][start+2]=SOLID;
    }

    for (int i=0; i<30; ++i)
    {
        grid.updateVectorField();
        Index index;
        for (index.row=0;index.row<grid.m_resolution-1;++index.row)
        {
            for (index.col=0;index.col<grid.m_resolution-1;++index.col)
            {
                int row = index.row;
                int col = index.col;
                if (types[row][col]==SOLID && types[row][col+1]==SOLID)
                    EXPECT_FLOAT_EQ(grid.m_x[row][col], 0.0f);
                if (types[row][col]==SOLID && types[row+1][col]==SOLID)
                    EXPECT_FLOAT_EQ(grid.m_y[row][col], 0.0f);
            }
        }
    }
    for (const auto &p: grid.m_particles)
    {
        std::cout << p << std::endl;
        EXPECT_FALSE(grid.isOutsideFluid(p));
    }
}
