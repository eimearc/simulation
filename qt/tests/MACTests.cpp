#include "MAC.h"

#include <gtest/gtest.h>
#include <set>

TEST(MAC, velocityAt)
{
    MAC grid(5);

//    EXPECT_EQ(ngl::Vec2(1,1), grid.velocityAt(1,1));

//    grid.m_x[0][0] = 2.0f;
//    grid.m_y[0][0] = 2.0f;
//    EXPECT_EQ(ngl::Vec2(2.0f,2.0f), grid.velocityAt(0,0));
//    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(1,1));
//    EXPECT_EQ(ngl::Vec2(1.25f,1.25f), grid.velocityAt(0.5,0.5));

//    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(2,2));
//    EXPECT_EQ(ngl::Vec2(1.0f,1.0f), grid.velocityAt(0,2));

    grid.updateVectorField(0.1f);
    grid.updateVectorField(0.1f);
    grid.updateVectorField(0.1f);
}

//TEST(MAC, pressure)
//{
//    MAC grid(4);
//    auto m = grid.constructCoefficientMatrix();
//}

//TEST(MAC, index)
//{
//    MAC grid(4);
//    size_t expected = 0;
//    EXPECT_EQ(expected, grid.index(0,0));
//    expected = 15;
//    EXPECT_EQ(expected, grid.index(3,3));
//    expected = 4;
//    EXPECT_EQ(expected, grid.index(1,0));
//}

//TEST(MAC, getType)
//{
//    MAC grid(3);
//    EXPECT_EQ("SOLID", grid.getType(0,0));
//    EXPECT_EQ("SOLID", grid.getType(2,2));
//    EXPECT_EQ("SOLID", grid.getType(0,2));
//    EXPECT_EQ("FLUID", grid.getType(1,1));
//}

//TEST(MAC, getNeighbours)
//{
//    MAC grid(3);
//    std::map<size_t, std::string> m;
//    for (size_t row = 0; row < 3; ++row)
//    {
//        for (size_t col = 0; col < 3; ++col)
//        {
//            m = grid.getNeighbourType(row, col);
//            for (const auto& e: m)
//            {
//                if (e.first == 4)
//                {
//                    EXPECT_EQ(e.second, "FLUID");
//                }
//                else
//                {
//                    EXPECT_EQ(e.second, "SOLID");
//                }
//            }
//        }
//    }
//}

//TEST(MAC, getNumNonLiquidNeighbours)
//{
//    MAC grid(5);
//    size_t expected = 0;
//    EXPECT_EQ(grid.getNumNonLiquidNeighbours(2,2), expected);
//    expected = 2;
//    EXPECT_EQ(grid.getNumNonLiquidNeighbours(0,0), expected);
//}

//TEST(MAC, constructDivergenceVector)
//{
//    MAC grid(5);
//    auto v = grid.constructDivergenceVector(0.01f);
//    grid.calculatePressure(0.1f);
//}

//TEST(MAC, getOwningCellIndex)
//{
//    MAC grid(4);
//    size_t row, col;
//    size_t expectRow=0, expectCol=0;
//    grid.positionToCellIndex(-0.5,-0.5, row, col);
//    EXPECT_EQ(expectRow, row);
//    EXPECT_EQ(expectCol, col);
//    expectRow=3;
//    expectCol=3;
//    grid.positionToCellIndex(0.49,0.49, row, col);
//    EXPECT_EQ(expectRow, row);
//    EXPECT_EQ(expectCol, col);
//    expectRow=1;
//    expectCol=0;
//    grid.positionToCellIndex(-0.26,-0.1, row, col);
//    EXPECT_EQ(expectRow, row);
//    EXPECT_EQ(expectCol, col);
//}

TEST(MAC, cellIndexToPosition)
{
    MAC grid(4);
    grid.gridWidth=1.0f;
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

    ngl::Vec2 v;
    grid.m_x[1][1] = 3.0f;
    grid.cellIndexToPositionX(0,0,x,y);
    v = grid.velocityAt(x,y);
    EXPECT_EQ(v.m_x, 0.0f);

    grid.m_x[0][0] = 3.0f;
    grid.cellIndexToPositionX(0,0,x,y);
    v = grid.velocityAt(x,y);
    EXPECT_EQ(v.m_x, 3.0f);

    grid.m_x[3][4] = 5.0f;
    grid.cellIndexToPositionX(3,4,x,y);
    v = grid.velocityAt(x,y);
    EXPECT_EQ(v.m_x, 5.0f);

    grid.m_y[0][0] = 7.0f;
    grid.cellIndexToPositionX(0,0,x,y);
    v = grid.velocityAt(x,y);
    EXPECT_EQ(v.m_y, 7.0f);

    grid.m_y[4][3] = 12.0f;
    grid.cellIndexToPositionX(4,3,x,y);
    v = grid.velocityAt(x,y);
    EXPECT_EQ(v.m_y, 12.0f);

    grid.m_y[1][3] = 12.0f;
    grid.cellIndexToPositionX(1,3,x,y);
    v = grid.velocityAt(x,y);
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

//TEST(MAC, bordersFluidCell)
//{
//    MAC grid(4);
//    grid.m_type[2][2] = "AIR";
//    EXPECT_EQ(grid.bordersFluidCellX(0,0), false);
//    EXPECT_EQ(grid.bordersFluidCellX(3,4), false);
//    EXPECT_EQ(grid.bordersFluidCellX(3,3), false);
//    EXPECT_EQ(grid.bordersFluidCellX(2,3), false);
//    EXPECT_EQ(grid.bordersFluidCellX(2,2), true);
//    EXPECT_EQ(grid.bordersFluidCellX(1,1), true);

//    EXPECT_EQ(grid.bordersFluidCellY(0,0), false);
//    EXPECT_EQ(grid.bordersFluidCellY(4,3), false);
//    EXPECT_EQ(grid.bordersFluidCellY(3,3), false);
//    EXPECT_EQ(grid.bordersFluidCellY(2,3), false);
//    EXPECT_EQ(grid.bordersFluidCellY(2,2), true);
//    EXPECT_EQ(grid.bordersFluidCellY(1,1), true);
//}
