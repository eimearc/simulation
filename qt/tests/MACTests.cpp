#include "MAC.h"

#include <gtest/gtest.h>

//TEST(MAC, ctor)
//{
//    MAC mac(3);
//    EXPECT_EQ(mac.m_resolution, 3);
//    MAC::Grid grid(3);
//    EXPECT_EQ(mac.m_pressure, grid);
//    EXPECT_EQ(mac.m_velocityX, grid);
//    EXPECT_EQ(mac.m_velocityY, grid);
//    EXPECT_EQ(mac.m_velocityZ, grid);
//}

TEST(MACGrid, ctor)
{
    MAC::Grid grid(3);
    std::vector<float> expected = std::vector<float>(9);
    EXPECT_EQ(grid.m_x, 3);
    EXPECT_EQ(grid.m_y, 3);
    EXPECT_EQ(grid.m_v, expected);

    grid = MAC::Grid(3,2);
    expected = std::vector<float>(6);
    EXPECT_EQ(grid.m_x, 3);
    EXPECT_EQ(grid.m_y, 2);
    EXPECT_EQ(grid.m_v, expected);
}

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
    MAC grid(3);
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
    EXPECT_EQ(0, grid.getType(0,0));
    EXPECT_EQ(0, grid.getType(2,2));
    EXPECT_EQ(0, grid.getType(0,2));
    EXPECT_EQ(1, grid.getType(1,1));
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
                    EXPECT_EQ(e.second, 1);
                }
                else
                {
                    EXPECT_EQ(e.second, 0);
                }
            }
        }
    }
}

//TEST(MACGrid, set)
//{
//    MAC::Grid grid(2);
//    grid.set(0,0,0,1.0f);

//    EXPECT_EQ(grid.m_v[0], 1.0f);
//    EXPECT_EQ(grid.at(0,0,0), 1.0f);
//    for (int i = 0; i < 2; ++i)
//    {
//        for (int j = 0; j < 2; ++j)
//        {
//            for (int k = 0; k < 2; ++k)
//            {
//                if (i!=0 && j!=0 && k!=0)
//                {
//                    EXPECT_EQ(grid.at(i,j,k), 0.0f);
//                }
//            }
//        }
//    }
//}

//TEST(MACGrid, iterator)
//{
//    MAC::Grid grid(2);
//    grid.set(0,0,0,1);
//    grid.set(0,0,1,2);
//    grid.set(0,1,0,3);
//    grid.set(0,1,1,4);
//    grid.set(1,0,0,5);
//    grid.set(1,0,1,6);
//    grid.set(1,1,0,7);
//    grid.set(1,1,1,8);
//    size_t j=0;
//    for (auto i = grid.begin(); i != grid.end() ; ++i)
//    {
//        EXPECT_EQ(++j, *i);
//    }
//}

//TEST(MACGrid, index)
//{
//    MAC::Grid grid(2);
//    size_t got = grid.index(1,1,1);
//    EXPECT_EQ(got, 7);
//    got = grid.index(0,1,0);
//    EXPECT_EQ(got, 2);
//    got = grid.index(0,0,0);
//    EXPECT_EQ(got, 0);

//    grid = MAC::Grid(3);
//    got = grid.index(0,0,0);
//    EXPECT_EQ(got, 0);
//    got = grid.index(0,0,1);
//    EXPECT_EQ(got, 1);
//    got = grid.index(2,2,2);
//    EXPECT_EQ(got, 26);
//}

//TEST(MAC, velocity)
//{
////    MAC mac(3);
//}
