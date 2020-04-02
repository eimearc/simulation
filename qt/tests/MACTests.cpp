#include "MAC.h"

#include <gtest/gtest.h>

TEST(MAC, ctor)
{
    MAC mac(3);
    EXPECT_EQ(mac.m_resolution, 3);
    MAC::Grid grid(3);
//    EXPECT_EQ(mac.m_pressure, grid);
//    EXPECT_EQ(mac.m_velocityX, grid);
//    EXPECT_EQ(mac.m_velocityY, grid);
//    EXPECT_EQ(mac.m_velocityZ, grid);
}

TEST(MACGrid, ctor)
{
    MAC::Grid grid(3);
    std::vector<std::vector<std::vector<float>>> expected =
    {
        { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } },
        { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } },
        { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }
    };

    EXPECT_EQ(grid.m_resolution, 3);
    EXPECT_EQ(grid.m_v, expected);
    EXPECT_EQ(grid.m_v[0].size(), 3);
    EXPECT_EQ(grid.m_v[0][0].size(), 3);
}

TEST(MACGrid, set)
{
    MAC::Grid grid(2);
    grid.set(0,0,0,1.0f);

    EXPECT_EQ(grid.m_v[0][0][0], 1.0f);
    EXPECT_EQ(grid.at(0,0,0), 1.0f);
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                if (i!=0 && j!=0 && k!=0)
                {
                    EXPECT_EQ(grid.at(i,j,k), 0.0f);
                }
            }
        }
    }
}

TEST(MACGrid, iterator)
{
    MAC::Grid grid(3);
    MAC::Grid::iterator i = grid.begin();
    size_t j;
    for (; i != grid.end() ; ++i)
    {
        std::cout << j++ << ":" << *i << std::endl;
    }
}

TEST(MAC, velocity)
{
//    MAC mac(3);
}
