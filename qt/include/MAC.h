#pragma once

#include <vector>
#include <gtest/gtest.h>
#include <ngl/Vec2.h>
#include <ngl/NGLStream.h>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <map>
#include <ngl/VAOFactory.h>

class MAC
{
public:
    MAC()=default;
    MAC(MAC&& other);
    MAC(const MAC& other);
    MAC& operator=(MAC&& other);
    MAC(size_t _resolution);
    ~MAC() noexcept=default;

    void draw(float _time);

    void updateVectorField(float _time);
    void applyConvection(float _time);
    void applyExternalForces(float _time);
    void applyViscosity(float _time);
    void calculatePressure(float _time);
    void applyPressure(float _time);
    void moveParticles(float _time);

    struct NeighbourTypes;

private:
    // Velocity Methods
    ngl::Vec2 velocityAt(const float x, const float y);
    ngl::Vec2 velocityAt(size_t row, size_t col);
    ngl::Vec2 traceParticle(float _x, float _y, float _time);
    void fixBorderVelocities();

    // Pressure Methods
    float calculateModifiedDivergence(size_t row, size_t col);
    std::vector<Eigen::Triplet<double>> constructNeighbourTriplets();
    Eigen::VectorXd constructDivergenceVector(float _time);
    Eigen::SparseMatrix<double> constructCoefficientMatrix();

    // Drawing Methods
    void setupVAO();
    void setupVBO();
    void updateVBO();

    // Helper Methods
    std::string getType(size_t row, size_t col);
    bool isSolidCell(size_t row, size_t col);
    bool isFluidCell(size_t row, size_t col);
    bool isAirCell(size_t row, size_t col);
    bool isOutsideGrid(ngl::Vec2 pos);
    size_t index(size_t row, size_t col);
    void coordinate(size_t index, size_t &row, size_t &col);
    void positionToCellIndex(float x, float y, size_t &row, size_t &col);
    void cellIndexToPositionX(size_t row, size_t col, float &x, float &y);
    void cellIndexToPositionY(size_t row, size_t col, float &x, float &y);
    void cellIndexToPosition(size_t row, size_t col, float &x, float &y);
    bool outOfBounds(size_t row, size_t col);
    std::map<size_t, std::string> getNeighbourType(size_t row, size_t col);
    size_t getNumNonLiquidNeighbours(size_t row, size_t col);
    size_t getNumNonSolidNeighbours(size_t row, size_t col);
    std::vector<std::pair<size_t, size_t>> getNeighbourIndices(size_t row, size_t col);
    ngl::Vec2 applyPressureToPoint(float x, float y, float _time);
    bool bordersSolidCellX(size_t row, size_t col);
    bool bordersSolidCellY(size_t row, size_t col);
    bool bordersFluidCellX(size_t row, size_t col);
    bool bordersFluidCellY(size_t row, size_t col);
    size_t numFluidCells();

    ngl::Vec2 calculatePressureGradient(size_t row, size_t col);
    void updateGrid();

    std::vector<std::vector<float>> m_x;
    std::vector<std::vector<float>> m_y;
    std::vector<std::vector<float>> m_pressure;
    std::vector<std::vector<std::string>> m_type;
    std::vector<std::vector<size_t>> m_numParticles;
    std::vector<ngl::Vec2> m_particles;
    std::vector<std::vector<int>> m_indices;
    size_t m_resolution;
    float gridWidth = 1;
    float cellWidth;
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec2> m_vbo;

    FRIEND_TEST(MAC, ctor);
    FRIEND_TEST(MAC, velocityAt);
    FRIEND_TEST(MAC, pressure);
    FRIEND_TEST(MAC, index);
    FRIEND_TEST(MAC, getType);
    FRIEND_TEST(MAC, getNeighbours);
    FRIEND_TEST(MAC, getNumNonLiquidNeighbours);
    FRIEND_TEST(MAC, constructDivergenceVector);
    FRIEND_TEST(MAC, getOwningCellIndex);
    FRIEND_TEST(MAC, cellIndexToPosition);
    FRIEND_TEST(MAC, constructCoefficientMatrix);
    FRIEND_TEST(MAC, bordersFluidCell);
    friend std::ostream& operator<<(std::ostream& os, MAC& mac);
};

std::ostream& operator<<(std::ostream& os, MAC& mac);
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<float>>& grid);
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<size_t>>& grid);
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<int>>& grid);
