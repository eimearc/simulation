#pragma once

#include <vector>
#include <gtest/gtest.h>
#include <ngl/Vec2.h>
#include <ngl/Vec3.h>
#include <ngl/NGLStream.h>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <map>
#include <ngl/VAOFactory.h>

struct Index
{
    int row;
    int col;
};

typedef ngl::Vec2 Position;
typedef ngl::Vec2 Velocity;

enum Dimension {x,y};

class MAC
{
public:
    MAC()=default;
    MAC(MAC&& other);
    MAC(const MAC& other);
    MAC& operator=(MAC&& other);
    MAC(size_t _resolution);
    ~MAC() noexcept=default;

    void update(float _time);
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
    Velocity velocityAtPosition(const Position p);
    Velocity velocityAtIndex(const Index index);
    Velocity traceParticle(const Position &p, float _time);
    void fixBorderVelocities();

    // Pressure Methods
    float calculateModifiedDivergence(size_t row, size_t col);
    std::vector<Eigen::Triplet<double>> constructNeighbourTriplets();
    Eigen::VectorXd constructDivergenceVector(float _time);
    Eigen::SparseMatrix<double> constructCoefficientMatrix();
    Velocity calculatePressureGradient(size_t row, size_t col);
    Velocity applyPressureToPoint(const Index &index, float _time, Dimension dimension);

    // Drawing Methods
    void setupVAO();
    void setupVBO();
    void updateVBO();

    // Helper Methods
    std::string getType(const Index &index);
    bool isSolidCell(const Index &index);
    bool isFluidCell(const Index &index);
    bool isAirCell(size_t row, size_t col);
    bool isOutsideGrid(const Position &p);
    size_t vectorIndex(size_t row, size_t col);
    void coordinate(size_t index, size_t &row, size_t &col);
    void positionToCellIndex(const Position &position, Index &index);
    void cellIndexToPositionX(Index index, Position &p);
    void cellIndexToPositionY(Index index, Position &p);
    void cellIndexToPosition(Index index, Position &p);
    bool outOfBounds(const Index &index);
    std::map<size_t, std::string> getNeighbourType(const Index &index);
    size_t getNumNonLiquidNeighbours(const Index &index);
    size_t getNumNonSolidNeighbours(const Index &index);
    std::vector<Index> getNeighbourIndices(const Index &index);
    float laplacian(Index index, float time, Dimension dimension);

    bool bordersSolidCellX(const Index &index);
    bool bordersSolidCellY(const Index &index);
    bool bordersFluidCellX(const Index &index);
    bool bordersFluidCellY(const Index &index);

    size_t numFluidCells();
    float calculateTimeStep();
    bool isOutsideFluid(const Position &p);
    void updateGrid();
    float interpolate(const Position p, Dimension dimension);

    std::vector<std::vector<float>> m_x;
    std::vector<std::vector<float>> m_y;
    std::vector<std::vector<float>> m_pressure;
    std::vector<std::vector<float>> m_density;

    std::vector<std::vector<std::string>> m_type;
    std::vector<std::vector<size_t>> m_numParticles;
    std::vector<Position> m_particles;
    std::vector<std::vector<int>> m_indices;
    size_t m_resolution;
    float gridWidth = 1;
    float cellWidth;
    std::unique_ptr<ngl::AbstractVAO> m_vao;
    std::vector<ngl::Vec3> m_vbo;

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
