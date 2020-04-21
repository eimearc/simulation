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

private:
    // Velocity Methods
    ngl::Vec2 velocityAt(const float x, const float y);
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
    size_t getType(size_t row, size_t col);
    bool isFluidCell(size_t row, size_t col);
    size_t index(size_t row, size_t col);
    void coordinate(size_t index, size_t &row, size_t &col);
    void positionToCellIndex(float x, float y, size_t &row, size_t &col);
    bool outOfBounds(size_t row, size_t col);
    std::map<size_t, size_t> getNeighbours(size_t row, size_t col);
    size_t getNumNonLiquidNeighbours(size_t row, size_t col);
    std::vector<std::pair<size_t, size_t>> getNeighbourIndices(size_t row, size_t col);

    std::vector<std::vector<float>> m_x;
    std::vector<std::vector<float>> m_y;
    std::vector<std::vector<std::string>> m_type;
    std::vector<std::vector<size_t>> m_numParticles;
    std::vector<ngl::Vec2> m_particles;
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
    friend std::ostream& operator<<(std::ostream& os, MAC& mac);
};

std::ostream& operator<<(std::ostream& os, MAC& mac);
