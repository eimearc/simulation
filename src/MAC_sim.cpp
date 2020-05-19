#include "MAC.h"

#include <gflags/gflags.h>

// =====================
// Simulation Methods
// =====================

DECLARE_int32(resolution);
DECLARE_int32(num_particles);
DECLARE_double(viscosity);
DECLARE_int32(obstacles);
DECLARE_double(time_step);

void MAC::updateVectorField()
{
    const float time_step = float(FLAGS_time_step);
    static float timeElapsed = 0.0f;
    float time = calculateTimeStep();
    if ((timeElapsed + time) > time_step)
    {
        time = time_step-timeElapsed;
    }
    timeElapsed+=time;
    if(timeElapsed>=time_step)
    {
        timeElapsed=0.0f;
        m_frame=true;
    }
    updateGrid();
    applyConvection(time);
    applyExternalForces(time);
    applyViscosity(time);
    calculatePressure(time);
    applyPressure(time);
    fixBorderVelocities();
    moveParticles(time);
}

float MAC::calculateTimeStep()
{
    float maxU = 0.000001f;
    for (int row = 0; row <= m_resolution; ++row)
    {
        for (int col = 0; col <= m_resolution; ++col)
        {
            if (row < m_resolution)
            {
                if ((m_x[row][col])*(m_x[row][col]) > maxU) maxU = (m_x[row][col])*(m_x[row][col]);
            }
            if (col < m_resolution)
            {
                if ((m_y[row][col])*(m_y[row][col]) > maxU) maxU = (m_y[row][col])*(m_y[row][col]);
            }
        }
    }
    return cellWidth/sqrt(maxU);
}

void MAC::updateGrid()
{
    Index index;
    for (index.row = 0; index.row < m_resolution-1; ++index.row)
    {
        for (index.col = 0; index.col < m_resolution-1; ++index.col)
        {
            if (!isSolidCell(index)) m_type[index.row][index.col] = AIR;
        }
    }

    for (const auto &p: m_particles)
    {
        Index index;
        positionToCellIndex(p, index);
        if (!outOfBounds(index) && !isSolidCell(index)) m_type[index.row][index.col] = FLUID;
    }
}

void MAC::applyConvection(float _time)
{
    Position p;
    Index index;
    Velocity updated;
    MAC tmp(m_resolution);
    for (index.row = 0; index.row <= m_resolution; ++index.row)
    {
        for (index.col = 0; index.col <= m_resolution; ++index.col)
        {
            size_t row = index.row;
            size_t col = index.col;
            if (bordersFluidCellX(index))
            {
                cellIndexToPositionX(index, p);
                updated = traceParticle(p, _time);
                tmp.m_x[row][col] = updated.m_x;
            }
            if (bordersFluidCellY(index))
            {
                cellIndexToPositionY(index, p);
                updated = traceParticle(p, _time);
                tmp.m_y[row][col] = updated.m_y;
            }
        }
    }
    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

void MAC::applyExternalForces(float _time)
{
    Velocity gravityVector = {0, -9.80665};
    gravityVector*=_time;
    Index index;
    for (int col = 0; col <= m_resolution; ++col)
    {
        index.col = col;
        for (int row = 0; row <= m_resolution; ++row)
        {
            index.row = row;
            if (bordersFluidCellY(index)) m_y[row][col] += gravityVector.m_y;
            if (bordersFluidCellX(index)) m_x[row][col] += gravityVector.m_x;
        }
    }
}

void MAC::applyViscosity(float _time)
{
    MAC tmp(m_resolution);
    tmp.m_x = m_x;
    tmp.m_y = m_y;

    Index index;
    for (index.row=0;index.row<=int(m_resolution);++index.row)
    {
        for (index.col=0;index.col<=int(m_resolution);++index.col)
        {
            if (bordersFluidCellX(index))
            {
                float l = laplacian(index, _time, Dimension::x);
                tmp.m_x[index.row][index.col] += l;
            }
            if (bordersFluidCellY(index))
            {
                float l = laplacian(index, _time, Dimension::y);
                tmp.m_y[index.row][index.col] += l;
            }
        }
    }

    m_x=tmp.m_x;
    m_y=tmp.m_y;
}

void MAC::calculatePressure(float _time)
{
    m_indices = std::vector<std::vector<int>>(m_resolution, std::vector<int>(m_resolution, -1));
    size_t i = 0;
    Index index;
    for (int row = 0; row < m_resolution ; ++row)
    {
        index.row=row;
        for (int col = 0; col < m_resolution; ++col)
        {
            index.col=col;
            if (isFluidCell(index))
            {
                m_indices[row][col] = i;
                ++i;
            }
        }
    }

    auto A = constructCoefficientMatrix();
    auto b = constructDivergenceVector(_time);

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd p = solver.solve(b);

    m_pressure = std::vector<std::vector<double>>(m_resolution, std::vector<double>(m_resolution, ATMOSPHERIC_PRESSURE));
    for (int row = 0; row < m_resolution; ++row)
    {
        index.row=row;
        for (int col = 0; col < m_resolution; ++col)
        {
            index.col=col;
            if (isFluidCell(index))
            {
                size_t j = m_indices[row][col];
                m_pressure[row][col] = p[j];
            }
        }
    }

    // Fix pressure for solid cells.
    for (int row = 0; row < m_resolution; ++row)
    {
        m_pressure[row][0] = m_pressure[row][1];
        m_pressure[row][m_resolution-1] = m_pressure[row][m_resolution-2];
    }
    for (int col = 0; col < m_resolution; ++col)
    {
        m_pressure[0][col] = m_pressure[1][col];
        m_pressure[m_resolution-1][col] = m_pressure[m_resolution-2][col];
    }
}

void MAC::applyPressure(float _time)
{
    MAC tmp(m_resolution);
    tmp.m_x = m_x;
    tmp.m_y = m_y;

    Index index;
    for (index.row = 0; index.row <= int(m_resolution); ++index.row)
    {
        for (index.col = 0; index.col <= int(m_resolution); ++index.col)
        {
            if (bordersFluidCellX(index) && !(bordersSolidCellX(index)))
            {
                Velocity v = applyPressureToPoint(index,_time,Dimension::x);
                tmp.m_x[index.row][index.col] = v.m_x;
            }

            if (bordersFluidCellY(index) && !(bordersSolidCellY(index)))
            {
                Velocity v = applyPressureToPoint(index,_time,Dimension::y);
                tmp.m_y[index.row][index.col] = v.m_y;
            }
        }
    }

    m_x = tmp.m_x;
    m_y = tmp.m_y;
}

void MAC::fixBorderVelocities()
{
    Index index;
    for (int row = 0; row < m_resolution; ++row)
    {
        for (int col = 0; col < m_resolution; ++col)
        {
            index.row=row;
            index.col=col;
            if(bordersSolidCellX(index))
            {
                if(bordersFluidCellX({index.row, index.col-1}))
                {
                    if(m_x[row][col]>0) m_x[row][col]=0;
                }
                else if(bordersFluidCellX({index.row, index.col+1}))
                {
                    if(m_x[row][col]<0) m_x[row][col]=0;
                }
                else
                {
                    m_x[row][col] = 0.0f;
                }
            }
            if(bordersSolidCellY(index))
            {
                if(bordersFluidCellY({index.row-1, index.col}))
                {
                    if(m_y[row][col]>0) m_y[row][col]=0;
                }
                if(bordersFluidCellY({index.row+1, index.col}))
                {
                    if(m_y[row][col]<0) m_y[row][col]=0;
                }
                else
                {
                    m_y[row][col] = 0.0f;
                }
            }
        }
    }
    // Top row y.
    for (float &v : m_y[m_resolution])
    {
        v = 0.0f;
    }
    // One from top row y.
    for (float &v : m_y[m_resolution-1])
    {
        if (v > 0) v = 0.0f;
    }
    // Bottom row y.
    for (float &v : m_y[0])
    {
        v = 0.0f;
    }
    // One from bottom row y.
    for (float &v : m_y[1])
    {
        if (v < 0) v = 0.0f;
    }
    // Right and left col y.
    for (auto &col : m_y)
    {
        col[0] = 0.0f;
        col[m_resolution-1] = 0.0f;
    }
    // Top row x.
    for (auto &v : m_x[m_resolution-1])
    {
        v = 0.0f;
    }
    // Bottom row x.
    for (auto &v : m_x[0])
    {
        v = 0.0f;
    }
    // Left and right solids x.
    for (auto &col : m_x)
    {
        col[m_resolution] = 0.0f;
        col[0] = 0.0f;
        if (col[1] < 0) col[1] = 0.0f;
        if (col[m_resolution-1] > 0) col[m_resolution-1] = 0.0f;
    }
}

void MAC::moveParticles(float _time)
{
    for (Position &p : m_particles)
    {
        Velocity velocity = traceParticle(p,_time);
        p += _time*velocity;
    }
}

// =========================
// Viscosity Helper Methods
// =========================

float MAC::laplacian(Index index, float time, Dimension dimension)
{
    float l = 0.0f;
    std::vector<std::vector<float>> *pm;
    std::function<bool(Index)> bordersFluidCell;

    switch (dimension)
    {
    case Dimension::x :
        pm = &m_x;
        bordersFluidCell = [&](Index i)
        {
            return this->bordersFluidCellX(i);
        };
        break;
    case Dimension::y :
        pm = &m_y;
        bordersFluidCell = [&](Index i)
        {
            return this->bordersFluidCellY(i);
        };
    }
    const std::vector<std::vector<float>> &m = *pm;

    float x1=0.0f,x2=0.0f;
    float y1=0.0f,y2=0.0f;

    int row=index.row;
    int col=index.col;
    if (bordersFluidCell({row, col-1})) x1=m[row][col-1];
    if (bordersFluidCell({row, col+1})) x2=m[row][col+1];
    if (bordersFluidCell({row-1, col})) y1=m[row-1][col];
    if (bordersFluidCell({row+1, col})) y2=m[row+1][col];

    l = x1 + x2 + y1 + y2;

    l = time*FLAGS_viscosity*l;

    return l;
}

// =========================
// Velocity Helper Methods
// =========================
Velocity MAC::velocityAtPosition(const Position p)
{
    // Separately bilinearly interpolate x and y.
    Velocity v;
    v.m_x = interpolate(p,Dimension::x);
    v.m_y = interpolate(p,Dimension::y);

    if (isOutsideGrid(p)||isOutsideFluid(p)) v={0,0};
    return v;
}

float MAC::interpolate(const Position p, Dimension dimension)
{
    float result = 0.0f;
    Index index;
    positionToCellIndex(p,index);
    std::vector<std::vector<float>> *pm;
    Position c;
    switch(dimension)
    {
    case Dimension::x :
    {
        pm = &m_x;
        cellIndexToPositionX(index,c);
        break;
    }
    case Dimension::y :
    {
        pm = &m_y;
        cellIndexToPositionY(index,c);
        break;
    }
    }
    const auto &m = *pm;
    const Position cellCenter = c;

    float q1=0.0f, q2 = 0.0f, q3 = 0.0f, q4 = 0.0f;

    if (p.m_y<cellCenter.m_y)
    {
        index.row--;
    }
    if (p.m_x<cellCenter.m_x)
    {
        index.col--;
    }

    if (index.row<0) index.row = 0;
    if (index.row>=m_resolution)
    {
        switch(dimension)
        {
        case Dimension::x :
        {
            index.row=m_resolution-2;
            break;
        }
        case Dimension::y :
        {
            index.row=m_resolution-1;
            break;
        }
        }
    }
    if (index.col<0) index.col = 0;
    if (index.col>=m_resolution)
    {
        switch(dimension)
        {
        case Dimension::x :
        {
            index.col=m_resolution-1;
            break;
        }
        case Dimension::y :
        {
            index.col=m_resolution-2;
            break;
        }
        }
    }

    const int tmpRow = index.row;
    const int tmpCol = index.col;

    Position p1;
    Position p2;
    if (dimension == Dimension::x)
    {
        cellIndexToPositionX({int(tmpRow),int(tmpCol)}, p1);
        cellIndexToPositionX({tmpRow+1, tmpCol+1}, p2);
    }
    else if (dimension == Dimension::y)
    {
        cellIndexToPositionY({tmpRow, tmpCol}, p1);
        cellIndexToPositionY({tmpRow+1, tmpCol+1}, p2);
    }
    q1 = m[tmpRow][tmpCol];
    if (tmpCol < int(m[0].size()-1)) q2 = m[tmpRow][tmpCol+1];
    if (tmpRow < int(m.size()-1)) q3 = m[tmpRow+1][tmpCol];
    if (tmpCol < int(m[0].size()-1) && tmpRow < int(m.size()-1)) q4 = m[tmpRow+1][tmpCol+1];

    const float &x1=p1.m_x;
    const float &x2=p2.m_x;
    const float &y1=p1.m_y;
    const float &y2=p2.m_y;
    const float &x=p.m_x;
    const float &y=p.m_y;

    // X direction.
    auto fx1 = (x2-x)/(x2-x1)*q1 + (x-x1)/(x2-x1)*q2;
    auto fx2 = (x2-x)/(x2-x1)*q3 + (x-x1)/(x2-x1)*q4;
    // Y direction.
    result = (y2-y)/(y2-y1)*fx1 + (y-y1)/(y2-y1)*fx2;
    return result;
}

Velocity MAC::velocityAtIndex(const Index index)
{
    Position p;
    cellIndexToPosition(index, p);
    return velocityAtPosition(p);
}

Velocity MAC::traceParticle(const Position &p, float _time)
{
    // Trace particle from point (_x, _y) using RK2.
    Velocity half_prev_v = velocityAtPosition(
                p-(0.5f*_time*velocityAtPosition(p))
    );
    return velocityAtPosition(p-_time*(velocityAtPosition(p - _time*half_prev_v)));
}

// =========================
// Pressure Helper Methods
// =========================

Velocity MAC::applyPressureToPoint(const Index &index, float _time, Dimension dimension)
{
    Velocity v;
    Velocity gradient;
    float density=0.0f;
    const int &row = index.row;
    const int &col = index.col;

    if (dimension==Dimension::x)
    {
        v = {m_x[row][col],0};
        float x1 = (m_pressure[row][col]+m_pressure[row][col-1])/2.0f;
        float x2 = (m_pressure[row][col-1]+m_pressure[row][col-2])/2.0f;
        gradient = {x1-x2,0.0f}; // Backwards difference.
        density=(m_density[row][col-1] + m_density[row][col])/2.0f;
    }
    else if (dimension==Dimension::y)
    {
        v = {0,m_y[row][col]};
        float y1 = (m_pressure[row][col]+m_pressure[row-1][col])/2.0f;
        float y2 = (m_pressure[row-1][col]+m_pressure[row-2][col])/2.0f;
        gradient = {0.0f,y1-y2}; // Backwards difference.
        density=(m_density[row-1][col] + m_density[row][col])/2.0f;
    }

    auto rhs = (_time/(density*cellWidth))*gradient;
    Velocity result = v-rhs;
    return result;
}

Eigen::SparseMatrix<double> MAC::constructCoefficientMatrix()
{
    size_t n = numFluidCells();
    Eigen::SparseMatrix<double> m(n,n);
    auto tripletList = constructNeighbourTriplets();
    m.setFromTriplets(tripletList.begin(), tripletList.end());
    return m;
}

std::vector<Eigen::Triplet<double>> MAC::constructNeighbourTriplets()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::Triplet<double> t;
    size_t i;
    Index index;
    for (int col = 0; col < int(m_resolution); ++col)
    {
        index.col = col;
        for (int row = 0; row < int(m_resolution); ++row)
        {
            index.row = row;
            if(isFluidCell(index))
            {
                i = m_indices[row][col];
                t = Eigen::Triplet<double>(i,i,-1*int(getNumNonSolidNeighbours(index)));
                tripletList.push_back(t);

                auto neighbours = getNeighbourIndices(index);
                for ( const auto &neighbour : neighbours)
                {
                    const int &n_row = neighbour.row;
                    const int &n_col = neighbour.col;
                    if(isFluidCell(neighbour))
                    {
                        size_t j = m_indices[n_row][n_col];
                        t = Eigen::Triplet<double>(i,j,1);
                        tripletList.push_back(t);
                    }
                }
            }
        }
    }
    return tripletList;
}

Eigen::VectorXd MAC::constructDivergenceVector(float _time)
{
    size_t n = numFluidCells();
    float maxParticlesPerCell = FLAGS_num_particles/100.0f;

    Eigen::VectorXd v(n);
    v.setZero();
    std::vector<std::vector<size_t>> numParticles(m_resolution, std::vector<size_t>(m_resolution, 0));
    Index index;
    for (const auto &p : m_particles)
    {
        positionToCellIndex(p, index);
        if (!outOfBounds(index))
        {
            numParticles[index.row][index.col]++;
        }
    }
    m_numParticles = numParticles;
    for (int col = 0; col < m_resolution; ++col)
    {
        for (int row = 0; row < m_resolution; ++row)
        {
            index = {row,col};
            m_density[row][col] = AIR_DENSITY;
            if (isFluidCell(index))
            {
                float f = m_numParticles[row][col]/maxParticlesPerCell;
                m_density[row][col] = (WATER_DENSITY*f)+(AIR_DENSITY*std::max((1-f),0.0f));
            }
        }
    }

    // Fix density for edge cells.
    for (int row = 0; row < m_resolution; ++row)
    {
        m_density[row][0] = m_density[row][1];
        m_density[row][m_resolution-1] = m_density[row][m_resolution-2];
    }
    for (int col = 0; col < m_resolution; ++col)
    {
        m_density[0][col] = m_density[1][col];
        m_density[m_resolution-1][col] = m_density[m_resolution-2][col];
    }

    for (int row = 0; row < m_resolution; ++row)
    {
        for (int col = 0; col < m_resolution; ++col)
        {
            index.row=row;
            index.col=col;
            if (isFluidCell(index))
            {
                size_t i = m_indices[row][col];
                float density = m_density[row][col];
                const float h = cellWidth;
                float divergence = calculateModifiedDivergence(row,col);

                size_t numNeighbourAirCells = 0;
                if (isAirCell(row,col-1)) numNeighbourAirCells++;
                if (isAirCell(row,col+1)) numNeighbourAirCells++;
                if (isAirCell(row-1,col)) numNeighbourAirCells++;
                if (isAirCell(row+1,col)) numNeighbourAirCells++;

                auto result = ((density*h)/_time)*divergence - (numNeighbourAirCells*ATMOSPHERIC_PRESSURE);

//                printf("i: %d result: %f density: %f h: %f _time: %f divergence: %f neighbourAirCells: %d ATMOSPHERIC_PRESSURE: %d\n",
//                       int(i), result, density, h, _time, divergence, int(numNeighbourAirCells), int(ATMOSPHERIC_PRESSURE));

                v[i] = result;
            }
        }
    }

    return v;
}

float MAC::calculateModifiedDivergence(size_t _row, size_t _col)
{
    const int &row = _row;
    const int &col = _col;
    // Div(u)
    // Velocity components between fluid cells and solid cells
    // are considered to be zero.
    float x1=0.0f,y1=0.0f,x2=0.0f,y2=0.0f;
    Index index = {row,col};
    if (!bordersSolidCellY(index)) y1 = m_y[row][col];
    if (!bordersSolidCellX(index)) x1 = m_x[row][col];

    index={row,col+1};
    if (!bordersSolidCellX(index)) x2 = m_x[row][col+1];
    index = {row+1,col};
    if (!bordersSolidCellY(index)) y2 = m_y[row+1][col];

    float xDiv = x2-x1;
    float yDiv = y2-y1;

    return xDiv + yDiv;
}
