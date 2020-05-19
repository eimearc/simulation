#include "MAC.h"

// =====================
// Helper Methods
// =====================

void MAC::cellIndexToPositionX(Index index, Position &p)
{
    cellIndexToPosition(index, p);
    p.m_x -= 0.5*cellWidth;
}

void MAC::cellIndexToPositionY(Index index, Position &p)
{
    cellIndexToPosition(index, p);
    p.m_y -= 0.5*cellWidth;
}

void MAC::cellIndexToPosition(Index index, Position &p)
{
    float l = -gridWidth/2.0f;
    p.m_x = l + index.col*(cellWidth) + 0.5*(cellWidth);
    p.m_y = l + index.row*(cellWidth) + 0.5*(cellWidth);
}

bool MAC::isOutsideGrid(const Position &p)
{
    Index index;
    positionToCellIndex(p, index);
    const int &row = index.row;
    const int &col = index.col;
    if (row > m_resolution-1 || row < 1 || col > m_resolution-1 || col < 1)
    {
        return true;
    }
    return false;
}

bool MAC::isOutsideFluid(const Position &p)
{
    if (isOutsideGrid(p)) return true;
    Index index;
    positionToCellIndex(p, index);
    if (isSolidCell(index)) return true;
    return false;
}

void MAC::positionToCellIndex(const Position &position, Index &index)
{
    float x = position.m_x;
    float y = position.m_y;
    x += gridWidth/2.0f;
    y += gridWidth/2.0f;
    index.col = x/cellWidth;
    index.row = y/cellWidth;
}

Type MAC::getType(const Index &index)
{
    if (outOfBounds(index))
    {
        return AIR; // Return air (non-solid)
    }
    return m_type[index.row][index.col];
}

size_t MAC::getNumNonLiquidNeighbours(const Index &index)
{
    auto neighbours = getNeighbourIndices(index);
    size_t num = 0;
    for (const auto& neighbour : neighbours)
    {
        if (!isFluidCell(neighbour))
        {
            num++;
        }
    }
    return num;
}

size_t MAC::getNumNonSolidNeighbours(const Index &index)
{
    auto neighbours = getNeighbourIndices(index);
    size_t num = 0;
    for (const auto& neighbour : neighbours)
    {
        if (!isSolidCell(neighbour))
        {
            num++;
        }
    }
    return num;
}

std::vector<Index> MAC::getNeighbourIndices(const Index &index)
{
    std::vector<Index> indices;
    const int &row = index.row;
    const int &col = index.col;

    if (row < int(m_resolution-1)) indices.push_back({row+1, col});
    if (row > 0) indices.push_back({row-1, col});
    if (col < int(m_resolution-1)) indices.push_back({row, col+1});
    if (col > 0) indices.push_back({row, col-1});

    return indices;
}

void MAC::coordinate(size_t index, size_t &row, size_t &col)
{
    col = index%m_resolution;
    row = index/m_resolution;
}

bool MAC::outOfBounds(const Index &index)
{
    const int &row = index.row;
    const int &col = index.col;
    const int &resolution = m_resolution;
    return ((row >= resolution) || (col >= resolution) || (row < 0) || (col < 0));
}

bool MAC::isFluidCell(const Index &index)
{
    if (outOfBounds(index))
    {
        return false;
    }
    return m_type[index.row][index.col] == FLUID;
}

bool MAC::isSolidCell(const Index &index)
{
    if (outOfBounds(index))
    {
        return false;
    }
    return m_type[index.row][index.col] == SOLID;
}

bool MAC::isAirCell(size_t row, size_t col)
{
    if (outOfBounds({int(row), int(col)}))
    {
        return false;
    }
    return m_type[row][col] == AIR;
}

size_t MAC::numFluidCells()
{
    size_t num = 0;
    for (int row = 0; row <m_resolution; ++row)
    {
        for (int col = 0; col < m_resolution; ++col)
        {
            if (isFluidCell({int(row), int(col)}))
            {
                num++;
            }
        }
    }
    return num;
}

bool MAC::bordersSolidCellX(const Index &_index)
{
    Index index = _index;
    if (outOfBounds(index)) return false;
    if (isSolidCell(index)) return true;
    index.col--;
    if(outOfBounds(index)) return false;
    if (isSolidCell(index)) return true;
    return false;
}

bool MAC::bordersSolidCellY(const Index &_index)
{
    Index index = _index;
    if (outOfBounds(index)) return false;
    if (isSolidCell(index)) return true;
    index.row--;
    if(outOfBounds(index)) return false;
    if (isSolidCell(index)) return true;
    return false;
}

bool MAC::bordersFluidCellX(const Index &_index)
{
    Index index = _index;
    if (outOfBounds(index)) return false;
    if (isFluidCell(index)) return true;
    index.col--;
    if(outOfBounds(index)) return false;
    if (isFluidCell(index)) return true;
    return false;
}

bool MAC::bordersFluidCellY(const Index &_index)
{
    Index index = _index;
    if (outOfBounds(index)) return false;
    if (isFluidCell(index)) return true;
    index.row--;
    if(outOfBounds(index)) return false;
    if (isFluidCell(index)) return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, MAC& mac)
{
    for (int i = mac.m_type.size()-1; i >= 0; --i)
    {
        os << i << ":\t";
        for (const auto &y: mac.m_type[i])
        {
            os << y << ' ';
        }
        os << '\n';
    }
    os << '\n';

    std::cout << "pressure\n" << mac.m_pressure << std::endl;

    std::cout << "density\n" << mac.m_density << std::endl;

    os << std::fixed << std::setprecision(4) << std::setfill('0') << std::showpos;
    for (int i = mac.m_y.size()-1; i >= 0 ; --i)
    {
        if (i < int(mac.m_x.size()))
        {
            os << "X" << i << "  ";
            for (const auto &x : mac.m_x[i])
            {
                os << x;
                os << "     ";
            }
        }
        os << "\n\n";
        os << "Y" << i << "  ";
        for (const auto &y : mac.m_y[i])
        {
            os <<"  *  ";
            os << y;
        }
        os <<"  *  ";

        os << "\n\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<float>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<double>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<size_t>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<std::vector<int>>& grid)
{
    os << std::fixed << std::setprecision(4) << std::setfill('0');
    for (int row = grid.size()-1; row >= 0 ; --row)
    {
        for (const auto &e : grid[row])
        {
            os << e;
            os << "    ";
        }

        os << "\n\n";
    }
    return os;
}
