// --- External Includes ---
#include <Eigen/Dense>
// --- Internal Includes ---
#include "MinimumDisc.hpp"

// --- STL Includes ---
#include <array>
#include <limits>
#include <stdexcept>
#include <vector>
#include <optional>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <random>


//// Create random generator
std::random_device randomDevice;
std::mt19937 randomGenerator(randomDevice());


namespace cie::geo {


// Type aliases -----------------------------------------------------------------
template <unsigned ArraySize>
using DoubleArray           = std::array<double,ArraySize>;

template <unsigned ArraySize>
using IntArray              = std::array<int, ArraySize>;

using IntVector             = std::vector<int>;
using UIntVector            = std::vector<unsigned>;

using PointSet2D            = std::vector<DoubleArray<2>>;
using PointSet2DPtr         = std::shared_ptr<PointSet2D>;
using ConstPointSet2DPtr    = std::shared_ptr<const PointSet2D>;


// Quick norm definition --------------------------------------------------------
double distance(const DoubleArray<2>& point1, const DoubleArray<2>& point2);


// Quick disc definition --------------------------------------------------------
struct Disc {
    Disc(const DoubleArray<2>& center, double radius);
    Disc(const DoubleArray<2>& point1, const DoubleArray<2>& point2);
    Disc(const DoubleArray<2>& point1, const DoubleArray<2>& point2, const DoubleArray<2>& point3);
    bool isInside(const DoubleArray<2>& point);

    double _radius2;
    DoubleArray<2> _center;
};


// Solver ----------------------------------------------------------------------
// Define class for encapsulating everything related to computing
// the minimum enclosing disc of a 2D point set
class MinimumEnclosingDisc {
public:
    MinimumEnclosingDisc();
    MinimumEnclosingDisc(const PointSet2D& pointSet);
    MinimumEnclosingDisc(const PointSet2DPtr& rp_points);
    int                build(double tolerance=1e-10);
    double          getRadius() const;
    DoubleArray<2>  getCenter() const;
    IntVector        getActiveIndices() const;

private:
    const DoubleArray<2>&   getPoint(size_t index) const;
    int                        checkEnclosure(double tolerance) const;
    int                        addActiveIndex(int index);
    void                    removeActiveIndex(int index);

    Disc                _disc;
    ConstPointSet2DPtr  _pointSetPtr;
    IntArray<4>         _activeIndices;
    UIntVector          _map;
};


// Define the necessary permutations to check, possible scenarios (neglecting 0 radius discs):
//    -    all possible discs defined by a set of 2 points:    1 possibility
//    -    all possible discs defined by a set of 3 points:    3 possibilities
//    -    all possible discs defined by a set of 4 points:    6 possibilities
std::vector<std::vector<IntVector>> permutationSets =
{
    {
        {0,1}
    },
    {
        {0,2},{1,2},
        {0,1,2}
    },
    {
        {0,3},{1,3},{2,3},
        {0,1,3},{0,2,3},{1,2,3}
    }
};


// Specialize and rename linalg::norm2 for easier use
double squareNorm(const DoubleArray<2>& vector)
{
    return vector[0] * vector[0] + vector[1] * vector[1];
}


// Array operators

DoubleArray<2> operator+(const DoubleArray<2>& r_lhs, const DoubleArray<2>& r_rhs)
{
    return DoubleArray<2> {r_lhs[0] + r_rhs[0], r_lhs[1] + r_rhs[1]};
}


DoubleArray<2> operator-(const DoubleArray<2>& r_lhs, const DoubleArray<2>& r_rhs)
{
    return DoubleArray<2> {r_lhs[0] - r_rhs[0], r_lhs[1] - r_rhs[1]};
}


DoubleArray<2> operator/(const DoubleArray<2>& r_lhs, double rhs)
{
    return DoubleArray<2> {r_lhs[0] / rhs, r_lhs[1] / rhs};
}


// Define the distance of two points - custom norm, no square root is needed
double distance(const DoubleArray<2>& point1, const DoubleArray<2>& point2)
{
    return squareNorm(point1 - point2);
}



Disc::Disc(const DoubleArray<2>& center, double radius) :
    _radius2(radius*radius),
    _center(center)
{
}


Disc::Disc( const DoubleArray<2>& point1,
            const DoubleArray<2>& point2) :
    _radius2(distance(point1, point2) / 4.0),
    _center((point1 + point2) / 2.0)
{
}


Disc::Disc( const DoubleArray<2>& point1,
            const DoubleArray<2>& point2,
            const DoubleArray<2>& point3 )
{
    // Check whether the definition is unique
    double tolerance    = 1e-15;
    double determinant  =   (point1[0]-point2[0])*(point1[1]-point3[1])
                            - (point1[0]-point3[0])*(point1[1]-point2[1]);

    if (std::abs(determinant)<tolerance)
    {
        // Definition is not unique, either 2 or more points are the same,
        // or the 3 points are co-linear

        double d12  = distance(point1,point2);
        double d13  = distance(point1,point3);
        double d23  = distance(point2,point3);

        Disc solution({0.0,0.0},0.0);

        if (d12<tolerance)
        {
            // Point1 and point2 are equal
            solution = Disc(point1,point3);
        }
        else if (d13<tolerance)
        {
            // Point1 and point3 are equal
            solution = Disc(point1,point2);
        }
        else if (d23<tolerance)
        {
            // Point2 and point3 are equal
            solution = Disc(point1,point2);
        }
        else
        {
            // The 3 points are co-linear
            // Suppose point1 and point2 are the defining points
            if ( d12>d23 && d12>d13 )
            {
                solution = Disc(point1,point2);
            }
            else if ( d13>d12 && d13>d23 )
            {
                // Nope, it was point1 and point3
                solution = Disc(point1,point3);
            }
            else
            {
                // Nope, it was point2 and point3
                solution = Disc(point2,point3);
            }

        }

        _center     = solution._center;
        _radius2    = solution._radius2;
        return;
    }

    // The disc definition is unique

    // Circle center defined by 3 points on its arc:
    //      --            -- ^(-1)    --                         --
    //   2* | x1-x2  y1-y2 |        * | (x1^2+y1^2) - (x2^2+y2^2) |
    //      | x1-x3  y1-y3 |          | (x1^2+y1^2) - (x3^2+y3^2) |
    //      --            --          --                         --

    // Calculate center
    Eigen::Matrix2d matrix;
    matrix << point1[1]-point3[1], point2[1]-point1[1],
              point3[0]-point1[0], point1[0]-point2[0];

    Eigen::Vector2d vector;
    vector << squareNorm(point1) - squareNorm(point2),
              squareNorm(point1) - squareNorm(point3);

    Eigen::Vector2d center = matrix * vector / determinant / 2.0;
    _center[0] = center[0];
    _center[1] = center[1];

    // Calculate radius based on the center
    _radius2 = distance(_center, point1);

    /*
    // DEBUG: check if the other two points are as far from the center as the first one
    double tolerance = 1e-16;
    if ( std::abs(distance(_center,point2)-_radius2)>tolerance || std::abs(distance(_center,point3)-_radius2)>tolerance )
        throw std::runtime_error("Disc center calculation failed!");
    */
}


bool Disc::isInside(const DoubleArray<2>& point)
{
    return distance(_center, point)<_radius2 ? true : false;
}


MinimumEnclosingDisc::MinimumEnclosingDisc()
    : MinimumEnclosingDisc(std::make_shared<PointSet2D>())
{
}


MinimumEnclosingDisc::MinimumEnclosingDisc(const PointSet2DPtr& rp_points)
    : _disc({0.0, 0.0}, 0.0),
      _pointSetPtr(rp_points),
      _activeIndices({-1,-1,-1,-1}),
      _map(rp_points->size())
{
    // Initialize map and shuffle it:
    //      from this point on, accessing the points in _pointSetPtr
    //      will be rerouted through _map. No element int pointSetPtr
    //      should be accessed with operator[], the private member
    //      function getPoint will be used instead.
    std::iota(_map.begin(), _map.end(), 0);
    std::shuffle(_map.begin(), _map.end(), randomGenerator); // <== disabled because of WASM

    // Fill the active set with the first point
    _disc._radius2        = 0.0;

    if (_pointSetPtr->size())
    {
        _activeIndices[0]   = 0;
        _disc._center        = getPoint(0);
    }
}



MinimumEnclosingDisc::MinimumEnclosingDisc(const PointSet2D& pointSet)
    : MinimumEnclosingDisc(std::make_shared<PointSet2D>(pointSet))
{
}


int MinimumEnclosingDisc::build(double tolerance)
{
    int outlierIndex;
    int activeSize;
    int restartCount = 0;

    for (restartCount=0;restartCount<(int)_pointSetPtr->size();++restartCount)
    {
        // Check whether all points are inside the current disc
        outlierIndex = checkEnclosure(tolerance);
        if (outlierIndex==-1)
            // Check succeeded, end the loop
            break;

        // Check failed, outlier identified
        // Add outlier ti the list of active indices
        activeSize = addActiveIndex(outlierIndex);

        // Find the points which defines a circle that encloses all other points
        // and pick the one with the least radius
        // 1) Define containers
        Disc trialDisc({0.0,0.0},0.0);
        IntArray<4> activeIndices;
        bool enclosed = true;
        std::pair<std::vector<Disc>,std::vector<IntVector>> candidates;

        // 2) Loop through possible permutations of the points
        //        (see the definition of permutationSets - the new point must be on the new circle)
        for (uint8_t permutation = 0; permutation < permutationSets[activeSize-2].size(); ++permutation)
            {
                auto indices = permutationSets[activeSize-2][permutation];

                // Define disc
                if (indices.size() == 2)
                {
                    trialDisc = Disc(    getPoint(_activeIndices[indices[0]]),
                                        getPoint(_activeIndices[indices[1]])    );
                }
                else if (indices.size() == 3)
                {
                    trialDisc = Disc(    getPoint(_activeIndices[indices[0]]),
                                        getPoint(_activeIndices[indices[1]]),
                                        getPoint(_activeIndices[indices[2]])    );
                }
                else
                    throw std::runtime_error("Invalid number of disc points!");

                // Check disc
                enclosed = true;
                for (uint8_t pIndex = 0; pIndex < activeSize; ++pIndex)
                {
                    if (    distance(trialDisc._center, getPoint(_activeIndices[pIndex])) > trialDisc._radius2+tolerance    )
                    {
                        enclosed = false;
                        break;
                    }
                }

                // If the current trialDisc encloses all other points from the active set,
                // store its data and the indices of the points in _activeIndices (yes, that's double index referencing)
                if (enclosed)
                {
                    candidates.first.push_back(trialDisc);
                    candidates.second.push_back(indices);
                }

        } // for permutation

        // Dump current state and throw error if no enclosing disc was found
        if (candidates.first.empty())
        {
            std::cout << "Dumping MinimumEnclosingDisc data:\n";
            std::cout << "\tActive set (" << activeSize << "):\n";
            for (uint8_t k = 0; k < _activeIndices.size(); ++k)
            {
                if (_activeIndices[k] == -1)
                    std::cout << "\t\tempty\n";
                else
                    std::cout << "\t\t" << getPoint(_activeIndices[k])[0] << ", " << getPoint(_activeIndices[k])[1] << "\n";
            }
            throw std::runtime_error("Failed to find a new disc!");
        }

        // Find best disc if there are candidates
        auto minDisc = std::min_element(    candidates.first.begin(),
                                            candidates.first.end(),
                                            [this](const Disc& lhs, const Disc& rhs) ->bool
                                                {return (lhs._radius2 < rhs._radius2 && lhs._radius2 >= this->_disc._radius2) ? true : false; });
        _disc = *minDisc;

        // Remove points from the active set that are enclosed (no longer active)
        activeIndices = _activeIndices;
        int discIndex = std::distance(candidates.first.begin(), minDisc);
        for (uint8_t k = 0; k < activeSize; ++k)
        {
            if (std::find(candidates.second[discIndex].begin(), candidates.second[discIndex].end(), k) == candidates.second[discIndex].end())
            {
                removeActiveIndex(activeIndices[k]);
            }
        }

    } // while

    return restartCount;
}


double MinimumEnclosingDisc::getRadius() const
{
    return std::sqrt(_disc._radius2);
}


DoubleArray<2> MinimumEnclosingDisc::getCenter() const
{
    return _disc._center;
}


IntVector MinimumEnclosingDisc::getActiveIndices() const
{
    // Returns the unshuffled indices of the active set
    IntVector indices(_activeIndices.size());
    uint8_t count = 0;
    for (uint8_t i = 0; i < _activeIndices.size(); ++i)
    {
        if (_activeIndices[i] != -1)
        {
            indices[count++] = (int)_map[_activeIndices[i]];
        }
    }
    indices.resize(count);
    return indices;
}


const DoubleArray<2>& MinimumEnclosingDisc::getPoint(size_t index) const
{
    // Reroute container indexing through the member variable: _map
    return (*_pointSetPtr)[_map[index]];
}


int MinimumEnclosingDisc::checkEnclosure(double tolerance) const
{
    int outlierIndex = -1;
    for (size_t i=0; i<_pointSetPtr->size();++i)
    {
        if (distance(_disc._center, getPoint(i)) > _disc._radius2+tolerance)
        {
            outlierIndex = (int)i;
            break;
        }
    }
    return outlierIndex;
}


int MinimumEnclosingDisc::addActiveIndex(int index)
{
    // Copy index to an empty position in _activeIndices, return the
    // number of non-empty active indices
    // (empty entries have a value of -1)
    size_t i = 0;
    for (i = 0; i < _activeIndices.size(); ++i)
    {
        if (_activeIndices[i] < 0)
        {
            _activeIndices[i] = index;
            break;
        }
        else if (i == _activeIndices.size() - 1)
        {
            std::cout << "Dumping MinimumEnclosingDisc data:\n"
                << "Active indices: "
                << _activeIndices[0] << ", "
                << _activeIndices[1] << ", "
                << _activeIndices[2] << ", "
                << _activeIndices[3] << "\n";
            std::cout << "Current disc:\n"
                << "\tcenter: " << _disc._center[0] << ", " << _disc._center[1] << "\n"
                << "\tradius: " << std::sqrt(_disc._radius2) << "\n";
            throw std::runtime_error("");
        }
    }

    return (int)++i;
}


void MinimumEnclosingDisc::removeActiveIndex(int index)
{
    // Set the entry of _activeIndices containing index to -1,
    // then compress
    bool swap = false;
    int temp;
    for (size_t i = 0; i < _activeIndices.size(); ++i)
    {
        if (!swap && _activeIndices[i] == index)
            swap = true;
        else if (!swap && i == _activeIndices.size() - 1)
        {
            std::cout << "Dumping MinimumEnclosingDisc data:\n"
                << "Active indices: "
                << _activeIndices[0] << ", "
                << _activeIndices[1] << ", "
                << _activeIndices[2] << ", "
                << _activeIndices[3] << "\n";
            std::cout << "Current disc:\n"
                << "\tcenter: " << _disc._center[0] << ", " << _disc._center[1] << "\n"
                << "\tradius: " << std::sqrt(_disc._radius2) << "\n";
            throw std::runtime_error("");
        }

        if (swap)
        {
            if (i < _activeIndices.size() - 1)
            {
                temp                    = _activeIndices[i];
                _activeIndices[i]        = _activeIndices[i + 1];
                _activeIndices[i + 1]    = temp;
            }
            else
            {
                _activeIndices[i] = -1;
            }
        }
    }
}



} // namespace cie::geo


namespace md {


struct MinimumDisc::Impl
{
    Impl()
        : _p_points(new cie::geo::PointSet2D),
          _disc()
    {}

    cie::geo::PointSet2DPtr _p_points;

    cie::geo::MinimumEnclosingDisc _disc;
}; // struct MinimumDisc::Impl


MinimumDisc::MinimumDisc()
    : _p_impl(new Impl)
{
}


MinimumDisc::~MinimumDisc()
{
}


void MinimumDisc::include(Reference<const Point> r_point)
{
    _p_impl->_p_points->push_back(r_point);
    _p_impl->_disc = cie::geo::MinimumEnclosingDisc(cie::geo::PointSet2DPtr {_p_impl->_p_points});
    _p_impl->_disc.build();
}


MinimumDisc::Point MinimumDisc::getCenter() const
{
    return _p_impl->_disc.getCenter();
}


MinimumDisc::Numeric MinimumDisc::getRadius() const
{
    return _p_impl->_disc.getRadius();
}


} // namespace md
