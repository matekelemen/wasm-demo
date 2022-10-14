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


namespace md {


double distanceSquared(const md::MinimumDisc::Point& r_lhs, const md::MinimumDisc::Point& r_rhs = {0.0, 0.0})
{
    double output = 0.0;
    for (unsigned dim=0; dim<2; ++dim)
    {
        const double tmp = r_lhs[dim] - r_rhs[dim];
        output += tmp * tmp;
    }
    return output;
}


class Disc
{
public:
    using Point = md::MinimumDisc::Point;

    using Numeric = md::MinimumDisc::Numeric;

public:
    Disc()
        : _center({0.0, 0.0}),
          _radiusSquared(0.0)
    {}

    Disc(const Point& r_point)
        : _center(r_point),
          _radiusSquared(0.0)
    {}

    Disc(const Point& r_point0, const Point& r_point1)
        : _center(),
          _radiusSquared(distanceSquared(r_point0, r_point1) / 4.0)
    {
        for (unsigned dim=0; dim<2; ++dim)
            _center[dim] = (r_point0[dim] + r_point1[dim]) / 2.0;
    }

    Disc(const Point& r_point0,
         const Point& r_point1,
         const Point& r_point2)
    {
        constexpr const double tolerance = 1e-14;
        const double determinant = (r_point0[0]-r_point1[0])*(r_point0[1]-r_point2[1])
                                   - (r_point0[0]-r_point2[0])*(r_point0[1]-r_point1[1]);

        // The definition is well formed.
        if (tolerance < std::abs(determinant))
        {
            // Using eigen for the sake of using an external dependency.
            Eigen::Matrix2d matrix;
            matrix << r_point0[1]-r_point2[1], r_point1[1]-r_point0[1],
                      r_point2[0]-r_point0[0], r_point0[0]-r_point1[0];

            Eigen::Vector2d vector;
            vector << distanceSquared(r_point0) - distanceSquared(r_point1),
                      distanceSquared(r_point0) - distanceSquared(r_point2);

            Eigen::Vector2d center = matrix * vector / determinant / 2.0;
        }

        // The definition is ill-formed: either 2 points coincide or
        // the 3 points are co-linear.
        else
        {
            // Put the points into an array to make working with them easier.
            std::array<std::reference_wrapper<const Point>, 3> points {
                r_point0,
                r_point1,
                r_point2
            };

            // Possible pair combinations and their distances.
            std::array<std::pair<std::pair<unsigned, unsigned>, Numeric>, 3> pairs {{
                {{0u, 1u}, distanceSquared(r_point0, r_point1)},
                {{0u, 2u}, distanceSquared(r_point0, r_point2)},
                {{1u, 2u}, distanceSquared(r_point1, r_point2)}
            }};
            std::sort(pairs.begin(), pairs.end(), [](const auto& r_lhs, const auto& r_rhs) {return r_lhs.second < r_rhs.second;});

            *this = Disc(points[pairs[0].first.first], points[pairs[0].first.second]);
        }
    }

    Disc(Disc&& r_rhs) noexcept = default;

    Disc(const Disc& r_rhs) = default;

    Disc& operator=(Disc&& r_rhs) noexcept = default;

    Disc& operator=(const Disc& r_rhs) = default;

    bool isInside(const Point& r_point) const
    {
        return distanceSquared(r_point, _center) <= _radiusSquared;
    }

    Numeric getRadiusSquared() const noexcept
    {
        return _radiusSquared;
    }

    void setRadiusSquared(Numeric value)
    {
        _radiusSquared = value;
    }

    const Point& getCenter() const noexcept
    {
        return _center;
    }

private:
    Point _center;

    md::MinimumDisc::Numeric _radiusSquared;
}; // class Disc


struct MinimumDisc::Impl
{
    Impl()
        : _center({0.0, 0.0}),
          _radiusSquared(0.0)
    {}

    Impl(Impl&& r_rhs) noexcept = default;

    Impl(const Impl& r_rhs) = default;

    static std::vector<std::vector<std::vector<unsigned>>> permutationSets;

    std::array<std::optional<MinimumDisc::Point>, 3> _activeSet;

    MinimumDisc::Point _center;

    MinimumDisc::Numeric _radiusSquared;
}; // struct MinimumDisc::Impl


std::vector<std::vector<std::vector<unsigned>>> MinimumDisc::Impl::permutationSets
{
    {
        {0u, 1u}
    },
    {
        {0u, 2u}, {1u, 2u},
        {0u, 1u, 2u}
    },
    {
        {0u, 3u}, {1u, 3u}, {2u, 3u},
        {0u, 1u, 3u}, {0u, 2u, 3u}, {1u, 2u, 3u}
    }
};


MinimumDisc::MinimumDisc()
    : _p_impl(new Impl)
{
}


MinimumDisc::~MinimumDisc()
{
}


void MinimumDisc::include(Reference<const Point> r_point)
{
    if (!this->isIncluded(r_point))
    {
        // Define candidates - the new point must be in the active set.
        std::array<std::optional<std::reference_wrapper<const Point>>, 4> candidates
        {
            r_point,
            _p_impl->_activeSet[0],
            _p_impl->_activeSet[1],
            _p_impl->_activeSet[2]
        };

        std::reference_wrapper<const std::vector<unsigned>> newSet = Impl::permutationSets[0][0];
        Disc disc;
        disc.setRadiusSquared(1000);

        const unsigned activeSetSize = 1 + std::accumulate(_p_impl->_activeSet.begin(),
                                                           _p_impl->_activeSet.end(),
                                                           0,
                                                           [](unsigned lhs, const auto& r_rhs)
                                                            {return r_rhs ? lhs + 1 : lhs;});

        const auto includesAll = [&candidates](const Disc& r_disc) -> bool
        {
            return std::all_of(
                candidates.begin(),
                candidates.end(),
                [&r_disc](const auto& r_optional){
                    return r_optional ? r_disc.isInside(r_optional.value()) : true;
                });
        };

        if (activeSetSize == 1) [[unlikely]]
        {
            disc = Disc(r_point);
            _p_impl->_activeSet[0].emplace(r_point);
        }
        else
        {
            for (const auto& r_permutationSet : Impl::permutationSets[activeSetSize - 2])
            {
                Disc candidate;
                switch (r_permutationSet.size())
                {
                    case 2:
                    {
                        candidate = Disc(candidates[r_permutationSet[0]].value(),
                                         candidates[r_permutationSet[1]].value());
                        break;
                    }
                    case 3:
                    {
                        candidate = Disc(candidates[r_permutationSet[0]].value(),
                                         candidates[r_permutationSet[1]].value(),
                                         candidates[r_permutationSet[2]].value());
                        break;
                    }
                    default:
                    {
                        throw std::runtime_error("Invalid permutation size");
                    }
                }

                if (includesAll(candidate) && candidate.getRadiusSquared() <= disc.getRadiusSquared())
                {
                    disc = candidate;
                    newSet = r_permutationSet;
                }
            }


            for (unsigned index=0; index<newSet.get().size(); ++index)
                _p_impl->_activeSet[index] = candidates[newSet.get()[index]];

            for (unsigned index=newSet.get().size(); index<_p_impl->_activeSet.size(); ++index)
                _p_impl->_activeSet[index].reset();
        }


        _p_impl->_center = disc.getCenter();
        _p_impl->_radiusSquared = disc.getRadiusSquared();
    }

    for (const auto& r_item : _p_impl->_activeSet)
    {
        r_item ? std::cout << "Active set: " << r_item.value()[0] << " " << r_item.value()[1] << " | " : std::cout << " | ";
    }
    std::cout << std::endl;
}


bool MinimumDisc::isIncluded(const Point& r_point) const noexcept
{
    return distanceSquared(r_point, this->getCenter()) <= _p_impl->_radiusSquared;
}


Reference<const MinimumDisc::Point> MinimumDisc::getCenter() const noexcept
{
    return _p_impl->_center;
}


MinimumDisc::Numeric MinimumDisc::getRadius() const
{
    return std::sqrt(_p_impl->_radiusSquared);
}


} // namespace md
