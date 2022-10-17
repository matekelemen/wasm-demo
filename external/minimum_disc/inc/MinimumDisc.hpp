#pragma once

// --- Internal Includes ---
#include "typedefs.hpp"

// --- STL Includes ---
#include <array>
#include <vector>
#include <memory>


namespace md {


class MinimumDisc
{
public:
    using Numeric = double;

    using Point = std::array<Numeric,2>;

    using PointContainer = std::vector<Point>;

public:
    MinimumDisc();

    ~MinimumDisc();

    void include(Reference<const Point> r_point);

    Point getCenter() const;

    Numeric getRadius() const;

private:
    struct Impl;
    std::unique_ptr<Impl> _p_impl; // <== pointer to implementation
}; // class MinimumDisc


} // namespace md
