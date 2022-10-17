// --- Internal Includes ---
#include "MinimumDisc.hpp"

// --- STL Includes ---
#include <iostream>


int main()
{
    std::cout << "Begin" << std::endl;
    md::MinimumDisc disc;
    disc.include({0, 0});
    disc.include({2, 2});

    std::cout << "Center: ";
    for (auto component : disc.getCenter())
        std::cout << component << " ";
    std::cout << std::endl
    << "Radius: " << disc.getRadius() << std::endl;
}
