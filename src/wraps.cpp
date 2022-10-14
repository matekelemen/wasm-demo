// --- External Includes ---
#include "MinimumDisc.hpp"
#include "emscripten/bind.h"

// --- Internal Includes ---
#include "common.hpp"

// --- STL Includes ---
#include <iterator>
#include <memory>
#include <vector>
#include <algorithm>


// Emscripten doesn't have bindings for std::array, so we need to
// avoid exposing them to JS. One way is to write decorators for
// the members that return or take std::arrays.
namespace MinimumDisc {


void include(md::MinimumDisc& r_this, double x, double y)
{
    r_this.include({x, y});
}


std::vector<double> getCenter(md::MinimumDisc& r_this)
{
    const auto& r_center = r_this.getCenter();
    return {r_center[0], r_center[1]};
}


} // namespace MinimumDisc


EMSCRIPTEN_BINDINGS(WASMModule)
{
    emscripten::register_vector<double>("VectorDouble");

    emscripten::class_<md::MinimumDisc>("MinimumDisc")
        .constructor<>()
        .function("include", &MinimumDisc::include)
        .function("getCenter", &MinimumDisc::getCenter)
        .function("getRadius", &md::MinimumDisc::getRadius)
        ;
}
