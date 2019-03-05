#include "../include/geometry.h"

namespace geometry
{
    double BumpFunction(double x)
    {
        return 0.0625 * std::exp(-25.0 * x * x);
    }



} // namespace geometry
