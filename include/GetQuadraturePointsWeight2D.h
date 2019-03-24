#ifndef QUAD2D_H
#define QUAD2D_H

/*
Dunavant points generated with .m code written by John Burkard

http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.html

1. David Dunavant,
    High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle,
    International Journal for Numerical Methods in Engineering,
    Volume 21, 1985, pages 1129-1148.
2. James Lyness, Dennis Jespersen,
    Moderate Degree Symmetric Quadrature Rules for the Triangle,
    Journal of the Institute of Mathematics and its Applications,
    Volume 15, Number 1, February 1975, pages 19-32.

*/

// Note, the coordinates in the x* variables are sequential pairs x,y
//       e.g. x2 = [x21 y21  x22 y22  x23 y23]


// Order 1 Dunavant Points

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/assign/std/vector.hpp>

void GetQuadraturePointsWeight2D(int p, int &n, std::vector<double>& x, std::vector<double>& w);

#endif //QUAD2D_H