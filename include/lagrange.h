#ifndef LAGRANGE
#define LAGRANGE

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace lagrange
{

    namespace ublas = boost::numeric::ublas;

    ublas::matrix<double> TriangleLagrange2D(int p);

    ublas::vector<ublas::vector<double> > MapReferenceToPhysicalLinear(ublas::vector<ublas::vector<double> > vertex, int p);

    ublas::vector<double> MapPhysicalToReferenceLinear(ublas::vector<ublas::vector<double> > vertex, ublas::vector<double> point, int p);

} // namespace lagrange

#endif