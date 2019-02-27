/*
The following code inverts the matrix input using LU-decomposition with backsubstitution of unit vectors. Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling & Flannery.
you can solve Ax=b using three lines of ublas code:
permutation_matrix<> piv;
lu_factorize(A, piv);
lu_substitute(A, piv, x);
*/
#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

/* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template <class T>
bool InvertMatrix(ublas::matrix<T> &input, ublas::matrix<T> &inverse);

#endif //INVERT_MATRIX_HPP