#include "../include/InvertMatrix.h"

template <class T>
bool InvertMatrix(ublas::matrix<T> &input, ublas::matrix<T> &inverse)
{
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input
	matrix<T> A(input);
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;
	// create identity matrix of "inverse"
	inverse.assign(ublas::identity_matrix<T>(A.size1()));
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	return true;
}
template bool InvertMatrix(ublas::matrix<double> &input, ublas::matrix<double> &inverse);