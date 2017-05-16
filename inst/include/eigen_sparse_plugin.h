
# ifndef CPPAD_EIGEN_SPARSE_PLUGIN_INCLUDED
# define CPPAD_EIGEN_SPARSE_PLUGIN_INCLUDED

// typedef Scalar value_type;

# endif

template<typename S, int R, int Opt, int MR, int MC>
inline void reserve(const Matrix<S,R,Opt,MR,MC> &reserveSizes)
{
  reserveInnerVectors(reserveSizes);
}


