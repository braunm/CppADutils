typedef Eigen::Triplet<double> TT;
typedef std::vector<TT> VectorTT;

template <class VectorBase, class VectorSet>
VectorTT SparseHessian(
		       const VectorBase& x,
		       const VectorBase& w,
		       const VectorSet& p,
		       CppAD::sparse_hessian_work& work) {
  
  size_t i, j, k;
  
  size_t n = x.size();
  
  VectorTT hes;
  
  typedef typename VectorSet::value_type Set_type;
  typedef typename CppAD::internal_sparsity<Set_type>::pattern_type Pattern_type;
  
  
  // arguments to SparseHessianCompute
  Pattern_type          s;
  CppAD::vector<size_t> row;
  CppAD::vector<size_t> col; 
  bool transpose = false;
  sparsity_user2internal(s, p, n, n, transpose);
  k = 0;
  for(i = 0; i < n; i++)
    {	s.begin(i);
      j = s.next_element();
      while( j != s.end() )
	{	row.push_back(i);
	  col.push_back(j);
	  k++;
	  j = s.next_element();
	}
    }
  size_t K = k;
  VectorBase H(K);
  
  // now we have folded this into the following case
  SparseHessianCompute(x, w, s, row, col, H, work);
  
  hes.resize(K); // a Triplet for each non-zero
  
  // now set the non-zero return values
  for(k = 0; k < K; k++) 
    hes[k] = TT(row[k], col[k], H[k]);
  
  return hes;
}
