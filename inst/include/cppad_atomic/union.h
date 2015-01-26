

#ifndef __MY_ATOMIC_UNION
#define __MY_ATOMIC_UNION

inline void my_union( 
	      std::set<size_t>&         result  , 
	      const std::set<size_t>&   left    , 
	      const std::set<size_t>&   right   ) 
{  
  std::set<size_t> temp; 
  std::set_union( 
		 left.begin()              , 
		 left.end()                , 
		 right.begin()             , 
		 right.end()               , 
		 std::inserter(temp, temp.begin()) 
		  ); 
  result.swap(temp); 
}

#endif
