#ifndef MB_BASE_PLUGIN
#define MB_BASE_PLUGIN <adtest_plugins.h>
#endif

/* plugins must be defined before including MB_Base, which is included
   in the model file
*/

#include <adtest.h>


RCPP_MODULE(adtest){  
  Rcpp::class_< MB_Base<adtest> >("adtest")
    
    .constructor<const List>()
    
    .method( "get.f", & MB_Base<adtest>::get_f)
    .method( "get.df", & MB_Base<adtest>::get_df)
    .method( "get.fdf", & MB_Base<adtest>::get_fdf)
    .method( "get.fdfh", & MB_Base<adtest>::get_fdfh)
    .method( "get.tape.stats", & MB_Base<adtest>::get_tape_stats)
    .method( "get.hessian", & MB_Base<adtest>::get_hessian)
    .method( "get.hessian.sparse", & MB_Base<adtest>::get_hessian_sparse)
    .method( "record.tape", & MB_Base<adtest>::record_tape)
    .method( "init.hessian", & MB_Base<adtest>::hessian_init)
    .method( "init.sparse.hessian", & MB_Base<adtest>::hessian_init)
    .method( "get.f.direct", & MB_Base<adtest>::get_f_direct) // plugin
    .method( "get.LL", & MB_Base<adtest>::get_LL) // plugin
    .method( "get.LLi", & MB_Base<adtest>::get_LLi) // plugin
    .method( "get.prior.i", & MB_Base<adtest>::get_prior_i) // plugin
    ;
}
