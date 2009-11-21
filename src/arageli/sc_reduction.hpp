/**
    \file sc_reduction.hpp
    \brief Shevchenko-Chirkov basis reduction.
*/

#ifndef _ARAGELI_sc_reduction_hpp_
#define _ARAGELI_sc_reduction_hpp_

#include "config.hpp"


#include "std_import.hpp"


namespace Arageli
{

/// Shevchenko-Chirkov basis reduction.
/**
    Function implements algorithm for finding c-reduced basis.

    Input matrix B must contain vectors of an initial basis (in columns).
    On output B will contain the vectors of reduced basis.
    H will be the corresponding transformation matrix such that B_new = B_old * H.
    c is a constant of reduction.

    If columns of B are not linear independent then the function return false.
    Otherwise it returns true.

    The entries of B and H must be belong to some subfieled of the field of real numbers.
    The entry c must be a rational constant greater than 4/3.

        @param B should be a matrix
        @param H should be a matrix
        @param c should be a rational number

*/
template <typename B_type, typename H_type, typename c_type>
bool sc_reduction(B_type &B, H_type &H, c_type c);

}

#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SC_REDUCTION
    #include "sc_reduction.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SC_REDUCTION
#endif


#endif  // #ifndef _ARAGELI_sc_reduction_hpp_

