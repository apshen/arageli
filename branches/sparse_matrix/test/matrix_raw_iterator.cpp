/*****************************************************************************

    test/matrix_raw_iterator.cpp

    This file is a part of the Arageli library test base.

    Copyright (C) 2005--2007 Sergey S. Lyalin
    University of Nizhni Novgorod, Russia

    The Arageli Library is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License version 2
    as published by the Free Software Foundation.

    The Arageli Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

    We are also open for dual licensing for the whole library or
    for its particular part. If you are interested to get the library
    in this way, i.e. not under the GNU General Public License,
    please contact Arageli Support Service support.arageli@gmail.com.

*****************************************************************************/

/**
    \file matrix_raw_iterator.cpp
    \brief This file includes test for matrix_raw_iterator.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "stdafx.hpp"

using namespace Arageli;


namespace
{

template <typename M, typename I, typename Iind>
struct Row_raw_adaptor
{
    typedef typename M::size_type size_type;
    typedef I iterator;
    typedef Iind ind_iterator;
    typedef typename M::value_type value_type;

    size_type n (M& a) const
    {
        return a.nrows();
    }

    size_type m (M& a) const
    {
        return a.ncols();
    }

    iterator begin (M& a, size_type i) const
    {
        return a.row_raw_begin(i);
    }

    ind_iterator ind_end (M& a, size_type i) const
    {
        return a.row_indraw_end(i);
    }

    ind_iterator ind_begin (M& a, size_type i) const
    {
        return a.row_indraw_begin(i);
    }

    iterator end (M& a, size_type i) const
    {
        return a.row_raw_end(i);
    }

    const value_type& el (M& a, size_type i, size_type j) const
    {
        return a(i, j);
    }

    size_type index (ind_iterator j) const
    {
        return j.icol();
    }
};


template <typename M, typename I, typename Iind>
struct Col_raw_adaptor
{
    typedef typename M::size_type size_type;
    typedef I iterator;
    typedef Iind ind_iterator;
    typedef typename M::value_type value_type;

    size_type n (M& a) const
    {
        return a.ncols();
    }

    size_type m (M& a) const
    {
        return a.nrows();
    }

    iterator begin (M& a, size_type i) const
    {
        return a.col_raw_begin(i);
    }

    ind_iterator ind_end (M& a, size_type i) const
    {
        return a.col_indraw_end(i);
    }

    ind_iterator ind_begin (M& a, size_type i) const
    {
        return a.col_indraw_begin(i);
    }

    iterator end (M& a, size_type i) const
    {
        return a.col_raw_end(i);
    }

    const value_type& el (M& a, size_type i, size_type j) const
    {
        return a(j, i);
    }

    size_type index (ind_iterator j) const
    {
        return j.irow();
    }
};


template <typename M, typename Adaptor>
bool iterate_for_all_cols_or_rows_helper_1 (M& a, Adaptor adaptor)
{
    typedef typename Adaptor::size_type size_type;

    for(size_type i = 0; i < adaptor.n(a); ++i)
    {
        // Check sizes and applicability of STL algorithms for element counting.

        size_type raw_size = std::distance(adaptor.begin(a, i), adaptor.end(a, i));
        if(raw_size > adaptor.m(a))
        {
            tout
                << "Filed at " << __LINE__ << '\n'
                << "i = " << i << '\n'
                << "raw_size = " << raw_size << '\n'
                << "size of the matrix = " << adaptor.m(a) << '\n'
                << "typeid(M).name() = " << typeid(M).name() << '\n';
            output_aligned(tout, a);
            tout << '\n';

            return false;
        }

        std::map<typename M::value_type, size_type> map1, map2;

        for(size_type j = 0; j < adaptor.m(a); ++j)
            if(!is_null(adaptor.el(a, i, j)))
                ++(map1[adaptor.el(a, i, j)]);

        for(typename Adaptor::iterator j = adaptor.begin(a, i); j != adaptor.end(a, i); ++j)
            if(!is_null(*j))
                ++(map2[*j]);

        if(map1 != map2)
        {
            tout
                << "Failed at " << __LINE__ << '\n'
                << "i = " << i << '\n'
                << "typeid(M).name() = " << typeid(M).name() << '\n';
            output_aligned(tout, a);
            tout << '\n';

            return false;
        }
    }

    return true;
}


template <typename M, typename Adaptor>
bool check_row_or_col_indices (M& a, Adaptor adaptor)
{
    typedef typename Adaptor::size_type size_type;

    for(size_type i = 0; i < adaptor.n(a); ++i)
    {
        // Iterate for all elements of the row/col by an index raw iterator.

        for(typename Adaptor::ind_iterator j = adaptor.ind_begin(a, i); j != adaptor.ind_end(a, i); ++j)
        {
            if(adaptor.el(a, i, adaptor.index(j)) != j.el())
            {
                tout
                    << "Filed at " << __LINE__ << '\n'
                    << "i = " << i << '\n'
                    << "adaptor.index(j) = " << adaptor.index(j) << '\n'
                    << "adaptor.el(a, i, adaptor.index(j)) = " << adaptor.el(a, i, adaptor.index(j)) << '\n'
                    << "j->el() = " << j.el() << '\n'
                    << "typeid(M).name() = " << typeid(M).name() << '\n';
                output_aligned(tout, a);
                tout << '\n';
                
                return false;
            }
        }
    }

    return true;
}

template <typename M>
bool iterate_for_all_cols_and_rows (M& a)
{
    return
        iterate_for_all_cols_or_rows_helper_1
        (
            a,
            Row_raw_adaptor
            <
                M,
                typename M::row_raw_iterator,
                typename M::row_indraw_iterator
            >()
        ) &&
        iterate_for_all_cols_or_rows_helper_1
        (
            (const M)a,
            Row_raw_adaptor
            <
                const M,
                typename M::const_row_raw_iterator,
                typename M::const_row_indraw_iterator
            >()
        ) &&
        check_row_or_col_indices
        (
            a,
            Row_raw_adaptor
            <
                M,
                typename M::row_raw_iterator,
                typename M::row_indraw_iterator
            >()
        ) &&
        check_row_or_col_indices
        (
            (const M)a,
            Row_raw_adaptor
            <
                const M,
                typename M::const_row_raw_iterator,
                typename M::const_row_indraw_iterator
            >()
        ) &&
        iterate_for_all_cols_or_rows_helper_1
        (
            a,
            Col_raw_adaptor
            <
                M,
                typename M::col_raw_iterator,
                typename M::col_indraw_iterator
            >()
        ) &&
        iterate_for_all_cols_or_rows_helper_1
        (
            (const M)a,
            Col_raw_adaptor
            <
                const M,
                typename M::const_col_raw_iterator,
                typename M::const_col_indraw_iterator
            >()
        ) &&
        check_row_or_col_indices
        (
            a,
            Col_raw_adaptor
            <
                M,
                typename M::col_raw_iterator,
                typename M::col_indraw_iterator
            >()
        ) &&
        check_row_or_col_indices
        (
            (const M)a,
            Col_raw_adaptor
            <
                const M,
                typename M::const_col_raw_iterator,
                typename M::const_col_indraw_iterator
            >()
        );
}

template <typename M>
bool concrete_test_1 ()
{
    {
        M a = "((1, 0, 2), (0, 0, 0), (0, 0, 1), (3, 2, 1))";
        if(!iterate_for_all_cols_and_rows(a))
            return false;
    }
    {
        M a = "((), (), (), ())";
        if(!iterate_for_all_cols_and_rows(a))
            return false;
    }
    {
        M a = "()";
        if(!iterate_for_all_cols_and_rows(a))
            return false;
    }
    {
        M a = "((1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12))";
        if(!iterate_for_all_cols_and_rows(a))
            return false;
    }
    {
        M a = "((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12))";
        if(!iterate_for_all_cols_and_rows(a))
            return false;
    }

    return true;
}

}


TEST
(
    matrix,
    work_with_raw_iterators,
    "Test for raw iterators for matrix class."
)
{
    bool is_ok = true;

    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_BEGIN
    {
        is_ok &= concrete_test_1<matrix<int> >();
        is_ok &= concrete_test_1<matrix<big_int> >();
        is_ok &= concrete_test_1<matrix<rational<int> > >();
        is_ok &= concrete_test_1<matrix<rational<big_int> > >();
    }
    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END

    return is_ok ? resOK : resFAIL;
}
