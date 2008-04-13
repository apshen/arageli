/*****************************************************************************

    sparse_matrix.hpp

    This file is a part of the Arageli library.

    Copyright (C) 2005--2008 Sergey S. Lyalin
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
    \file sparse_matrix.hpp
    \brief Sparse matrix and supporting tools declaration.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#ifndef _ARAGELI_sparse_matrix_hpp_
#define _ARAGELI_sparse_matrix_hpp_

#include "config.hpp"
#include "misc.hpp"
#include "sparse_vector.hpp"

// REFERENCE ADDITIONAL HEADERS HERE


namespace Arageli
{


namespace cnfg
{

template <typename T, typename RepType, bool REFCNT = true>
struct sparse_matrix_default
{
    static const bool refcounting = REFCNT;
    typedef RepType line_representation_tag;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, RepType, false> line_configurator;
};

}


template <typename T, typename RepType, bool REFCNT>
struct type_traits<cnfg::sparse_matrix_default<T, RepType, REFCNT> > :
    public type_traits_default<cnfg::sparse_matrix_default<T, RepType, REFCNT> >
{
    static const bool is_specialized = true;
    typedef type_category::configurator category_type;
    static const category_type category_value;
};


namespace spct
{

/// Default spectator for sparse_matrix class.
struct sparse_matrix_idler
{};

}


template <>
struct type_traits<spct::sparse_matrix_idler> :
    public type_traits_default<spct::sparse_matrix_idler>
{
    static const bool is_specialized = true;
    typedef type_category::spectator category_type;
    static const category_type category_value;
};


namespace _Internal
{


template <typename T, typename Param, typename Category>
sparse_matrix_param_extractor_helper_1;


/// Gives Param as configurator and set the default spectator.
template <typename T, typename Param>
sparse_matrix_param_extractor_helper_1<T, Param, type_category::configurator>
{
    /// Choosen configurator.
    typedef Param configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};


/// Gives Param as spectator and set the default configurator.
template <typename T, typename Param>
sparse_matrix_param_extractor_helper_1<T, Param, type_category::spectator>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T> configurator;

    /// Choosen spectator.
    typedef Param spectator;
};


/// Uses Param as a tag to choose the representation type.
template <typename T, typename Param>
sparse_matrix_param_extractor_helper_1<T, Param, type_category::tag>
{
    /// Default configurator with choosen type of representation.
    typedef cnfg::sparse_matrix_default<T, Param> configurator;

    /// Choosen spectator.
    typedef spct::sparse_matrix_idler spectator;
};


template <typename T, typename Param>
struct sparse_matrix_param_extractor
{
    /// Determined configurator.
    typedef typename sparse_matrix_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::configurator
        configurator;

    /// Determined spectator.
    typedef typename sparse_matrix_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::spectator
        spectator;
};


template <typename T>
struct sparse_matrix_param_extractor<T, default_tag>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};


}


/// Sparse matrix with elements of type T.
/** Way of implementation and other configuration parameters are defined
    by the second template parameter (Param). */
template
<
    typename T,
    typename Param = default_tag
>
class sparse_matrix
{


public:
};


} // namesapce Arageli


#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_sparse_matrix
    #include "sparse_matrix.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_sparse_matrix
#endif

#endif    // #ifndef _ARAGELI_sparse_matrix_hpp_
