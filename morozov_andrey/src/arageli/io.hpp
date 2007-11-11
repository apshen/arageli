/*****************************************************************************

    io.hpp

    This file is a part of Arageli library.

    Copyright (C) 1999--2006 Nikolai Yu. Zolotykh
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


#ifndef _ARAGELI_io_hpp_
#define _ARAGELI_io_hpp_

#include "config.hpp"

#include <iostream>
#include <complex>


namespace Arageli
{


/// Default ouputting method for subexpression as the first coefficient in a polynomial.
template <typename T, typename Ch, typename ChT>
inline std::basic_ostream<Ch, ChT>& output_polynom_first_default
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x
)
{
    return out << x;
}

/// Default ouputting method for subexpression as an internal coefficient in a polynomial.
template <typename T, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& output_polynom_internal_default
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x
);

/// Default ouputting method for subexpression as a degree of variable in a polynomial.
template <typename T, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& output_pow_default
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x
);

/// Default inputting method for subexpression as the first coefficient in a polynomial.
template <typename T, typename Ch, typename ChT>
std::basic_istream<Ch, ChT>& input_polynom_first_default
(
    std::basic_istream<Ch, ChT>& in,
    T& x
);

/// Default inputting method for subexpression as an internal coefficient in a polynomial.
template <typename T, typename Ch, typename ChT>
std::basic_istream<Ch, ChT>& input_polynom_internal_default
(
    std::basic_istream<Ch, ChT>& in,
    T& x
);

/// Default inutting method for subexpression as a degree of variable in a polynomial.
template <typename T, typename Ch, typename ChT>
inline std::basic_istream<Ch, ChT>& input_pow_default
(
    std::basic_istream<Ch, ChT>& in,
    T& x
)
{
    return in >> x;
}


template <typename T, typename Ch, typename ChT>
inline std::basic_ostream<Ch, ChT>& output_polynom_first
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x
)
{
    return output_polynom_first_default(out, x);
}

template <typename T, typename Ch, typename ChT>
inline std::basic_ostream<Ch, ChT>& output_polynom_internal
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x)
{
    return output_polynom_internal_default(out, x);
}

template <typename T, typename Ch, typename ChT>
inline std::basic_ostream<Ch, ChT>& output_pow
(
    std::basic_ostream<Ch, ChT>& out,
    const T& x
)
{
    return output_pow_default(out, x);
}

template <typename T, typename Ch, typename ChT>
inline std::basic_istream<Ch, ChT>& input_polynom_first
(
    std::basic_istream<Ch, ChT>& in,
    T& x
)
{
    return input_polynom_first_default(in, x);
}

template <typename T, typename Ch, typename ChT>
inline std::basic_istream<Ch, ChT>& input_polynom_internal
(
    std::basic_istream<Ch, ChT>& in,
    T& x
)
{
    return input_polynom_internal_default(in, x);
}

template <typename T, typename Ch, typename ChT>
inline std::basic_istream<Ch, ChT>& input_pow
(
    std::basic_istream<Ch, ChT>& in,
    T& x
)
{
    return input_pow_default(in, x);
}


// Specializations for std::complex.


/// Default ouputting method for subexpression as an internal coefficient in a polynomial.
template <typename T, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& output_polynom_internal_default
(
    std::basic_ostream<Ch, ChT>& out,
    const std::complex<T>& x
);

/// Default ouputting method for subexpression as a degree of variable in a polynomial.
template <typename T, typename Ch, typename ChT>
inline std::basic_ostream<Ch, ChT>& output_pow_default
(
    std::basic_ostream<Ch, ChT>& out,
    const std::complex<T>& x
)
{
    out << x;
}


// Binary serialization generic functions.
// They suppose the object can read/write in raw mode,
// i.e. the corresponding type have only trivial constructor/destructor.

/// Binary seft-delimeted serialization for raw objects -- store.
template <typename T, typename Ch, typename ChT>
inline void output_binary (std::basic_ostream<Ch, ChT>& out, const T& x)
{
    out.write(reinterpret_cast<const char*>(&x), sizeof(T));
}


/// Binary seft-delimeted serialization for raw objects -- load.
template <typename T, typename Ch, typename ChT>
inline void input_binary (std::basic_istream<Ch, ChT>& in, T& x)
{
    in.read(reinterpret_cast<char*>(&x), sizeof(T));
}


/// Binary seft-delimeted serialization for arrays of raw objects -- store.
template <typename T, typename Ch, typename ChT>
inline void output_binary
(
    std::basic_ostream<Ch, ChT>& out,
    const T* x,
    std::size_t n
)
{
    if(n)
        out.write(reinterpret_cast<const char*>(x), n*sizeof(T));
}


/// Binary seft-delimeted serialization for arrays of raw objects -- load.
template <typename T, typename Ch, typename ChT>
inline void input_binary
(
    std::basic_istream<Ch, ChT>& in,
    T* x,
    std::size_t n
)
{
    if(n)
        in.read(reinterpret_cast<char*>(x), n*sizeof(T));
}


} // namespace Arageli


#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_IO
    #include "io.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_IO
#endif


#endif  //  #ifndef _ARAGELI_io_hpp_
