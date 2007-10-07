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

/// Stores a raw object state to a binary stream. Seft-delimeted binary serialization.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename Stream, typename T>
inline Stream& output_binary_stream (Stream& out, const T& x)
{
    out.write(reinterpret_cast<const char*>(&x), sizeof(T));
    return out;
}


/// Loads a raw object state from a binary stream. Compatible with output_binary_stream.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename Stream, typename T>
inline Stream& input_binary_stream (Stream& in, T& x)
{
    in.read(reinterpret_cast<char*>(&x), sizeof(T));
    return in;
}


/// Stores an array of a raw object states to a binary stream. Seft-delimeted binary serialization.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename Stream, typename T>
inline Stream& output_binary_stream
(
    Stream& out,
    const T* x,
    std::size_t n
)
{
    if(n)
        out.write(reinterpret_cast<const char*>(x), n*sizeof(T));
    return out;
}


/// Loads an array of a raw object states from a binary stream. Compatible with output_binary_stream.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename Stream, typename T>
inline Stream& input_binary_stream
(
    Stream& in,
    T* x,
    std::size_t n
)
{
    if(n)
        in.read(reinterpret_cast<char*>(x), n*sizeof(T));
    return in;
}


/// Calculates the number of chars required to store a given object in The Simple Binary form.
/** This function calculates precise number of chars that will emit
    any function outputs in The Simple Binary format for one object,
    for example, output_binary_mem function. */
template <typename T>
inline std::size_t calc_binary (const T& x)
{
    return sizeof(T);
}

/// Calculates the number of chars required to store a given array of objects in The Simple Binary form.
/** This function calculates precise number of chars that will emit
    any function outputs in The Simple Binary format for an array of
    objects, for example, output_binary_mem function. */
template <typename T>
inline std::size_t calc_binary (const T* x, std::size_t n)
{
    return n*sizeof(T);
}


/// Stores a raw object state to a memory location. Seft-delimeted binary serialization.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename T>
inline char* output_binary_mem (char* out, const T& x)
{
    return
        std::copy
        (
            reinterpret_cast<const char*>(&x),
            reinterpret_cast<const char*>(&x) + sizeof(T),
            out
        );
}


/// Loads a raw object state from a memory location. Compatible with output_binary_stream.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename T>
inline const char* input_binary_mem (const char* in, T& x)
{
    const char* in_end = in + sizeof(T);

    std::copy
    (
        in,
        in_end,
        reinterpret_cast<char*>(&x)
    );

    return in_end;
}


/// Stores an array of a raw object states to a memory location. Seft-delimeted binary serialization.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename T>
inline char* output_binary_mem
(
    char* out,
    const T* x,
    std::size_t n
)
{
    return
        std::copy
        (
            reinterpret_cast<const char*>(x),
            reinterpret_cast<const char*>(x) + n*sizeof(T),
            out
        );
}


/// Loads an array of a raw object states from a memory location. Compatible with output_binary_stream.
/** This function supposes T doesn't have a non trivial constructor/destructor pair.
    The format that this function uses is call The Simple Binary format. */
template <typename T>
inline const char* input_binary_mem
(
    const char* in,
    T* x,
    std::size_t n
)
{
    const char* in_end = in + n*sizeof(T);

    std::copy
    (
        in,
        in_end,
        reinterpret_cast<char*>(x)
    );

    return in_end;
}


} // namespace Arageli


#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_IO
    #include "io.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_IO
#endif


#endif  //  #ifndef _ARAGELI_io_hpp_
