/*****************************************************************************

    big_int.hpp

    This file is a part of the Arageli library.

    WARNIG. This file has no complete implementation.

    Copyright (C) 1999 -- 2005 Nikolai Yu. Zolotykh
    Copyright (C) 2006 -- 2010 Sergey S. Lyalin

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

#ifndef __ARAGELI_big_int_h__
#define __ARAGELI_big_int_h__

#include "config.hpp"

#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cctype>
#include <cmath>
#include <limits>

#include "frwrddecl.hpp"

#include "type_traits.hpp"
#include "type_pair_traits.hpp"
#include "exception.hpp"
#include "factory.hpp"
#include "io.hpp"
#include "powerest.hpp"
#include "bigar.hpp"
#include "gcd.hpp"

#include "std_import.hpp"

/**
    \file big_int.hpp.
    \brief Big integer number class implementation.

    This module implements a class big_int for representing
    Big Integer Numbers, i.e. `arbitrary precision' integers.
    Implements some common library and standard C++ functions
    specialization and reloaded versions.
*/

namespace Arageli
{

namespace _Internal
{
    // a = b / c; res = b % c
    void xdivide (big_int& a, const big_int& b, const big_int& c, big_int& res);
}


// The following declaration is needed only for one additional definition
// in big_int class.  That definition is created for GCC compatibility only.
// See template <typename TT2> big_int::big_int (const rational<big_int>& x).

// template <typename T> class rational;    // see rational.hpp


template <>
class io_binary<big_int>;   // see implementation below in this file


template <typename T_factory>
big_int gcd (const big_int& a, const big_int& b, const T_factory& tfctr);


template <typename T_factory>
big_int gcd (const big_int& a, const big_int& b, const T_factory& tfctr);


/// Big integer number.
/** Integer number with arbitrary big size. Only available memory
    limits the maximum value of a number. */
class big_int
{
    // A type for one elementary digit.
    // Number consist of sequence of objects of this type.
    typedef _Internal::digit digit;

    friend class io_binary<big_int>;

    template <typename T_factory>
    friend big_int gcd (const big_int& a, const big_int& b, const T_factory& tfctr);

    friend std::size_t magnitude (const big_int& x);

public:

    /// Exception of class big_int.
    /** All other exception class inherit from this class. */
    class exception : public virtual Arageli::exception
    {
        std::string msg_m;
    public:

        /// Initialization with empty message (empty string).
        exception () {}

        /// Initialization with a particular message string.
        /** @param msg_a a message. */
        exception (const std::string& msg_a) :
            msg_m(msg_a)
        {}

        /// Returns the message of exception.
        /** The result is equal to initialization */
        virtual std::string msg () const
        {
            return msg_m;
        }
    };

    // WARNING!  A message storage is duplicated in the incorrec_string class.
    // Realy two verstions of this storage are presented in the class.
    // TODO: Resolve it.

    /// Exeption. The throwing situation is incorrect format of string.
    /** When exception will throwed: an input string
        for converting to big_int has invalid format.*/
    class incorrect_string : public exception, public Arageli::incorrect_string
    {
    public:
        incorrect_string (const std::string& msg_a) :
            Arageli::incorrect_string(msg_a) {}

        virtual std::string msg () const
        {
            return Arageli::incorrect_string::msg();
        }
    };


    /// Constructs a zero.
    /** After this initialization call of is_null returns true. */
    big_int () = default;

    /// Converts a string to a big_int.
    /** Conversions is based on operator>> for std::istream.
        See rules for the string there.
        @param str a string for initialization */
    big_int (const char* str);

    /// Copy constructor
    big_int (const big_int& b)
    {
        ARAGELI_ASSERT_0(this != &b);

        if(b.number.sign != 0)
        {
            number = big_struct(b.number.sign, b.number.len, b.number.len);
            copy_data(number.digits(), b.number.digits(), b.number.len);
        }
    }

    /// Move constructor
    big_int (big_int&& b) = default;

    /// Additional constructor.  For GCC compatibility only!!!
    /** I do not really know why this is needed for GCC, but... it is needed!!!
        If you know this, please, email me on
        <a href="mailto:sergey.lyalin@gmail.com">sergey.lyalin@gmail.com</a>.
        Without this constructor compilation of construct from rational
        will be failed. */
    template <typename T>
    big_int (const rational<T>& x)
    {
        *this = x.numerator()/x.denominator();
    }

    /**    @name Initializations from built-in types. */
    //@{

    template<typename I,
             typename std::enable_if<std::is_integral<I>::value, bool>::type = true>
    big_int (I x)
    {
        from_native_int(x);
    }

    template<typename F,
             typename std::enable_if<std::is_floating_point<F>::value, bool>::type = true>
    big_int (F x)
    {
        from_native_float(x);
    }

    //@}


    /**    @name Convertions to built-in types. */
    //@{

    template<typename I,
             typename std::enable_if<std::is_integral<I>::value, bool>::type = true>
    operator I () const
    {
        return to_native_int<I>();
    }

    operator bool () const
    {
        return !is_null();
    }

    template<typename F,
             typename std::enable_if<std::is_floating_point<F>::value, bool>::type = true>
    operator F () const
    {
        return to_native_float<F>();
    }

    //@}


    /// Returns true if the number is not zero and false otherwise.
    bool operator! () const
    {
        return is_null();
    }

    ~big_int () = default;

    /// Copy assignment
    /** This is assignment operator
     *
     *   @param b a new value
    */
    big_int& operator= (const big_int & b);

    /// Move assignment
    /** This is move-assignment operator
     *
     *   @param b a new value
    */
    big_int& operator= (big_int && b) = default;

    big_int& operator= (const char* s)
    {
        return *this = big_int(s);
    }

    template<typename I,
             typename std::enable_if<std::is_integral<I>::value, bool>::type = true>
    big_int& operator= (I x)
    {
        from_native_int(x);
        return *this;
    }

    template<typename F,
             typename std::enable_if<std::is_floating_point<F>::value, bool>::type = true>
    big_int& operator= (F x)
    {
        from_native_float(x);
        return *this;
    }

    template <typename T>
    big_int& operator= (const rational<T>& x)
    {
        return *this = big_int(x);
    }

    /// Returns the number of bits in number.
    /** Complexity: O(log(n)), where n is the current value of number. */
    std::size_t length () const;

    /// Access to specific bit in this number.
    /**    Requirements: 0 <= k < length().
        @param k a position of bit, where k = 0 points to lower bit.*/
    bool operator[] (std::size_t k) const;

    ///    Sign of a number.
    /**    Returns
        -   0    if x = 0,
        -  -1    if x < 0,
        -  +1    if x > 0. */
    int sign () const
    {
        return number.sign;
    }

    /// Return true if this number is zero.
    inline bool is_null () const;

    /// Return true if this number is unit.
    inline bool is_unit () const;

    /// Return true if this number is opposite unit.
    inline bool is_opposite_unit () const;

    /// Unary plus.
    /** Returns a reference on this object.  It does not anything. */
    const big_int& operator+ () const
    {
        return *this;
    }

    /// Unary minus.
    big_int operator- () const;

    /// Bitwise inversion. NOT IMPLIMENTED.
    big_int operator~ () const
    {
        ARAGELI_ASSERT_ALWAYS(!"big_int::operator~ isn't implemented yet");
        return big_int();
    }

    /// Prefix form of an increment.
    /** Increments this number by 1 and returns a reference on this object. */
    big_int& operator++ ()
    {
        return *this = *this + big_int(1);
    }

    /// Prefix form of a decrement.
    /** Decrements this number by 1 and returns a reference on this object. */
    big_int& operator-- ()
    {
        return *this = *this - big_int(1);
    }

    /// Postfix form of an increment.
    /** Increments this number by 1 and returns an old value of this object. */
    big_int operator++ (int)
    {
        big_int t = *this;
        operator++();
        return t;
    }

    /// Postfix form of a decrement.
    /** Decrements this number by 1 and returns an old value of this object. */
    big_int operator-- (int)
    {
        big_int t = *this;
        operator--();
        return t;
    }

    friend big_int operator+ (const big_int & b, const big_int & c);
    friend big_int operator- (const big_int & b, const big_int & c);
    friend big_int operator* (const big_int & b, const big_int & c);
    friend big_int operator/ (const big_int & b, const big_int & c);
    friend big_int operator% (const big_int & b, const big_int & c);
    friend big_int operator& (const big_int & b, const big_int & c);
    friend big_int operator| (const big_int & b, const big_int & c);
    friend big_int operator^ (const big_int & b, const big_int & c);
    friend std::ostream& operator<< (std::ostream& s, const big_int& x);
    friend std::istream& operator>> (std::istream& s, big_int& x);
    friend int cmp (const big_int & a, const big_int & b);

    friend void _Internal::xdivide
    (
        big_int& a,
        const big_int& b,
        const big_int& c,
        big_int& res
    );

    /// Returns pseudo-random number with exact len bits.
    /** Generating is based on std::rand function.
        Random number in range [2^(len - 1), 2^len - 1].
        Each number from range is equiprobable. */
    static big_int random_with_length (std::size_t len)
    {
        return
            len ?
            random_with_length_or_less(len - 1) + (big_int(1) << (len - 1)) :
            big_int();
    }

    /// Returns pseudo-random number in range [0, max].
    /** Generating is based on std::rand function.
        Each number from range [0, max] is equiprobable. */
    static big_int random_in_range (const big_int& max);

    /// Returns pseudo-random number with len or less then len bits.
    /** Generating is based on std::rand function.
        Random number in range [0, 2^len - 1].
        Each number from range is equiprobable. */
    static big_int random_with_length_or_less (std::size_t len);

    friend big_int operator<< (const big_int& a, std::size_t n);
    friend big_int operator>> (const big_int& a, std::size_t n);
    friend big_int& operator>>= (big_int& a, std::size_t n);


    /// Returns the k-th bit of the number
    bool bit (std::size_t k) const
    {
        return operator[](k);
    }

    /// Returns true if even, false if odd
    bool is_even () const
    {
        if(number.len == 0)
            return true;
        else
            return !(operator[](0));
    }


    /// Returns true if odd, false if even
    bool is_odd () const
    {
        if(!number.len)
            return false;
        else
            return operator[](0);
    }

    /// Swaps two numbers.
    void swap (big_int& x)
    {
        std::swap(number, x.number);
    }

    const digit* _digits () const
    {
        return number.digits();
    }

private:

    friend class big_float;    // see big_float.hpp source file

    // the following type is used inside the big_int unit and implements
    // the storage for a Big Integer Number

    struct digit_storage
    {
        digit *mbeg;        // the storage for digits
        digit *mend;        // the end of storage
        std::size_t len;    // the number of used digits

        digit_storage() : mbeg(nullptr), mend(nullptr), len(0)
        {
        }

        digit_storage(std::size_t s, std::size_t l) :
            mbeg(allocate(s)),
            mend(mbeg + s),
            len(l)
        {
            ARAGELI_ASSERT_1(s >= l);
            ARAGELI_ASSERT_1(l > 0);
        }

        ~digit_storage()
        {
            deallocate(mbeg);
        }

        digit_storage(digit_storage &&x) : digit_storage()
        {
            swap(x);
        }

        digit_storage &operator = (digit_storage &&x)
        {
            swap(x);

            return *this;
        }

        digit &operator[](std::size_t i) const
        {
            return mbeg[i];
        }

        void swap(digit_storage &x)
        {
            std::swap(mbeg, x.mbeg);
            std::swap(mend, x.mend);
            std::swap(len, x.len);
        }

    private:
        digit *allocate(size_t s)
        {
            digit *p = static_cast<digit *>(std::malloc(s * sizeof(digit)));
            if(!p)
                big_int::exception("Big arith error: the heap overflow");
            return p;
        }

        void deallocate(digit *m)
        {
            if(m)
                std::free(m);
        }

    };

    struct big_struct : public digit_storage
    {
        int sign;           // the sign: 0, 1 or -1

        big_struct() : sign(0)
        {
        }

        big_struct(big_struct &&) = default;
        big_struct &operator = (big_struct &&) = default;

        big_struct(int s, std::size_t c, std::size_t l) :
            digit_storage(c, l),
            sign(s)
        {
        }

        digit *digits() const
        {
            return mbeg;
        }

        void set_zero()
        {
            *this = big_struct();
        }
    } number;

    static void copy_data
    (
        digit* dest,
        const digit* source,
        std::size_t newnitems
    );

    big_int(int s, std::size_t c, std::size_t l) : number(s, c, l)
    {
    }

    // Erases leading zeros.
    static void optimize (big_struct &numner);

    template <typename T>
    void from_native_int (const T& x);

    template <typename T>
    void from_native_int_helper (const T& x, true_type);

    template <typename T>
    void from_native_int_helper (const T& x, false_type);

    template <typename T>
    void from_native_float (const T& x);

    template <typename T>
    T to_native_int () const;

    template <typename T>
    T to_native_float () const;

    template <typename T>
    T to_native_int_without_sign () const;

};


/// Returns the number of bits in number. The same as big_int::length.
inline std::size_t nbits (const big_int& x)
{
    return x.length();
}


/// Reads a number from a string notation.
std::ostream & operator<< (std::ostream & s, const big_int & x);

/// Writes a number to a string notation.
std::istream & operator>> (std::istream & s, big_int & x);


/// Serialization in the Simple Binary format for big_int.
template <>
class io_binary<big_int> : public io_binary_base<big_int>
{
public:

    using io_binary_base<big_int>::output_stream;
    using io_binary_base<big_int>::input_stream;
    using io_binary_base<big_int>::calc;
    using io_binary_base<big_int>::input_mem;
    using io_binary_base<big_int>::output_mem;

    /// Stores big_int object state to a binary stream. Seft-delimeted binary serialization.
    /** This functions uses the following format:
            SIGN [LEN DIGITS]
        where SIGN is in {-1,0,+1};
        LEN (optional) is a number of limbs used by a given big_int object;
        DIGITS is a dump of all limbs.
        All fields of the format are stored in binary format without any delimeters.
        The optional part LEN DIGITS is emitted only if SIGN != 0.
        This format is The Simple Binary format for big_int.
    */
    template <typename Stream>
    static Stream& output_stream (Stream& out, const big_int& x);


    /// Loads big_int object state from a binary stream. Compatible with output_stream.
    /** See output_stream(stream, big_int) function for detailes on the format.
        If the function fails to read some bytes from a stream, an old value of x
        may be lost (but a given big_int object x remains in the correct state).
        So, do not relay on the value of x when a given stream is not in a good state
        after call returns.

        The function takes input in The Simple Binary format.

        \todo What is a correct way to handle errors in format when we read from a binary stream?
    */
    template <typename Stream>
    static Stream& input_stream (Stream& in, big_int& x);


    /// Calculates the number of chars required to store a given big_int object in The Simple Binary form.
    /** This function calculates precise number of chars that will emit
        any function outputs in The Simple Binary format for one big_int object,
        for example, output_binary_mem function. */
    static std::size_t calc (const big_int& x);
};


/// Compares two big integers
/**
    Returns
    -  0  if a = b,
    -  -1 if a < b,
    -  1  if a > b
*/
int cmp (const big_int & a, const big_int & b);


template <typename TT>
inline int sign (const big_int& x)
{
    return x.sign();
}


/**    @name Standard integer mathematical operations. */
//@{

big_int operator+ (const big_int & b, const big_int & c);
big_int operator- (const big_int & b, const big_int & c);
big_int operator* (const big_int & b, const big_int & c);
big_int operator/ (const big_int & b, const big_int & c);
big_int operator% (const big_int & b, const big_int & c);

big_int operator<< (const big_int& a, std::size_t n);
big_int operator>> (const big_int& a, std::size_t n);

// WARNING.  The following operators have not implemented yet.
inline big_int operator| (const big_int& a, const big_int& b)
{
    ARAGELI_ASSERT_ALWAYS(!"operator| for big_int isn't implemented yet");
    return big_int();
}

inline big_int operator^ (const big_int& a, const big_int& b)
{
    ARAGELI_ASSERT_ALWAYS(!"operato^ for big_int isn't implemented yet");
    return big_int();
}

big_int operator& (const big_int& a, const big_int& b);



/** Requirements: b is representable as std::size_t. */
template <typename T>
inline big_int operator<< (const big_int& a, const T& b)
{
    return a << std::size_t(b);
}

/** Requirements: b is representable as std::size_t. */
template <typename T>
inline big_int operator>> (const big_int& a, const T& b)
{
    return a >> std::size_t(b);
}

template <typename T>
inline big_int& operator<<= (big_int& a, const T& b)
{
    return a = a << b;
}

inline big_int& operator>>= (big_int& a, std::size_t b)
{
    if(a.number.len == 1)
    {
        if(b >= _Internal::bits_per_digit || !(a.number[0] >>= b))
            a = big_int();
    }
    else
        a = a >> b;

    return a;
}

template <typename T>
inline big_int& operator>>= (big_int& a, const T& b)
{
    return a >>= std::size_t(b);
}


//@}



/// Divizion with quotient and remainder.
/** Divides a by b and stores quotient in q and remainder in r.
    @param a a dividend
    @param b a divizor, b must not be a zero
    @param q a quotient
    @param r a remainder  */
inline void divide (const big_int& a, const big_int& b, big_int& q, big_int& r)
{
    return _Internal::xdivide(q, a, b, r);
}


#define _ARAGELI_big_int_CMP1(OPER)    \
    inline bool operator OPER (const big_int& a, const big_int& b)    \
    {    \
        return cmp(a, b) OPER 0;    \
    }

/**    @name Standard comparision operations. */
//@{

_ARAGELI_big_int_CMP1(==)
_ARAGELI_big_int_CMP1(!=)
_ARAGELI_big_int_CMP1(<)
_ARAGELI_big_int_CMP1(>)
_ARAGELI_big_int_CMP1(<=)
_ARAGELI_big_int_CMP1(>=)

//@}


/***************************/
/*                         */
/*   Combined operators    */
/*                         */
/***************************/

#define _ARAGELI_big_int_COMBINED_OPER(OPER)    \
    inline big_int& operator OPER##= (big_int& b, const big_int& c)    \
    {    \
        return b = b OPER c;    \
    }

/**    @name Computed assignment operators. */
//@{

inline big_int& operator += (big_int& b, const big_int& c)
{
    return b = b + c;
}

inline big_int& operator -= (big_int& b, const big_int& c)
{
    return b = b - c;
}

inline big_int& operator *= (big_int& b, const big_int& c)
{
    return b = b * c;
}

inline big_int& operator /= (big_int& b, const big_int& c)
{
    return b = b / c;
}

inline big_int& operator %= (big_int& b, const big_int& c)
{
    return b = b % c;
}

//inline big_int& operator >>= (big_int& b, const big_int& c)
//{
//    return b = b >> c;
//}
//
//inline big_int& operator << (big_int& b, const big_int& c)
//{
//    return b = b << c;
//}

//@}

#define _ARAGELI_big_int_MIXED_AOPER1(NATIVE, AOPER)    \
    inline big_int& operator AOPER (big_int& b, NATIVE c)    \
    {    \
        return b AOPER big_int(c);    \
    }    \
    \
    inline NATIVE& operator AOPER (NATIVE& b, const big_int& c)    \
    {    \
        return b AOPER (NATIVE)c; /* WARNING */    \
    }

#define _ARAGELI_big_int_MIXED_BOPER1(NATIVE, BOPER)    \
    inline big_int operator BOPER (NATIVE b, const big_int& c)    \
    {    \
        return big_int(b) BOPER c;    \
    }    \
    \
    inline big_int operator BOPER (const big_int& b, NATIVE c)    \
    {    \
        return b BOPER big_int(c);    \
    }


#define _ARAGELI_big_int_MIXED_COMARE1(NATIVE, BOPER)    \
    inline bool operator BOPER (NATIVE b, const big_int& c)    \
    {    \
        return big_int(b) BOPER c;    \
    }    \
    \
    inline bool operator BOPER (const big_int& b, NATIVE c)    \
    {    \
        return b BOPER big_int(c);    \
    }


#define _ARAGELI_big_int_MIXED_COMPARE3(NATIVE)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, ==)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, !=)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, <=)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, >=)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, <)    \
    _ARAGELI_big_int_MIXED_COMARE1(NATIVE, >)


#ifdef ARAGELI_INT64_SUPPORT

#define _ARAGELI_big_int_MIXED_INT_AOPER2(AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed __int64, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned __int64, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(bool, AOPER)

#define _ARAGELI_big_int_MIXED_INT_AOPER2(AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed __int64, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned __int64, AOPER)

#define _ARAGELI_big_int_MIXED_INT_BOPER2(BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed short, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned short, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed int, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned int, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed long, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned long, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed __int64, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned __int64, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(bool, BOPER)

#else

#define _ARAGELI_big_int_MIXED_INT_AOPER2(AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(bool, AOPER)

#define _ARAGELI_big_int_MIXED_INT_AOPER2_WO_BOOL(AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned char, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned short, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned int, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(signed long, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(unsigned long, AOPER)

#define _ARAGELI_big_int_MIXED_INT_BOPER2(BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned char, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed short, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned short, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed int, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned int, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(signed long, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(unsigned long, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(bool, BOPER)

#endif

#define _ARAGELI_big_int_MIXED_FLOAT_AOPER2(AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(float, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(double, AOPER)    \
    _ARAGELI_big_int_MIXED_AOPER1(long double, AOPER)

#define _ARAGELI_big_int_MIXED_FLOAT_BOPER2(BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(float, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(double, BOPER)    \
    _ARAGELI_big_int_MIXED_BOPER1(long double, BOPER)


_ARAGELI_big_int_MIXED_INT_AOPER2(+=)
_ARAGELI_big_int_MIXED_INT_AOPER2(-=)
_ARAGELI_big_int_MIXED_INT_AOPER2_WO_BOOL(*=)
_ARAGELI_big_int_MIXED_INT_AOPER2(/=)

_ARAGELI_big_int_MIXED_FLOAT_AOPER2(+=)
_ARAGELI_big_int_MIXED_FLOAT_AOPER2(-=)
_ARAGELI_big_int_MIXED_FLOAT_AOPER2(*=)
_ARAGELI_big_int_MIXED_FLOAT_AOPER2(/=)

_ARAGELI_big_int_MIXED_INT_AOPER2(%=)

//_ARAGELI_big_int_MIXED_INT_AOPER2(<<=)
//_ARAGELI_big_int_MIXED_INT_AOPER2(>>=)

_ARAGELI_big_int_MIXED_INT_BOPER2(+)
_ARAGELI_big_int_MIXED_INT_BOPER2(-)
_ARAGELI_big_int_MIXED_INT_BOPER2(*)
_ARAGELI_big_int_MIXED_INT_BOPER2(/)

_ARAGELI_big_int_MIXED_FLOAT_BOPER2(+)
_ARAGELI_big_int_MIXED_FLOAT_BOPER2(-)
_ARAGELI_big_int_MIXED_FLOAT_BOPER2(*)
_ARAGELI_big_int_MIXED_FLOAT_BOPER2(/)

_ARAGELI_big_int_MIXED_INT_BOPER2(%)
//_ARAGELI_big_int_MIXED_INT_BOPER2(<<)
//_ARAGELI_big_int_MIXED_INT_BOPER2(>>)

//_ARAGELI_big_int_MIXED_INT_BOPER2(==)
//_ARAGELI_big_int_MIXED_INT_BOPER2(!=)
//_ARAGELI_big_int_MIXED_INT_BOPER2(<)
//_ARAGELI_big_int_MIXED_INT_BOPER2(>)
//_ARAGELI_big_int_MIXED_INT_BOPER2(<=)
//_ARAGELI_big_int_MIXED_INT_BOPER2(>=)

//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(==)
//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(!=)
//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(<)
//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(>)
//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(<=)
//_ARAGELI_big_int_MIXED_FLOAT_BOPER2(>=)

_ARAGELI_big_int_MIXED_COMPARE3(char)
_ARAGELI_big_int_MIXED_COMPARE3(signed char)
_ARAGELI_big_int_MIXED_COMPARE3(unsigned char)
_ARAGELI_big_int_MIXED_COMPARE3(signed short)
_ARAGELI_big_int_MIXED_COMPARE3(unsigned short)
_ARAGELI_big_int_MIXED_COMPARE3(signed int)
_ARAGELI_big_int_MIXED_COMPARE3(unsigned int)
_ARAGELI_big_int_MIXED_COMPARE3(signed long)
_ARAGELI_big_int_MIXED_COMPARE3(unsigned long)
_ARAGELI_big_int_MIXED_COMPARE3(bool)

_ARAGELI_big_int_MIXED_COMPARE3(float)
_ARAGELI_big_int_MIXED_COMPARE3(double)
_ARAGELI_big_int_MIXED_COMPARE3(long double)


#ifdef ARAGELI_INT64_SUPPORT

    _ARAGELI_big_int_MIXED_COMPARE3(signed __int64)
    _ARAGELI_big_int_MIXED_COMPARE3(unsigned __int64)

#endif


inline bool big_int::is_null () const
{
    return number.sign == 0;
}

inline bool big_int::is_unit () const
{
    return
        number.len == 1 &&
        number.sign == +1 &&
        number[0] == 1;
}

inline bool big_int::is_opposite_unit () const
{
    return
        number.len == 1 &&
        number.sign == -1 &&
        number[0] == 1;
}


/// Specialization of factory template class for big_int.
template <>
struct factory<big_int>
{
private:
    typedef big_int T;
public:

    /// True iff the functions of this class has a meaning.
    static const bool is_specialized = true;

    /// Unit element (1).
    static const T& unit ()
    {
        static const T unit_s = T(1);
        return unit_s;
    }

    /// Unit element (1) by pattern.
    static const T& unit (const T& x)
    {
        return unit();
    }

    /// Minus unit element (-1).
    static const T& opposite_unit ()
    {
        static const T opposite_unit_s = T(-1);
        return opposite_unit_s;
    }

    /// Minus unit element (-1) by pattern.
    static const T& opposite_unit (const T& x)
    {
        return opposite_unit();
    }

    /// Null element (0).
    static const T& null ()
    {
        static const T null_s;
        return null_s;
    }

    /// Null element by pattern (0).
    static const T& null (const T& x)
    {
        return null();
    }
};


inline bool is_null (const big_int& x)
{
    return x.is_null();
}

inline bool is_unit (const big_int& x)
{
    return x.is_unit();
}

inline bool is_opposite_unit (const big_int& x)
{
    return x.is_opposite_unit();
}

inline bool is_positive (const big_int& x)
{
    return x.sign() > 0;
}

inline bool is_negative (const big_int& x)
{
    return x.sign() < 0;
}

inline bool is_even (const big_int& x)
{
    return x.is_even();
}


inline bool is_odd (const big_int& x)
{
    return x.is_odd();
}


////////////////////////////////////////////////////////////////////////////////
// WARNING! This function prototype is placed here because of header dependences.
// WARNING! The problem appears only if ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE macro is undefined.
// TODO: Resolve this problem and delete.
template <typename T>
T intsqrt (const T& a);
////////////////////////////////////////////////////////////////////////////////


inline big_int sqrt (const big_int& a)
{
    return intsqrt(a);
}

/// Returns absolutly value of x.
inline big_int abs (const big_int& x)
{
    return (x.sign() < 0) ? -x : x;
}

/// Raise x in to a power p.
/** Just uses Arageli::power. */
inline big_int pow (const big_int& x, int p)
{
    return power(x, p);
}

/// Raise x in to a power p.
/** Just uses Arageli::power. */
inline big_int pow (const big_int& x, long int p)
{
    return power(x, p);
}

/// Raise x in to a power p.
/** Just uses Arageli::power. */
inline big_int pow (const big_int& x, const big_int& p)
{
    return power(x, p);
}

/// Returns the base-2 logarithm of abs(x). If x == 0, returns -1.
inline big_int log2 (const big_int& x)
{
    return x.length() - 1;
}


/// Specialization of gcd function for big_int.
/** Allows more efficient computation with small big_int objects
    than generic gcd function. */
template <typename T_factory>
inline big_int gcd (const big_int& a, const big_int& b, const T_factory& tfctr)
{
    if(a.number.len == 1 && b.number.len == 1)
    {
        // As this is Euclid's algorithm for two integers
        // we ignore the signs of the numbers.

        return euclid(a.number[0], b.number[0]);
    }
    else
        return euclid(a, b, tfctr);
}


/// Retruns a number of digits in a given number with specified radix.
std::size_t ndigits (big_int x, std::size_t r);


/// Retruns platform specific magnitude of a given big integer.
/** Exact values returned by this function are implementation defined
    and may differ under different platforms. But this function is satisfied
    the following criteria:

        - magnitude(x) retruns only non negative integer numbers
        for any input x.

        - magnitude(x_1) <= magnitude(x_2) for a given pair x_1 and x_2,
        abs(x_1) <= abs(x_2).

        - C1*nbits(x) - C2 <= magnitude(x) <= C1*nbits(x) + C2,
        where C1 and C2 are platform specific positive constants
        that are the same for any x.
*/
inline std::size_t magnitude (const big_int& x)
{
    return x.number.len;
}

namespace _Internal
{
    template <typename T> struct Unsigned : std::make_unsigned<T> {};
    template <> struct Unsigned<bool> { typedef bool type; };
}

template <typename T>
void big_int::from_native_int_helper (const T &x, true_type)
{
    number = big_struct(1, 1, 1);
    number[0] = _Internal::digit(x);
}

template <typename T>
void big_int::from_native_int_helper (const T &x, false_type)
{

    T xx = x;

    std::size_t n = 0;
    while(xx)
    {
        xx >>= _Internal::bits_per_digit;
        ++n;
    }

    number = big_struct(1, n, n);

    xx = x;

    for(std::size_t i = 0; i < n; ++i, xx >>= _Internal::bits_per_digit)
        number[i] = xx & _Internal::max_digit;

    ARAGELI_ASSERT_1(Arageli::is_null(xx));
}


template <typename T>
void big_int::from_native_int (const T& x)
{
    typedef std::numeric_limits<T> Nl;
    ARAGELI_ASSERT_1(Nl::is_specialized);
    ARAGELI_ASSERT_1(Nl::is_integer);

    if(Arageli::is_null(x))
    {
        number.set_zero();
        return;
    }
    else if(is_negative(x))
    {
        // WARNING. The following expression is not portable.
        from_native_int(static_cast<typename _Internal::Unsigned<T>::type>(-x));
        number.sign = -1;
    }
    else
    {
        typename bool_type<Nl::digits <= _Internal::bits_per_digit>::type dummy;
        from_native_int_helper(x, dummy);
    }
}


template <typename T>
void big_int::from_native_float (const T& x)
{
    typedef std::numeric_limits<T> Nl;

    ARAGELI_ASSERT_0(!(Nl::has_infinity && x == Nl::infinity()));
    ARAGELI_ASSERT_0(!(Nl::has_denorm && x == Nl::denorm_min()));

    if(x > opposite_unit(x) && x < unit(x))
    {
        number.set_zero();
        return;
    }
    else if(is_negative(x))
    {
        from_native_float(-x);
        number.sign = -1;
    }
    else
    {
        T xx = std::floor(x);
        int expon;
        xx = std::frexp(xx, &expon);

        ARAGELI_ASSERT_1(expon > 0);

        xx *= std::pow(2 * unit(x), Nl::digits + 1);
        expon -= Nl::digits + 1;

        ARAGELI_ASSERT_1(std::floor(xx) == xx);

        std::size_t n =
            (Nl::digits + 1) / _Internal::bits_per_digit +
            bool((Nl::digits + 1) % _Internal::bits_per_digit);

        number = big_struct(1, n, n);

        T module = T(_Internal::max_digit);
        module += unit(module);

        for(std::size_t i = 0; i < n; ++i)
        {
            T remaind = std::fmod(xx, module);

            ARAGELI_ASSERT_1(std::floor(remaind) == remaind);
            ARAGELI_ASSERT_1(_Internal::digit(remaind) == remaind);

            number[i] = _Internal::digit(remaind);
            xx = (xx - remaind) / module;

            ARAGELI_ASSERT_1(std::floor(xx) == xx);
        }

        ARAGELI_ASSERT_1(Arageli::is_null(xx));

        big_int::optimize(number);
        ARAGELI_ASSERT_1(number.len != 0);

        if(expon > 0)
            *this <<= expon;
        else
            *this >>= -expon;
    }
}


template <typename T>
T big_int::to_native_int () const
{
    ARAGELI_ASSERT_1(std::numeric_limits<T>::is_specialized);
    ARAGELI_ASSERT_1(std::numeric_limits<T>::is_integer);

    ARAGELI_ASSERT_0(big_int(std::numeric_limits<T>::min()) <= *this && *this <= big_int(std::numeric_limits<T>::max()));

    if(is_null())return factory<T>::null();
    else if(sign() < 0)
    {
        // WARNING. The following expression is not portable.
        return
            -static_cast<T>
            (to_native_int_without_sign<typename _Internal::Unsigned<T>::type>());
    }
    else return to_native_int_without_sign<T>();
}


template <typename T>
T big_int::to_native_int_without_sign () const
{
    ARAGELI_ASSERT_1(std::numeric_limits<T>::is_specialized);
    ARAGELI_ASSERT_1(std::numeric_limits<T>::is_integer);
    ARAGELI_ASSERT_1(*this <= big_int(std::numeric_limits<T>::max()));
    ARAGELI_ASSERT_1(!is_null());

    T res = number.digits()[0];

    for(std::size_t i = 1; i < number.len; ++i)
        res |= number.digits()[i] << (i * _Internal::bits_per_digit);

    return res;
}


template <typename T>
T big_int::to_native_float () const
{
    typedef std::numeric_limits<T> Nl;

    ARAGELI_ASSERT_1(Nl::is_specialized);
    ARAGELI_ASSERT_1(!Nl::is_integer);

    if(is_null())
        return factory<T>::null();
    if(*this < big_int(-Nl::max()))
    {
        ARAGELI_ASSERT_0(Nl::has_infinity);
        return -Nl::infinity();
    }
    else if(*this > big_int(Nl::max()))
    {
        ARAGELI_ASSERT_0(Nl::has_infinity);
        return Nl::infinity();
    }
    else
    {
        big_int t = *this;
        int expon = 0;
        std::size_t blen = t.length();

        if(blen - 1 > Nl::digits)
        {
            ARAGELI_ASSERT_1(size_t(int(blen - 1 - Nl::digits)) == blen - 1 - Nl::digits);
            expon = static_cast<int>(blen - 1 - Nl::digits);
            t >>= expon;
        }

        ARAGELI_ASSERT_1(t.length() - 1 <= Nl::digits);

        T res = factory<T>::null();

        T module = T(_Internal::max_digit);
        module += unit(module);
        T curscale = unit(module);

        for(std::size_t i = 0; i < t.number.len; ++i)
        {
            res += t.number.digits()[i]*curscale;
            curscale *= module;
        }

        res *= std::pow(2*unit(res), expon);
        if(sign() < 0)opposite(&res);

        ARAGELI_DEBUG_EXEC_1(big_int backres = res);
        ARAGELI_ASSERT_1(backres.length() == length());

        ARAGELI_ASSERT_1
        (
            (length() <= Nl::digits && backres - *this == 0)  ||
            (
                length() > Nl::digits &&
                (backres - *this).length() <= length() - Nl::digits
            )
        );

        return res;
    }
}


template <typename Stream>
Stream& io_binary<big_int>::output_stream (Stream& out, const big_int& x)
{
    int sign = x.number.sign;
    output_binary_stream(out, sign);
    if(sign)
    {
        std::size_t len = x.number.len;    // length in limbs
        output_binary_stream(out, len);
        output_binary_stream(out, x.number.digits(), len);
    }

    return out;
}


template <typename Stream>
Stream& io_binary<big_int>::input_stream (Stream& in, big_int& x)
{
    // The following first reads of SIGN and LEN can't break x value.

    int sign;
    if(!input_binary_stream(in, sign))
        return in;
    ARAGELI_ASSERT_ALWAYS(sign == 0 || sign == -1 || sign == +1);

    if(sign)
    {
        // The number isn't zero. Read LEN and DIGITS.

        std::size_t len;
        if(!input_binary_stream(in, len))
            return in;
        ARAGELI_ASSERT_ALWAYS(len > 0);

        if(x.number.len == len)
            x.number.sign = sign;
        else
            x = big_int(sign, len, len);

        // Load DIGITS.
        if(!input_binary_stream(in, x.number.digits(), len))
        {
            // A new value load fails and an old value is lost.
            // Make sure that x object is in correct state.
            x.number[len - 1] = 1;    // kills all leading zeros
        }
    }
    else
    {
        // The number is zero.
        x = big_int();  // WARNING! replace by big_int::assign_null()
    }

    return in;
}


} // namespace Arageli


namespace std
{

inline Arageli::big_int sqrt (const Arageli::big_int& a)
{
    return Arageli::sqrt(a);
}

/// Returns absolutly value of x.
inline Arageli::big_int abs (const Arageli::big_int& x)
{
    return Arageli::abs(x);
}

/// Raise x to a power p.
/** Just uses Arageli::power. */
template <typename T>
inline Arageli::big_int pow (const Arageli::big_int& x, const T& p)
{
    return Arageli::pow(x, p);
}

/// Raise x to a power p.
/** Just uses Arageli::power. */
inline Arageli::big_int pow
(
    const Arageli::big_int& x,
    const Arageli::big_int& p
)
{
    return Arageli::pow(x, p);
}

/// Specialization of std::numeric_limits for Arageli::big_int.
template <> class numeric_limits<Arageli::big_int>
{
    typedef Arageli::big_int TTT;

public:

    static const bool is_specialized = true;

    static TTT min () throw()
    {
        return TTT();
    }

    static TTT max () throw()
    {
        return TTT();
    }

    static const int digits = 0;
    static const int digits10 = 0;
    static const bool is_signed = true;
    static const bool is_integer = true;
    static const bool is_exact = true;
    static const int radix = 2;

    static TTT epsilon () throw()
    {
        return TTT();
    }

    static TTT round_error () throw()
    {
        return TTT();
    }

    static const int min_exponent = 0;
    static const int min_exponent10 = 0;
    static const int max_exponent = 0;
    static const int max_exponent10 = 0;
    static const bool has_infinity = false;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN = false;
    static const float_denorm_style has_denorm = denorm_absent;
    static const bool has_denorm_loss = false;

    static TTT infinity () throw()
    {
        return TTT();
    }

    static TTT quiet_NaN () throw()
    {
        return TTT();
    }

    static TTT signaling_NaN () throw()
    {
        return TTT();
    }

    static TTT denorm_min () throw()
    {
        return TTT();
    }

    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;
    static const bool traps = true;
    static const bool tinyness_before = false;
    static const float_round_style round = round_toward_zero;
};


/// Swaps two numbers without actually copying.
inline void swap (Arageli::big_int& a, Arageli::big_int& b)
{
    a.swap(b);
}


}


namespace Arageli
{

#if 0

/// Specialization of type_pair_traits for big_int and unknown type.
template <typename T2>
struct type_pair_traits<big_int, T2>
{
    static const bool is_specialized = type_traits<T2>::is_specialized;
    static const bool is_convertible = type_traits<T2>::is_number;
    static const bool is_safe_convertible = false;
    static const bool is_assignable = is_convertible;
    static const bool is_safe_assignable = is_safe_convertible;
    static const bool is_initializable = is_convertible;
    static const bool is_safe_initializable = is_safe_convertible;
};


/// Specialization of type_pair_traits for big_int and unknown type.
template <typename T1>
struct type_pair_traits<T1, big_int>
{
    static const bool is_specialized = type_traits<T1>::is_specialized;
    static const bool is_convertible = type_traits<T1>::is_number;

    static const bool is_safe_convertible =
        type_traits<T1>::is_integer_number &&
        std::numeric_limits<T1>::is_specialized &&
        std::numeric_limits<T1>::is_bounded;

    static const bool is_assignable = is_convertible;
    static const bool is_safe_assignable = is_safe_convertible;
    static const bool is_initializable = is_convertible;
    static const bool is_safe_initializable = is_safe_convertible;
};


/// Specialization of type_pair_traits for two big_int.
template <>
struct type_pair_traits<big_int, big_int> :
    public type_pair_traits_for_the_same<big_int> {};

#endif

}



#endif
