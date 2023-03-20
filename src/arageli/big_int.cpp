/*****************************************************************************

    big_int.cpp -- see the big_int.hpp file for description.

    This file is a part of the Arageli library.

    WARNIG. This file has no complate implementation.

    Copyright (C) 1999 -- 2005, 2008 Nikolai Yu. Zolotykh

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
    \file big_int.cpp
    \brief The big_int.hpp file stuff implementation.
*/

#include "config.hpp"

#include <cstdlib>
#include <malloc.h>
#include <sstream>
#include <limits>
#include <cctype>

#include "big_int.hpp"
#include "rational.hpp"


namespace
{
typedef Arageli::_Internal::digit digit;
}


namespace Arageli
{

void big_arith_error(const char *s)
{
    throw big_int::exception(std::string("Big arith error: ") + s);
}


/***************************/
/*                         */
/*    low level memory     */
/*   managment routines    */
/*                         */
/***************************/

void big_int::copy_data
(digit* dest, const digit* source, std::size_t newnitems)
{
    memmove(dest, source, newnitems * sizeof(digit));
}

/***************************/
/*                         */
/*      calc_bdn_radix     */
/*                         */
/***************************/

void calc_bdn_radix (digit radix, digit& bdn_radix, std::size_t& chars_per_digit)
{
    //bdn_radix = maximal power of radix fit in digit;

    bdn_radix = 1;
    chars_per_digit = 0;
    digit t = _Internal::max_digit/radix;

    while(t)
    {
        t /= radix;
        bdn_radix *= radix;
        chars_per_digit++;
    }
}


/***************************/
/*                         */
/*         optimize        */
/*                         */
/***************************/

void big_int::optimize (big_struct &number)
{
    // deleting leading zeros

    std::size_t new_len = _Internal::do_optimize(number.digits(), number.len);

    if(new_len != number.len)
    {
        if(new_len)
        {
            big_struct new_num(number.sign, new_len, new_len);
            copy_data(new_num.digits(), number.digits(), new_len);
            number = std::move(new_num);
        }
        else
        {
            number.set_zero();
        }
    }
}

/***************************/
/*                         */
/*      constructors       */
/*    and destructors      */
/*                         */
/***************************/

big_int::big_int (const char* str)
{
    std::istringstream s(str);
    big_int b;
    s >> b;
    if(!s && !s.eof())
        throw incorrect_string(str);

    *this = b;
}


/***************************/
/*                         */
/*        operator =       */
/*                         */
/***************************/


big_int& big_int::operator= (const big_int & b)
{
    ARAGELI_ASSERT_1(b.number.sign == -1 || b.number.sign == 0 || b.number.sign == 1);

    if(b.number.sign == 0)
    {
        number.set_zero();
    }
    else
    {
        number = big_struct(b.number.sign, b.number.len, b.number.len);
        copy_data(number.digits(), b.number.digits(), b.number.len);
    }

    return *this;
}

/***************************/
/*                         */
/*      Unary operators    */
/*                         */
/***************************/


big_int big_int::operator- () const    // unary minus
{
    if(number.sign != 0)
    {
        big_int a(-number.sign, number.len, number.len);
        copy_data(a.number.digits(), number.digits(), number.len);
        return a;
    }

    return big_int();
}

/***************************/
/*                         */
/*        operator +       */
/*                         */
/***************************/

big_int operator+ (const big_int& b, const big_int& c)
{
    ARAGELI_ASSERT_1(b.number.sign == -1 || b.number.sign == 0 || b.number.sign == 1);
    ARAGELI_ASSERT_1(c.number.sign == -1 || c.number.sign == 0 || c.number.sign == 1);

    big_int a;
    int bsign = b.number.sign;
    int csign = c.number.sign;
    std::size_t blen = b.number.len;
    std::size_t clen = c.number.len;
    digit *bdata = b.number.digits();
    digit *cdata = c.number.digits();
    const big_int::big_struct *u, *v;

    if(bsign == 0)
        a = c;
    else if(csign == 0)
        a = b;
    else
    {
        if(bsign == csign)
        {
            if (blen >= clen)
            {
                u = &b.number;
                v = &c.number;
            }
            else
            {
                u = &c.number;
                v = &b.number;
            }
            big_int::big_struct sum(bsign, u->len + 1, u->len + 1);
            big_int::copy_data(sum.digits(), u->digits(), u->len);
            sum.len = _Internal::do_add(sum.digits(), v->digits(), u->len, v->len);
            a.number = std::move(sum);
        }
        else // bsign != csign
        {
            if(blen >= clen)
            {
                big_int::big_struct sum(bsign, blen, blen);
                big_int::copy_data(sum.digits(), bdata, blen);
                bool was_borrow = _Internal::do_sub(sum.digits(), cdata, blen, clen);
                if(was_borrow)
                {
                    big_int::copy_data(sum.digits(), cdata, clen);
                    _Internal::do_sub(sum.digits(), bdata, clen, blen);
                    sum.sign = csign;
                }
                big_int::optimize(sum);
                a.number = std::move(sum);
            }
            else
            {
                big_int::big_struct sum(csign, clen, clen);
                big_int::copy_data(sum.digits(), cdata, clen);
                _Internal::do_sub(sum.digits(), bdata, clen, blen);
                big_int::optimize(sum);
                a.number = std::move(sum);
            }
        }
    }

    return a;
}


/***************************/
/*                         */
/*        operator -       */
/*                         */
/***************************/


big_int operator- (const big_int& b, const big_int& c)
{
    return b + (-c);
}


/***************************/
/*                         */
/*        operator *       */
/*                         */
/***************************/


big_int operator* (const big_int& b, const big_int& c)
{
    ARAGELI_ASSERT_1(b.number.sign == -1 || b.number.sign == 0 || b.number.sign == 1);
    ARAGELI_ASSERT_1(c.number.sign == -1 || c.number.sign == 0 || c.number.sign == 1);

    int bsign = b.number.sign;
    int csign = c.number.sign;

    if((bsign != 0) && (csign != 0))
    {
        std::size_t blen = b.number.len;
        std::size_t clen = c.number.len;

        std::size_t resultlen = blen + clen;
        big_int a(bsign * csign, resultlen, resultlen);

        resultlen = _Internal::do_mult
        (
            b.number.digits(),
            c.number.digits(),
            a.number.digits(),
            blen,
            clen
        );

        a.number.len = resultlen;
        return a;
    }
    else
        return big_int();
}


/***************************/
/*                         */
/*         xdivide          */
/*                         */
/***************************/


void _Internal::xdivide (big_int& a, const big_int& b, const big_int& c, big_int& res)
{
    // a = b / c; res = b % c

    digit runint;

    std::size_t alen;
    std::size_t rlen;
    std::size_t blen = b.number.len;
    std::size_t clen = c.number.len;
    int bsign = b.number.sign;
    int csign = c.number.sign;

    if(csign == 0)
        big_arith_error("divide by zero");
    else if(bsign == 0)
    {
        a = big_int();
        res = big_int();
    }
    else if (blen < clen)
    {
        a = big_int();
        res = b;
    }
    else if (clen == 1)
    {
        big_int::big_struct q(bsign * csign, blen, blen);
        runint =
            _Internal::do_divide_by_digit
            (
                b.number.digits(),
                q.digits(),
                blen,
                c.number[0]
            );

        if(q[blen - 1])
            alen = blen;
        else
            alen = blen - 1;

        if(alen)
        {
            a.number = std::move(q);
            a.number.len = alen;
        }
        else
            a.number.set_zero();

        if(runint != 0)
        {
            res = big_int(bsign, 1, 1);
            res.number[0] = runint;
        }
        else
            res.number.set_zero();

    }
    else
    {
        big_int::big_struct u(1, blen + 1, blen + 1);
        big_int::big_struct v(1, clen, clen);
        big_int::big_struct q(bsign * csign, blen - clen + 1, blen - clen + 1);
        big_int::copy_data(u.digits(), b.number.digits(), blen);
        big_int::copy_data(v.digits(), c.number.digits(), clen);

        rlen = _Internal::do_divide(u.digits(), v.digits(), q.digits(), blen, clen);
        if(q[blen - clen])
            alen = blen - clen + 1;
        else
            alen = blen - clen;

        if(rlen == 0)
            res.number.set_zero();
        else
        {
            res = big_int(bsign, rlen, rlen);
            big_int::copy_data(res.number.digits(), u.digits(), rlen);
        }

        if(alen != 0)
        {
            a.number = std::move(q);
            a.number.len = alen;
        }
        else
        {
            a.number.set_zero();
        }
    }

    ARAGELI_ASSERT_1(b == c*a + res);
}


big_int operator% (const big_int& b, const big_int& c)
{
    big_int a, r;
    _Internal::xdivide(a, b, c, r);
    return r;
}


big_int operator/ (const big_int& b, const big_int& c)
{
    big_int a, r;
    _Internal::xdivide(a, b, c, r);
    return a;
}


/***************************/
/*                         */
/*    Compare functions    */
/*                         */
/***************************/


int cmp (const big_int & a, const big_int & b)
{
    //sign(a-b)

    int result;

    int asign, bsign;
    std::size_t alen, blen;
    big_int c;

    asign = a.number.sign;
    bsign = b.number.sign;
    alen = a.number.len;
    blen = b.number.len;

    if(asign < bsign)
        result = -1;
    else if(asign > bsign)
        result = 1;
    else if(asign == 0)
        result = 0;            // asign == bsign == 0
    else if(alen < blen)
        result = -asign;    // asign == bsign != 0
    else if(alen > blen)
        result = asign;
    else
    { // asign == bsign != 0, alen == blen
        c = a - b;
        result = c.number.sign;
    }

    return result;
}


/***************************/
/*                         */
/*      Random numbers     */
/*                         */
/***************************/

digit random_digit ()
{
    if(_Internal::max_digit <= RAND_MAX)return std::rand();

    static const std::size_t rand_max_len = big_int(RAND_MAX).length();    // WARNING! It is not efficient.
    const std::size_t n =
        _Internal::bits_per_digit / rand_max_len +
        bool(_Internal::bits_per_digit % rand_max_len);

    digit res = 0;
    for(std::size_t i = 0; i < n; ++i)
        (res <<= rand_max_len) |= std::rand();    // WARNING! Lower bits from rand is
                                            // placed to higher bits of the result.
    return res;
}

big_int big_int::random_in_range (const big_int& max)
{
    if(max.is_null())
        return max;

    std::size_t len = max.length();

    for(;;)
    {
        big_int res = random_with_length_or_less(len);
        if(res <= max)
            return res;
    }
}


big_int big_int::random_with_length_or_less (std::size_t length) //bits
{
    int bits_in_highest = length % _Internal::bits_per_digit;
    std::size_t len = length/_Internal::bits_per_digit + 1;

    big_int a(1, len, len);

    for(std::size_t i = 0; i < len - 1; i++)
        a.number[i] = random_digit();

    if(bits_in_highest)
        a.number[len - 1] = random_digit() >>
           (_Internal::bits_per_digit - bits_in_highest);
    else
        a.number[len - 1] = 0;

    big_int::optimize(a.number);

    return a;
}

/***************************/
/*                         */
/*     Stream operators    */
/*                         */
/***************************/

digit stream_radix(std::ios & s)
{
    if(s.flags() & std::ios::dec)
        return 10;
    if(s.flags() & std::ios::hex)
        return 16;
    if(s.flags() & std::ios::oct)
        return 8;
    else
        return 10; // it is very strange
}


std::ostream& operator<< (std::ostream& s, const big_int& x)
{
    if(s.flags() & std::ios_base::showpos && x.sign() > 0)
        s << '+';

    if(!x.number.sign)return s << "0";

    digit bdn_radix;
    std::size_t chars_per_block;
    calc_bdn_radix(stream_radix(s), bdn_radix, chars_per_block);

    std::size_t bdnlen = x.number.len;
    big_int::big_struct bdn(1, 2 * bdnlen, 2 * bdnlen); //bdnlen + bdnlen/chars_per_block

    big_int::big_struct numberdata(1, x.number.len, x.number.len);

    big_int::copy_data(numberdata.digits(), x.number.digits(), x.number.len);
    bdnlen = _Internal::do_big_int_to_bdn(numberdata.digits(), bdn.digits(), x.number.len, bdn_radix);

    _Internal::auto_stream_state _ass(s, s.flags() & ~std::iostream::showpos);

    if(x.number.sign == -1)
        s << '-';
    s << bdn[bdnlen - 1];

    for (std::size_t i = bdnlen - 1; i > 0; i--)
    {
        s.width(chars_per_block);
        s.fill('0');
        s << bdn[i - 1];
    }

    return s;
}


inline void set_stream_radix (std::ios& s, digit radix)
{
    switch(radix)
    {
        case 10:
            s.setf(std::ios::dec, std::ios::basefield);
            break;
        case 16:
            s.setf(std::ios::hex, std::ios::basefield);
            break;
        case 8:
            s.setf(std::ios::oct, std::ios::basefield);
            break;
        default:
            throw big_int::exception("Value for radix is invalid");
    }
}


std::istream& operator>> (std::istream& s, big_int& x)
{
#if 0    // old version

    // WARNING!!! This function is incorrect!

    char ch;
    int sign = 1;
    digit radix = 10;
    bool brackets = false;

    while(s.get(ch) && isspace(ch));    // pass leading spaces

    if(!s)
        big_arith_error("empty string for reading as big int");

    if(ch == '-')
    {
        sign = -1;
        s.get(ch);
    }
    else if(ch == '+')
        s.get(ch);
    else if(ch == '(')
    {
        s.get(ch);
        brackets = true;
    }

    if(!s)big_arith_error("invalid string format for reading as big int");

    // now in ch the first digit of number notation

    int zero = 0;
    if (ch == '0')
    {
        zero = 1;
        s.get(ch);
        switch (ch)
        {
            case 'o':
                radix = 8;
                break;
            case 'x':
                radix = 16;
                break;
            default:
                s.putback(ch);
        }
    }
    else
        s.putback(ch);

    do    // pass non significant first zeros
    {
        s.get(ch);
    }while (ch == '0');
    s.putback(ch);

    std::ostringstream buffer;
    std::size_t tellp = 0;

    while(s.get(ch))
    {
        if
        (
            radix == 10 && isdigit(ch) ||
            radix == 16 && isxdigit(ch) ||
            radix == 8  && isdigit(ch) && ch < '8'
        )
        {
            buffer << ch;
            ++tellp;
        }
        else
        {
            s.putback(ch);
            break;
        }
    }

    if(tellp)
    {
        if(!zero)
            big_arith_error("empty string for reading as big int");
        x = 0;
        return s;
    }

    std::string buffer_cpp_str = buffer.str();
    const char* buffer_str = buffer_cpp_str.c_str();

    digit bdn_radix;
    std::size_t chars_per_block;
    calc_bdn_radix(radix, bdn_radix, chars_per_block);

    std::size_t length_of_str = buffer_cpp_str.length();
    std::size_t number_of_blocks = (length_of_str - 1)/chars_per_block + 1;
    digit* bdn = big_int::get_mem_for_data(number_of_blocks);
    digit* data = 0;
    std::size_t len = 0;

    try
    {
        std::size_t length_of_block = length_of_str % chars_per_block;
                                // for the first block
        if(!length_of_block)
            length_of_block = chars_per_block;

        const char * buffer_str_cur = buffer_str;
        for(std::size_t j = number_of_blocks; j > 0; j--)
        {
            std::istringstream
                stream_block
                    (std::string(buffer_str_cur).substr(0, length_of_block));
            set_stream_radix(stream_block, radix);
            stream_block >> bdn[j-1];
            buffer_str_cur += length_of_block;
            length_of_block = chars_per_block;  // for the last blocks
        }

        data = big_int::get_mem_for_data(number_of_blocks);
        len = _Internal::do_bdn_to_big_int(bdn, data, number_of_blocks, bdn_radix);
    }
    catch(...)
    {
        big_int::free_data(bdn); big_int::free_data(data);
        throw;
    }

    big_int::free_data(bdn);

    try
    {
        data = big_int::realloc_data(data, len);

        if(data)
            x.free_mem_and_alloc_number(sign, data, len);
        else
            x.free_mem_and_alloc_zero();
    }
    catch(...)
    {
        big_int::free_data(data);
        throw;
    }

    if(brackets)
    {
        s >> ch;
        if(ch != ')')
            big_arith_error("invalid string format for reading as big int");
    }

    s.clear();
    return s;

#else

    char ch = 0;
    int sign = 1;
    digit radix = 10;
    bool brackets = false;

    std::size_t len;

    //do s.get(ch); while (isspace(ch));
    s >> ch;    // pass leading spaces (if this feature is turn on for s)
                // and read first non space character
    if(!s)
    {
        s.clear(std::ios_base::badbit);
        return s;
    }

    if(ch == '(')
    {
        brackets = true;
        s >> ch;
    }

    if(ch == '-')sign = -1;
    else if(ch != '+')s.putback(ch);

    s.get(ch);

    int zero = 0;
    if(ch == '0')
    {
        zero = 1;
        s.get(ch);

        switch (ch)
        {
            case 'o':
                radix = 8;
                break;
            case 'x':
                radix = 16;
                break;
            default:
                s.putback(ch);
        }
    }
    else s.putback(ch);

    s.clear();
    do
    {
        s.get(ch);
    }while (ch == '0' && s);

    if(s)
        s.putback(ch);
    else
    {
        x = big_int();
        s.clear();
        return s;
    }

    std::string buffer;

    while(s.get(ch))
    {
        if
        (
            (radix == 10 && isdigit(ch)) ||
            (radix == 16 && isxdigit(ch)) ||
            (radix == 8  && isdigit(ch) && ch < '8')
        )
            buffer += ch;
        else
        {
            s.putback(ch);
            break;
        }
    }

    if(buffer.empty())
    {
        if(!zero)
            throw big_int::incorrect_string(std::string(1, s.peek()));

        if(brackets && !_Internal::read_literal(s, ")"))
            throw big_int::incorrect_string(std::string(1, s.peek()));

        x = big_int();
        s.clear();
        return s;
    }

    digit bdn_radix;
    std::size_t chars_per_block;
    calc_bdn_radix(radix, bdn_radix, chars_per_block);

    std::size_t length_of_str = buffer.length();
    std::size_t number_of_blocks = (length_of_str - 1)/chars_per_block + 1;
    big_int::big_struct bdn(1, number_of_blocks, number_of_blocks);
    std::size_t length_of_block = length_of_str % chars_per_block;
                            // for the first block
    if(!length_of_block)
        length_of_block = chars_per_block;

    std::size_t buffer_cur = 0;
    for(std::size_t j = number_of_blocks; j > 0; j--)
    {
        std::string block = buffer.substr(buffer_cur, length_of_block);
        std::istringstream stream_block(block);
        set_stream_radix(stream_block, radix);
        stream_block >> bdn[j-1];
        buffer_cur += length_of_block;
        length_of_block = chars_per_block;  // for the last blocks
    }

    big_int::big_struct data(1, number_of_blocks, number_of_blocks);

    len = _Internal::do_bdn_to_big_int
    (
        bdn.digits(),
        data.digits(),
        number_of_blocks,
        bdn_radix
    );

    if(len)
    {
        x = big_int(sign, len, len);
        big_int::copy_data(x.number.digits(), data.digits(), len);
    }
    else
        x.number.set_zero();

    if(brackets && !_Internal::read_literal(s, ")"))
        throw big_int::incorrect_string(std::string(1, s.peek()));

    s.clear();
    return s;

#endif
}


std::size_t big_int::length () const
{
    if(!number.len)
        return 0;
    std::size_t l = (number.len - 1) * _Internal::bits_per_digit;
    digit highest = number.digits()[number.len - 1];
    while(highest >>= 1)
        l++;
    return l + 1;
}


bool big_int::operator[] (std::size_t k) const
{
    ARAGELI_ASSERT_0(k < length());
    return (number.digits()[k / _Internal::bits_per_digit] >>
        (k % _Internal::bits_per_digit)) % 2;
}

big_int operator<< (const big_int& a, std::size_t n)
{
    std::size_t a_len = a.number.len;
    if(!n || !a_len)
        return a;
    std::size_t l = a.length() + n;
    std::size_t res_len = l/_Internal::bits_per_digit;
    if (l%_Internal::bits_per_digit)
        res_len++;
    std::size_t m = n/_Internal::bits_per_digit;
    n %= _Internal::bits_per_digit;
    digit t = 0;
    std::size_t k = _Internal::bits_per_digit - n;
    big_int res(a.number.sign, res_len, res_len);
    const digit* data = a.number.digits();

    std::fill_n(res.number.digits(), m, digit(0));

    if(n)
    {
        for(std::size_t i = 0; i < a_len; i++)
        {
            res.number[m + i] = (data[i] << n) | t;
            t = data[i] >> k;
        }
        if(t)
            res.number[m + a_len] = t;
    }
    else
    {
        std::copy(data, data + a_len, res.number.digits() + m);
    }

    return res;
}


big_int operator>> (const big_int& a, std::size_t n)
{
    std::size_t a_len = a.number.len;
    if (!n || !a_len)
        return a;
    std::size_t l = a.length();
    if (l <= n)
        return big_int();
    l -= n;
    std::size_t res_len = l/_Internal::bits_per_digit;
    if (l % _Internal::bits_per_digit)
        res_len++;
    std::size_t m = n / _Internal::bits_per_digit;
    n %= _Internal::bits_per_digit;
    std::size_t k = _Internal::bits_per_digit - n;
    digit t = 0;

    big_int res(a.number.sign, res_len, res_len);
    const digit* data = a.number.digits();

    if(n)
    {
        if(a_len > res_len + m)
            t = data[a_len - 1] << k;
        for(std::size_t i = res_len; i > 0; i--)
        {
            res.number[i - 1] = (data[i + m - 1] >> n) | t;
            t = data[i + m - 1] << k;
        }
    }
    else
    {
        std::copy(data + m, data + m + res_len, res.number.digits());
    }

    return res;
}


big_int operator& (const big_int& a, const big_int& b)
{
    std::size_t reslen = std::min(a.number.len, b.number.len);
    if(reslen == 0)
        return big_int();
    big_int::big_struct res(1, reslen, reslen);

    // the current size of nonzero part of the number in res
    std::size_t lastnz = 0;

    for(std::size_t i = 0; i < reslen; ++i)
    {
        res[i] = a.number.digits()[i] & b.number.digits()[i];
        if(res[i])
            lastnz = i+1;
    }

    if(lastnz)
    {
        big_int resbi(1, lastnz, lastnz);
        big_int::copy_data(resbi.number.digits(), res.digits(), lastnz);
        return resbi;
    }
    else
    {
        return big_int();
    }
}


std::size_t ndigits (big_int x, std::size_t r)
{
    ARAGELI_ASSERT_0(r >= 2);

    if(is_negative(x))
        opposite(&x);   // Is it really necessary?

    std::size_t res = 0;
    big_int br = r; // as we do not have optimized operations with the built-ins
    while(!is_null(x))
    {
        x /= br;
        ++res;
    }

    return res;
}


std::size_t io_binary<big_int>::calc (const big_int& x)
{
    if(x.number.sign)
        return
            sizeof(int) +
            sizeof(std::size_t) +
            x.number.len * sizeof(big_int::digit);
    else
        return sizeof(int);
}


}
