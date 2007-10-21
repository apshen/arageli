/*****************************************************************************

    test/karatsuba.cpp

    This file is a part of the Arageli library test base.

    Copyright (C) 1999--2006 Nikolai Yu. Zolotykh
    Copyright (C) 2006 Aleksey Bader
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
    \file karatsuba.cpp
    \brief This file includes test for karatsuba multiplication algorithm.

    This test multiplies two random integers approximately with equal lengths
    by karatsuba method and compare result with trivial algorithm.
*/

//#define CHECK_TIME

#include "stdafx.hpp"

using namespace Arageli;

TEST_FUNCTION(do_mult_karatsuba, "Test Karatsuba algorithm for multiplication.")
{
    bool is_ok = true;

    try
    {
        int i;
        const unsigned int num_lengths = 200000;
        // generate two random integers approximately with equal lengths
        big_int a = big_int::random_with_length(num_lengths),
            b = big_int::random_with_length(num_lengths);
        // compute digit length
        const unsigned digit_len = sizeof(_Internal::digit)*8;
        unsigned long mask = 0;             // mask for one digit
        for (i = 0; i < digit_len; ++i)
        {
            mask |= 1 << i;
        }
        unsigned a_len = 0;
        unsigned b_len = 0;
        // conpute number of digits in a and b
        big_int t = a;
        for (; t != 0; ++a_len)
        {
            t >>= digit_len;
        }
        t = b;
        for (; t != 0; ++b_len)
        {
            t >>= digit_len;
        }
        _Internal::digit *a_digits = new _Internal::digit[a_len];
        _Internal::digit *b_digits = new _Internal::digit[b_len];
        t = a;
        for (i = 0; i < a_len; ++i)
        {
            a_digits[i] = _Internal::digit(t & big_int(mask));
            t >>= digit_len;
        }
        t = b;
        for (i = 0; i < b_len; ++i)
        {
            b_digits[i] = _Internal::digit(t & big_int(mask));
            t >>= digit_len;
        }
        // alloc memory for output number
        big_int out = 0;
        _Internal::digit *w_digits = new _Internal::digit[a_len+b_len];
        _Internal::digit *t_digits = new
            _Internal::digit[9*(a_len+b_len)];
        w_digits[a_len+b_len-1] = 0;
#ifdef CHECK_TIME
        timer s;
        s.start();
#endif
        unsigned w_len =
            do_mult_karatsuba<_Internal::digit,unsigned>(a_digits,
                    b_digits, w_digits, &t_digits, a_len, b_len);
#ifdef CHECK_TIME
        s.stop();
        std::cerr << "Karatsuba time: " << s.time() << "\n";
#endif
        out = 0;
        for (i = 0; i < w_len; ++i)
        {
            out += big_int(w_digits[i])<<(i*digit_len);
        }
        // compare result with trivial multiplication
        if (out != a*b)
        {
            is_ok = false;
        }
        delete [] a_digits;
        delete [] b_digits;
        delete [] w_digits;
        delete [] t_digits;
    }
    catch(const Arageli::exception& e)
    {
        tout << e;
        return resEXCEPT;
    }
    catch(const std::exception& e)
    {
        tout << e.what();
        return resEXCEPT;
    }

    return is_ok ? resOK : resFAIL;

    return resOK;
}

