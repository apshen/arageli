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

#include "stdafx.hpp"

using namespace Arageli;

TEST_FUNCTION(do_mult_karatsuba, "Test Karatsuba algorithm for multiplication.")
{
    bool test1 = true;
    bool test2 = true;

    try
    {
        int i;
        const unsigned int num_lengths = 50000;
        // generate two random integers approximately with equal lengths
        big_int a = big_int::random_with_length(num_lengths),
                b = big_int::random_with_length(num_lengths);
        _Internal::digit *w_digits = new _Internal::digit[magnitude(a)+magnitude(b)];
        w_digits[magnitude(a)+magnitude(b)-1] = 0;
        _Internal::digit *r_digits = new _Internal::digit[magnitude(a)+magnitude(b)];
        _Internal::digit *t_digits = new _Internal::digit[3*(magnitude(a)+magnitude(b))];
        unsigned w_len =
            do_mult_karatsuba<_Internal::digit,unsigned>(a._digits(),
                    b._digits(), w_digits, t_digits, magnitude(a), magnitude(b));
        // compare karatsuba method result with classic algorithm
        if (_Internal::do_mult_classic(a._digits(), b._digits(), r_digits, magnitude(a), magnitude(b)) != w_len)
        {
            test1 = false;
        }
        for (i = 0; i < w_len; ++i)
        {
            if(w_digits[i] != r_digits[i])
            {
                test1 = false;
            }
        }
        if (!test1)
        {
            tout << "\nThe first karatsuba test failed!\n";
        }
        if (magnitude(a) == magnitude(b))
        {
            w_digits[magnitude(a)*2-1]=0;
            w_len = do_mult_karatsuba<_Internal::digit,unsigned>(w_digits, a._digits(),
                        b._digits(), t_digits, magnitude(a));
            // compare karatsuba method result with classic algorithm
            if (_Internal::do_mult_classic(a._digits(), b._digits(), r_digits, magnitude(a), magnitude(b)) != w_len)
            {
                test2 = false;
            }
            for (i = 0; i < w_len; ++i)
            {
                if(w_digits[i] != r_digits[i])
                {
                    tout << i << '\t' << "diff: " << r_digits[i] - w_digits[i] << '\t';
                    test2 = false;
                }
            }
            if (!test2)
            {
                tout << "The second karatsuba test failed!\n";
            }
        }
        else
        {
            tout << "Generated two numbers with different lengths! Try once more!\n";
        }
        delete [] w_digits;
        delete [] r_digits;
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

    return test1 && test2 ? resOK : resFAIL;
}

