/*****************************************************************************

    karatsuba.cpp

    This file is a part of the Arageli library.

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
    \brief The karatsuba.hpp file stuff implementation
*/


#include "config.hpp"
#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) ||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_karatsuba)
#include "karatsuba.hpp"
#include "bigar.hpp"
#include "_utility.hpp"
#include <iostream>

namespace Arageli
{

/**
Preliminary conditions: w[n+m-1] = 0!
                        m >= n
*/
template <typename N,typename T>
T do_mult_karatsuba(const N *u, const N *v, N *w, N *t, T m, T n)
{
    ARAGELI_ASSERT_0(m>=n);

    T k = m >> 1;

    if (!is_null(k) && n > k && n > ARAGELI_KARATSUBA_THRESHOLD)
    {
        T k2 = k << 1;
        N *W0 = w+k2;
        N *W1 = w+k;
        N *W2 = w;
        const N *U0 = u+k;
        const N *U1 = u;
        const N *V0 = v+k;
        const N *V1 = v;

        T Clen, C1len, C2len, UV0len, UV1len;

        // allocate temprary space
        // N *C1 = new N[k+2];
        N *C1 = t;
        // N *C2 = new N[k+2];
        N *C2 = t+k+2;
        N *C = t+m+4;
        t += 2*m+6;

        UV0len = do_mult_karatsuba(U0, V0, W0, t, m-k, n-k); // mult U0*V0
        UV1len = do_mult_karatsuba(U1, V1, W2, t, k, k);     // mult U1*V1

        memcpy(C1, U0, sizeof(N)*(m-k));
        if (_Internal::do_add(C1, U1, m-k, k))              // U0+U1
        {
            C1len = m - k + 1;
            C1[m-k] = 1;
        }
        else
        {
            C1len = m - k; // U0+U1
        }
        if (n>=k2)    // V0+V1
        {
            // here n = m or n = m-1 and (n-k) = k or (n-k) = k+1
            memcpy(C2, V0, sizeof(N)*(n-k));
            if(_Internal::do_add(C2, V1, n-k, k))           // V0+V1
            {
                C2len = n - k + 1;
                C2[n-k] = 1;
            }
            else
            {
                C2len = n - k;
            }
        }
        else
        {
            memcpy(C2, V1, sizeof(N)*k);
            if (_Internal::do_add(C2, V0, k, n-k))           // V0+V1
            {
                C2len = k + 1;
                C2[k] = 1;
            }
            else
            {
                C2len = k;
            }
        }

        if (C2len > C1len)              // (U0+U1)*(V0+V1)
        {
            Clen = do_mult_karatsuba(C2, C1, C, t, C2len, C1len); // mult C = C1*C2
        }
        else
        {
            Clen = do_mult_karatsuba(C1, C2, C, t, C1len, C2len); // mult C = C1*C2
        }

        if (_Internal::do_sub(C, W0, Clen, UV0len))     // U1*V1+U0*V1+U1*V0
        {
            C[Clen] = 1;
        }
        if (_Internal::do_sub(C, W2, Clen, UV1len))     // U0*V1+V0*U1
        {
            C[Clen] = 1;
        }

        if (_Internal::do_add(W1, C, m+n-k, Clen))      // W += C
        {
            W1[m+n-k] = 1;
        }

        t -= 2*m+6;
    }
    else
    {
        // no sence to use karatsuba algorithm in this case
        return _Internal::do_mult_classic(u, v, w, m, n);
    };

    return (w[m+n-1] != 0) ? (m+n): (m+n-1);
};

/**
 * Preliminary conditions: w[2*n-1] = 0!
 */
// This implementation was taken from GMP library.
template <typename N,typename T>
T do_mult_karatsuba(N *r, const N *u, const N *v, N *t, T n)
{
    ARAGELI_ASSERT_0(n >= 1);

    T i;
    N w, w0, w1;
    const N *x, *y;
    T n2 = n >> 1;  // floor(n/2);
    ARAGELI_ASSERT_0(n2 > 0);
    int sign = 0;

    if (n & 1)
    {
        // Odd length
        // Computing u0-u1. Result stores in r.
        T n3 = n - n2;
        w = u[n2];
        if (w)
        {   // here u0 is greater than u1, i.e. (u0 - u1) > 0
            w -= _Internal::do_sub(r, u, u+n3, n2, n2);
        }
        else
        {
            i = n2;
            do{
                --i;
                w0 = u[i];
                w1 = u[n3+i];
            }while ((w0 == w1) && (i != 0));
            if (w0 < w1)
            {   // here u1 is greater than u0, so compute (u1 - u0) and store sign.
                x = u+n3;
                y = u;
                sign = ~0;
            }
            else
            {
                x = u;
                y = u+n3;
            }
            _Internal::do_sub(r, x, y, n2, n2);
        }
        r[n2] = w;

        // Computing v0-v1. Result stores in r+n3.
        w = v[n2];
        if (w)
        {   // here v0 is greater than v1, i.e. (v0 - v1) > 0
            w -= _Internal::do_sub(r+n3, v, v+n3, n2, n2);
        }
        else
        {
            i = n2;
            do{
                --i;
                w0 = v[i];
                w1 = v[n3+i];
            }while ((w0 == w1) && (i != 0));
            if (w0 < w1)
            {   // here v1 is greater than v0, so compute (v1 - v0) and store sign.
                x = v+n3;
                y = v;
                sign = ~sign;
            }
            else
            {
                x = v;
                y = v+n3;
            }
            _Internal::do_sub(r+n3, x, y, n2, n2);
        }
        r[n] = w;

        T n1 = n+1;
        if (n2 < ARAGELI_KARATSUBA_THRESHOLD)
        {
            if (n3 < ARAGELI_KARATSUBA_THRESHOLD)
            {
                _Internal::do_mult_classic(r, r+n3, t, n3, n3);
                _Internal::do_mult_classic(u, v, r, n3, n3);
            }
            else
            {
                // TODO:Check if we get into this branch!!
                do_mult_karatsuba(t, r, r+n3, t+n1, n3);
                do_mult_karatsuba(r, u, v, t+n1, n3);
            }
            _Internal::do_mult_classic(u+n3, v+n3, r+n1, n2, n2);
        }
        else
        {
            do_mult_karatsuba(t, r, r+n3, t+n1, n3);
            do_mult_karatsuba(r, u, v, t+n1, n3);
            do_mult_karatsuba(r+n1, u+n3, v+n3, t+n1, n2);
        }
        if (sign)
        {
            _Internal::do_add(t, r, n1, n1);
        }
        else
        {
            _Internal::do_sub(t, r, t, n1, n1);
        }

        T nm1 = n-1;
        if (_Internal::do_add(t, r+n1, nm1, nm1))
        {
            N x = (t[nm1] + 1) & _Internal::max_digit;
            t[nm1] = x;
            if (x == 0)
            {
                t[n] = (t[n] + 1) & _Internal::max_digit;
            }
        }
        if (_Internal::do_add(r+n3, t, n1, n1))
        {
            int cy = 1;
            for (int i = 0; cy; ++i)
            {
                N t = r[n1 + n3 + i];
                r[n1 + n3 + i] += cy;
                cy = r[n1 + n3 + i] < t;
            }
        }
    }
    else
    {
        // Even length
        i = n2;
        do
        {
            --i;
            w0 = u[i];
            w1 = u[n2+i];
        }
        while ((w0 == w1) && (i != 0));
        if (w0 < w1)
        {
            x = u+n2;
            y = u;
            sign = ~0;
        }
        else
        {
            x = u;
            y = u+n2;
        }
        _Internal::do_sub(r, x, y, n2, n2);

        i = n2;
        do
        {
            --i;
            w0 = v[i];
            w1 = v[n2+i];
        }
        while ((w0 == w1) && (i != 0));
        if (w0 < w1)
        {
            x = v+n2;
            y = v;
            sign = ~sign;
        }
        else
        {
            x = v;
            y = v+n2;
        }
        _Internal::do_sub(r+n2, x, y, n2, n2);

        if (n2 < ARAGELI_KARATSUBA_THRESHOLD)
        {
            _Internal::do_mult_classic(r, r+n2, t, n2, n2);
            _Internal::do_mult_classic(u, v, r, n2, n2);
            _Internal::do_mult_classic(u+n2, v+n2, r+n, n2, n2);
        }
        else
        {
            do_mult_karatsuba(t, r, r+n2, t+n, n2);
            do_mult_karatsuba(r, u, v, t+n, n2);
            do_mult_karatsuba(r+n, u+n2, v+n2, t+n, n2);
        }

        if (sign)
        {
            w = _Internal::do_add(t, r, n, n);
        }
        else
        {
            w = -1*_Internal::do_sub(t, r, t, n, n);
        }
        w += _Internal::do_add(t, r+n, n, n);
        w += _Internal::do_add(r+n2, t, n, n);
        int i = 0;
        do
        {
            N t = r[n + n2 + i];
            r[n + n2 + i] += w;
            w = r[n + n2 + i] < t;
            ++i;
        }while(w && (i < 2*n - n + n2));
    }
    return (r[2*n-1]) ? 2*n : 2*n - 1;
};

}
#endif