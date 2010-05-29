/*****************************************************************************

    algebraic/algorithm_class.hpp

    This file is a part of the Arageli library.

    Copyright (C) 2010 Natalia Klemenova

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


#ifndef ALGORITHM_CLASS_H
#define ALGORITHM_CLASS_H

#include "../arageli.hpp"

using namespace Arageli;

#if 1

// BINARY SEARCH METHOD
template <typename T, typename StopCriteria>
class sign_binarypolnum_alg
{
public:

    int operator ()(const T &POL)const
    {
        StopCriteria StopPolNum;

        if ( POL.Pol().is_null() ) return 0;
        else
        {
            int j = 0;
            rational< big_int> v = rational< big_int>(POL.BasisPol->M(), 1);
            rational< big_int> u = 0;
            rational< big_int> x = (u+v)/big_int(2);
            while ( (j < 1000) && (!StopPolNum(POL, x, j)) )
            {
                if ( sign(POL.BasisPol->BasisPol().subs(x)) == sign(POL.BasisPol->BasisPol().subs(v)) ) v = x;
                else u = x;
                j++;
                x = (u + v)/((big_int)2);
            }
            return sign(POL.Pol().subs(x)); // j *  //j - will be deleted; for debuging info
        }
    };
};

//bool is_it_out = true; //extern

// CHAIN FRACTION METHOD
template <typename T, typename StopCriteria>
class sign_chainfraction_alg
{
public:

    int operator ()(const T &POL) const
    {
        StopCriteria StopPolNumCFM;

        if ( POL.Pol().is_null() ) return 0;
        else
        {
            sparse_polynom< big_int> pol = sparse_polynom< big_int>(POL.BasisPol->BasisPol());
            int j = 1;

            rational< big_int> p_q_1 = 0;
            rational< big_int> p_q_2 = 0;

            calc_fraction_3(pol, p_q_1, p_q_2);
            rational< big_int> p_q_j = p_q_1;

#if 0
            std::cout << "p_q_1 = " << p_q_1 << std::endl;
            std::cout << "p_q_2 = " << p_q_2 << std::endl;
            std::cout << "p_q_j = " << p_q_j << std::endl;
#endif

            while ( (j < 1000) && (!StopPolNumCFM(POL, p_q_j, j)) )
            {
                p_q_j = p_q_2;
                p_q_2 = calc_fraction(pol, p_q_j, p_q_1);
                p_q_1 = p_q_j;
                j++;

#if 0
                std::cout << "j = " << j << std::endl;
                std::cout << "p_q_1 = " << p_q_1 << std::endl;
                std::cout << "p_q_2 = " << p_q_2 << std::endl;
                std::cout << "p_q_j = " << p_q_j << std::endl;
#endif
            }

            return j * sign(POL.Pol().subs(p_q_j));
        }
    };

// fanction's for method of suitable fractions -------------------------------------------------
    int Dih(const sparse_polynom< big_int> & pol, const rational< big_int> & M) const
    {

#if 0
        std::cout << std::endl << "Dih" << std::endl;
#endif

        rational< big_int> u(0, 1), v = M + 1;
        rational< big_int> x = Arageli::floor((u+v)/2); //big_int x = big_int((u+v)/2);

#if 0
        std::cout << "u = " << u << "  v = " << v << "  x = " << x << std::endl;
#endif

        int i = 0;
        while ( x != u ) // x != u && i < 19
        {
            i++;
            if ( pol.subs(x)*pol.subs(u) < 0 ) v = rational< big_int>(x,1);
            else u = rational< big_int>(x,1);
            x = Arageli::floor((u+v)/2); //x = (big_int)(u+v)/2;

#if 0
            std::cout << "i = " << i << "  u = " << u << "  v = " << v << "  x = " << x << std::endl;
#endif

        }
        return (int)x;
    };

// ---------------------------------------------------------------------------------------------
    big_int calc_H(const sparse_polynom< big_int> & Pol) const  //Pol
    {

#if 0
        std::cout << std::endl << "calc_H" << std::endl;
#endif

        sparse_polynom< big_int> pol = sparse_polynom< big_int>(Pol);
        big_int cmax = pol.leading_coef();
        for (sparse_polynom< big_int>::coef_iterator ci = pol.coefs_begin(), cj = pol.coefs_end(); ci != cj; ++ci)
            if (cmax < Arageli::abs(*ci) ) cmax = Arageli::abs(*ci);

#if 0
        std::cout << "pol = " << pol << "  cmax = " << cmax << std::endl << std::endl;
#endif

        return cmax;
    };

// ---------------------------------------------------------------------------------------------
    rational< big_int> calc_M(const sparse_polynom< big_int> & pol) const
    {

#if 0
        std::cout << std::endl << "calc_M" << std::endl;
#endif

        big_int coef = pol.leading_coef();

#if 0
        std::cout << "coef = " << coef << std::endl;
#endif

        return (1 + (rational< big_int>)calc_H(pol)/Arageli::abs(coef));
    };

// fanction's for method of suitable fractions -------------------------------------------------
    void calc_fraction_3(sparse_polynom< big_int> & pol,
                        rational< big_int> & p_q_1,
                        rational< big_int> & p_q_2) const //const
    {

#if 0
        std::cout << std::endl << "calc_fraction_3" << std::endl;
#endif

        big_int p0 = 1;
        big_int q0 = 0;
        rational< big_int> M = calc_M(pol);

#if 0
        std::cout << "M = " << M << std::endl;
#endif

        int c = Dih(pol, M);

#if 0
        std::cout << "c = " << c << std::endl;
#endif

        int n = pol.degree();

#if 0
        std::cout << "n = " << n << std::endl;
        std::cout << "p0 = " << p0 << "  q0 = " << q0 << "  pol = " << pol << "  M = " << M << "  c = " << c << "  n = " << n << std::endl;
#endif

        sparse_polynom< big_int> monom_n = sparse_polynom< big_int>::monom(1, n);
        sparse_polynom< big_int> monom_x = sparse_polynom< big_int>::monom(1, -1);

        big_int p1 = c;
        big_int q1 = 1;

        monom_x += sparse_polynom< big_int>::monom(c, 0);
        pol = pol.subs(monom_x);
        pol *= monom_n;

#if 0
        std::cout << "pol = " << pol << std::endl;
#endif

        M = calc_M(pol);

#if 0
        std::cout << "M = " << M << std::endl;
#endif

        c = Dih(pol, M);

#if 0
        std::cout << "c = " << c << std::endl;
        std::cout << "p1 = " << p1 << "  q1 = " << q1 << "  pol = " << pol << "  M = " << M << "  c = " << c << std::endl;
#endif

        p_q_1 = rational< big_int>(p1, q1);
        p_q_2 = rational< big_int>(c*p1 + p0, c*q1 + q0);

#if 0
        std::cout << "p_q_1 = " << p_q_1 << "  p_q_2 = " << p_q_2 << std::endl << std::endl;
#endif

    };

// ---------------------------------------------------------------------------------------------
    rational< big_int> calc_fraction(sparse_polynom< big_int> & pol,
                                    const rational< big_int> & p_q_j1,
                                    const rational< big_int> & p_q_j2) const //const
    {
#if 0
        std::cout << std::endl << "calc_fraction" << std::endl;
#endif

        big_int p_j1 = p_q_j1.numerator();
        big_int q_j1 = p_q_j1.denominator();
        rational< big_int> M = calc_M(pol);
        int n = pol.degree();
        int c = Dih(pol, M);

#if 0
        std::cout << "p_j1 = " << p_j1 << "  q_j1 = " << q_j1 << "  pol = " << pol << "  M = " << M << "  c = " << c << "  n = " << n << std::endl;
#endif

        sparse_polynom< big_int> monom_n = sparse_polynom< big_int>::monom(1, n);
        sparse_polynom< big_int> monom_x = sparse_polynom< big_int>::monom(1, -1);
        big_int p_j2 = p_q_j2.numerator();
        big_int q_j2 = p_q_j2.denominator();

        monom_x += sparse_polynom< big_int>::monom(c,0);
        pol = pol.subs(monom_x);
        pol *= monom_n;
        M = calc_M(pol);
        c = Dih(pol,M);
        big_int p_j = c*p_j1 + p_j2;
        big_int q_j = c*q_j1 + q_j2;

#if 0
        std::cout << "p_j2 = " << p_j2 << "  q_j2 = " << q_j2 << "  pol = " << pol << "  M = " << M << "  c = " << c << std::endl;
        std::cout << "p_j = " << p_j << "  q_j = " << q_j << std::endl << std::endl;
#endif

        return rational< big_int>(p_j, q_j);
    };


    /*
    //fanction's for method of suitable fractions ---------------------------
    void calc_p_q(big_int* p, big_int* q, rational< big_int>* p_q)
    {
        p_q[0] = rational< big_int>(1,1);
        for (int i = 1; i < 200; i++ )
        {
            p_q[i] = rational< big_int>(p[i],q[i]);
        }
    }
    //extern big_int* p;
    //extern big_int* q;
    //
    //extern rational< big_int>* p_q;
    */
};

#endif


#endif