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


#ifndef ALGEBRAIC_algorithm_class_hpp
#define ALGEBRAIC_algorithm_class_hpp

#include "../sparse_polynom.hpp"
#include "../rational.hpp"
#include "../big_int.hpp"

#include "../std_import.hpp"

namespace Arageli
{

/// BINARY SEARCH METHOD
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



/// CHAIN FRACTION METHOD
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

            while ( (j < 1000) && (!StopPolNumCFM(POL, p_q_j, j)) )
            {
                p_q_j = p_q_2;
                p_q_2 = calc_fraction(pol, p_q_j, p_q_1);
                p_q_1 = p_q_j;
                j++;
            }

            return j * sign(POL.Pol().subs(p_q_j));
        }
    };

// fanction's for method of suitable fractions -------------------------------------------------
    int Dih(const sparse_polynom< big_int> & pol, const rational< big_int> & M) const
    {
        rational< big_int> u(0, 1), v = M + 1;
        rational< big_int> x = Arageli::floor((u+v)/2); //big_int x = big_int((u+v)/2);

        int i = 0;
        while ( x != u ) // x != u && i < 19
        {
            i++;
            if ( pol.subs(x)*pol.subs(u) < 0 ) v = rational< big_int>(x,1);
            else u = rational< big_int>(x,1);
            x = Arageli::floor((u+v)/2); //x = (big_int)(u+v)/2;
        }
        return (int)x;
    };

// ---------------------------------------------------------------------------------------------
    big_int calc_H(const sparse_polynom< big_int> & Pol) const  //Pol
    {
        sparse_polynom< big_int> pol = sparse_polynom< big_int>(Pol);
        big_int cmax = pol.leading_coef();
        for (sparse_polynom< big_int>::coef_iterator ci = pol.coefs_begin(), cj = pol.coefs_end(); ci != cj; ++ci)
            if (cmax < Arageli::abs(*ci) ) cmax = Arageli::abs(*ci);
        return cmax;
    };

// ---------------------------------------------------------------------------------------------
    rational< big_int> calc_M(const sparse_polynom< big_int> & pol) const
    {
        big_int coef = pol.leading_coef();
        return (1 + (rational< big_int>)calc_H(pol)/Arageli::abs(coef));
    };

// fanction's for method of suitable fractions -------------------------------------------------
    void calc_fraction_3(sparse_polynom< big_int> & pol,
                        rational< big_int> & p_q_1,
                        rational< big_int> & p_q_2) const //const
    {
        big_int p0 = 1;
        big_int q0 = 0;
        rational< big_int> M = calc_M(pol);

        int c = Dih(pol, M);
        int n = pol.degree();

        sparse_polynom< big_int> monom_n = sparse_polynom< big_int>::monom(1, n);
        sparse_polynom< big_int> monom_x = sparse_polynom< big_int>::monom(1, -1);

        big_int p1 = c;
        big_int q1 = 1;

        monom_x += sparse_polynom< big_int>::monom(c, 0);
        pol = pol.subs(monom_x);
        pol *= monom_n;

        M = calc_M(pol);
        c = Dih(pol, M);

        p_q_1 = rational< big_int>(p1, q1);
        p_q_2 = rational< big_int>(c*p1 + p0, c*q1 + q0);
    };

// ---------------------------------------------------------------------------------------------
    rational< big_int> calc_fraction(sparse_polynom< big_int> & pol,
                                    const rational< big_int> & p_q_j1,
                                    const rational< big_int> & p_q_j2) const //const
    {
        big_int p_j1 = p_q_j1.numerator();
        big_int q_j1 = p_q_j1.denominator();
        rational< big_int> M = calc_M(pol);
        int n = pol.degree();
        int c = Dih(pol, M);

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

        return rational< big_int>(p_j, q_j);
    };

};


} //- end namespace Arageli --------------------------------------------------------
#endif