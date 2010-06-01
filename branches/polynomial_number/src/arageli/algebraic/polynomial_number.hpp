/*****************************************************************************

    algebraic/polynomial_number.hpp

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


#ifndef POLYNOMIAL_NUMBER_HPP
#define POLYNOMIAL_NUMBER_HPP

#include <math.h>
#include "../arageli.hpp"
#include "algorithm_class.hpp"
#include "stop_indication_algorithm_class.hpp"
#include "matrix_frobeniuss_form.hpp"


using namespace Arageli;

typedef Arageli::sparse_polynom<big_int > s_p_int;
typedef Arageli::sparse_polynom<rational< big_int> > s_p_rat;


class basis_field
{
private:
    s_p_int BasisPol_m;

    big_int LBasisPol_m; // have not get-method, but have culculate-method
    big_int HBasisPol_m; // have not get-method, but have culculate-method
    rational< big_int> M_m;

    big_int CalculateLBasisPol();
    big_int CalculateHBasisPol();
    rational< big_int> CalculateM(); //use M1 for BasisPol
public:

    basis_field();
    basis_field(const s_p_int & pol);
    basis_field(const char* pol);
    basis_field(const basis_field &BP);

    const s_p_int & BasisPol() const;
    rational< big_int> M() const;
    big_int LBasisPol() const;
    big_int HBasisPol() const;
};

class PolynomialNumber
{
private:
    rational< big_int> L_m; // have not get-method, but have culculate-method
    rational< big_int> H_m; // have not get-method, but have culculate-method
    rational< big_int> F_m;    //
    rational< big_int> U_m;    //
    rational< big_int> SepL_m;    //

    big_int j_1_B;
    big_int j_2_B;
    big_int j_3_B;

    big_int j_1_C;
    big_int j_2_C;
    big_int j_3_C;

    s_p_rat Pol_m;

    rational< big_int> CalculateU();
    rational< big_int> CalculateF(); //de-fact function is return int, cause GetM1 -> int
    rational< big_int> CalculateSepL();

    void Calculate_j_B();
    void Calculate_j_C();
public:
    basis_field * BasisPol;

public:
    PolynomialNumber();
    PolynomialNumber(basis_field & BP);
    PolynomialNumber(basis_field & BP, rational< big_int> m); //monom< rational< big_int>> m);
    PolynomialNumber(const PolynomialNumber &POL);

    PolynomialNumber(rational<big_int> x);
    PolynomialNumber(int x);

    //access to fields in class
    void Pol(const char* pol);
    void Pol(const s_p_rat& pol);
    void Pol(const rational< big_int> x);
    void Pol(int x);
    const s_p_rat & Pol() const;

    rational< big_int> F() const;
    rational< big_int> U() const;
    rational< big_int> SepL() const;

    big_int Get_j_B(int j) const;
    big_int Get_j_C(int j) const;

    //arithmetics operation
    PolynomialNumber operator -() const;
    friend PolynomialNumber abs(PolynomialNumber& POL);

    //PolynomialNumber abs(PolynomialNumber& POL);

    PolynomialNumber & operator = (const PolynomialNumber& POL);
    PolynomialNumber operator +(const PolynomialNumber& POL) const;
    PolynomialNumber operator -(const PolynomialNumber& POL) const;
    PolynomialNumber& operator +=(const PolynomialNumber& POL);
    PolynomialNumber& operator -=(const PolynomialNumber& POL);

    const PolynomialNumber InversePolNum();
    const PolynomialNumber operator *(const PolynomialNumber &POL) const;
    const PolynomialNumber operator /(const PolynomialNumber &POL);
    const PolynomialNumber operator %(const PolynomialNumber &POL) const;
    PolynomialNumber& operator *=(const PolynomialNumber &POL);
    PolynomialNumber& operator /=(const PolynomialNumber &POL);
    PolynomialNumber& operator %=(const PolynomialNumber &POL);

    bool operator ==(const PolynomialNumber &POL) const;
    int operator <(const PolynomialNumber &POL) const;
    int operator >(const PolynomialNumber &POL) const;
    int operator <=(const PolynomialNumber &POL) const;
    int operator >=(const PolynomialNumber &POL) const;

    bool is_positive(const PolynomialNumber &POL) const;
    bool is_negative(const PolynomialNumber &POL) const;
    bool is_null(const PolynomialNumber &POL) const;
    bool is_null() const;
    bool is_unit() const;

    int EquivPolNum(const PolynomialNumber &POL); // +1 - equiv BasisPol & Pol; 0 - equiv BasisPol, not equiv Pol; -1 - other
    rational< big_int> GetLPol();
    rational< big_int> GetHPol();

    int sign() const; //use the sign_binarypolnum_abs_stop_alg< PolynomialNumber> algorithm

    int signPOL_abs(const PolynomialNumber & POL);
    int signPOL_u(const PolynomialNumber & POL);
    int signPOL_sepl(const PolynomialNumber & POL);
    int signPOL_u_sepl(const PolynomialNumber & POL);
    int signPOL_j_1_b(const PolynomialNumber & POL);
    int signPOL_j_2_b(const PolynomialNumber & POL);
    int signPOL_j_3_b(const PolynomialNumber & POL);

    int signCFMPOL_abs(const PolynomialNumber & POL);
    int signCFMPOL_u(const PolynomialNumber & POL);
    int signCFMPOL_sepl(const PolynomialNumber & POL);
    int signCFMPOL_u_sepl(const PolynomialNumber & POL);
    int signCFMPOL_j_1_c(const PolynomialNumber & POL);
    int signCFMPOL_j_2_c(const PolynomialNumber & POL);
    int signCFMPOL_j_3_c(const PolynomialNumber & POL);

    int signums(const s_p_rat f, const s_p_int g); //Shevchenko - Gruzdev method's

};


/// Specialization of common factory template for algebraic number (PolynomialNumber)
template <>
struct factory<PolynomialNumber >
{
private:

    typedef PolynomialNumber T;

public:

    static const bool is_specialized = true;

    static const T& unit ()
    {
        static const T unit_s = T(rational<big_int>(1, 1));
        return unit_s;
    }

    static const T& unit (const T& x)
    {
        static const T unit_s(*(x.BasisPol), rational<big_int>(1, 1)); 
        return unit_s;
    }

    static const T& opposite_unit ()
    {
        static const T opposite_unit_s = T(rational<big_int>(-1, 1));
        return opposite_unit_s;
    }

    static const T& opposite_unit (const T& x)
    {
        static const T opposite_unit_s = T(*(x.BasisPol), rational<big_int>(-1, 1));
        return opposite_unit_s;
    }

    static const T& null ()
    {
        static const T null_s = T(rational<big_int>(0, 1));
        return null_s;
    }


    static const T& null (const T& x)
    {
        static const T null_s(*(x.BasisPol), rational<big_int>(0, 1)); 
        return null_s; 
    }

};


#endif