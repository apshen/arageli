/*****************************************************************************

    algebraic/stop_indication_algorithm_class.hpp

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


#ifndef ALGEBRAIC_stop_indication_algorithm_class_hpp
#define ALGEBRAIC_stop_indication_algorithm_class_hpp

#include "../config.hpp"
#include "../rational.hpp"
#include "../big_int.hpp"
#include "../sparse_polynom.hpp"
#include "../powerest.hpp"

namespace Arageli
{

// STOP INDICATIONS FOR BINARY SEARCH METHOD

// STOP INDICATION WITH ABS(X) ======================================================================
template <typename T>
class sign_binarypolnum_abs_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(POL.BasisPol->M(),Arageli::power(big_int(2),j)); //SignBinaryPolNum   //+ 1
        if ( Arageli::abs(POL.Pol().subs(x)) >= POL.F()*delta ) return true;                      //+ 1
        return false;
    };

};

template <typename T>
class sign_binarypolnum_abs_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_abs_stop_criteria<T> >
{
};
//=====================================================================================================


// StopPolNum for SignBinaryPolNum 1,B,t.e. U(Pol,BasisPol) >= F*delta ================================
template <typename T>
class sign_binarypolnum_u_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(POL.BasisPol->M(),power((big_int)2,j)); //+ 1,U
        if ( POL.U() >= delta*POL.F() ) return true;   //+ 1,U
        return false;
    };

};

template <typename T>
class sign_binarypolnum_u_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_u_stop_criteria<T> >
{
};
//=====================================================================================================


// StopPolNum for SignBinaryPolNum 3,B,t.e. SepL(Pol,BasisPol) >= delta ===============================
template <typename T>
class sign_binarypolnum_sepl_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(POL.BasisPol->M(),power(big_int(2),j));  //+ 3,S
        if ( delta <= POL.SepL() ) return true;                      //+ 3,S
        return false;
    };

};

template <typename T>
class sign_binarypolnum_sepl_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_sepl_stop_criteria<T> >
{
};
//=====================================================================================================



// StopPolNum for SignBinaryPolNum, delta <= SepL(Pol,BasisPol) || F*delta <= U(Pol,BasisPol) =========
template <typename T>
class sign_binarypolnum_u_sepl_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(POL.BasisPol->M(),power((big_int)2,j));
        if ( (delta <= POL.SepL()) || (POL.U() >= POL.F()*delta) ) return true;
        return false;
    };

};

template <typename T>
class sign_binarypolnum_u_sepl_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_u_sepl_stop_criteria<T> >
{
};
//=====================================================================================================


// StopPolNum for SignBinaryPolNum, j >= j_1_B ========================================================
template <typename T>
class sign_binarypolnum_j_1_b_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_B(1) ) return true;
        return false;
    };

};

template <typename T>
class sign_binarypolnum_j_1_b_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_j_1_b_stop_criteria<T> >
{
};
//=====================================================================================================


// StopPolNum for SignBinaryPolNum, j >= j_2_B ========================================================
template <typename T>
class sign_binarypolnum_j_2_b_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_B(2) ) return true;
        return false;
    };

};

template <typename T>
class sign_binarypolnum_j_2_b_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_j_2_b_stop_criteria<T> >
{
};
//=====================================================================================================


// StopPolNum for SignBinaryPolNum, j >= j_3_B ========================================================
template <typename T>
class sign_binarypolnum_j_3_b_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_B(3) ) return true;
        return false;
    };

};

template <typename T>
class sign_binarypolnum_j_3_b_stop_alg:public sign_binarypolnum_alg<T, sign_binarypolnum_j_3_b_stop_criteria<T> >
{
};
//=====================================================================================================







// STOP INDICATIONS FOR FRACTION SEARCH METHOD

// STOP INDICATION WITH ABS(X) ======================================================================
template <typename T>
class sign_chainfraction_abs_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(1, power(x.denominator(),2)); //SignChainFractionMetodPolNum
        if ( Arageli::abs(POL.Pol().subs(x)) >= POL.F()*delta ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_abs_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_abs_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum 1,B,t.e. U(Pol,BasisPol) >= F*delta =================
template <typename T>
class sign_chainfraction_u_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(1, power(x.denominator(),2)); //SignChainFractionMetodPolNum
        if ( POL.U() >= POL.F()*delta ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_u_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_u_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum 3,B,t.e. SepL(Pol,BasisPol) >= delta ==================
template <typename T>
class sign_chainfraction_sepl_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(1, power(x.denominator(),2)); //SignChainFractionMetodPolNum
        if ( delta <= POL.SepL() ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_sepl_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_sepl_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum, delta <= SepL(Pol,BasisPol) || F*delta <= U(Pol,BasisPol) ====
template <typename T>
class sign_chainfraction_u_sepl_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        rational< big_int> delta = rational< big_int>(1, power(x.denominator(),2)); //SignChainFractionMetodPolNum
        if ( (POL.U() >= POL.F()*delta) || (delta <= POL.SepL()) ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_u_sepl_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_u_sepl_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum, j >= j_1_C ===============================================
template <typename T>
class sign_chainfraction_j_1_c_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_C(1) ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_j_1_c_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_j_1_c_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum, j >= j_2_C ===============================================
template <typename T>
class sign_chainfraction_j_2_c_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_C(2) ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_j_2_c_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_j_2_c_stop_criteria<T> >
{};
//=====================================================================================================



// StopPolNumCFM for SignChainFractionMetodPolNum, j >= j_3_C ===============================================
template <typename T>
class sign_chainfraction_j_3_c_stop_criteria
{
public:
    bool operator() (const T &POL, const rational< big_int> &x, int j)
    {
        if ( j >= POL.Get_j_C(3) ) return true;        //+ 1
        return false;
    };

};

template <typename T>
class sign_chainfraction_j_3_c_stop_alg:public sign_chainfraction_alg<T, sign_chainfraction_j_3_c_stop_criteria<T> >
{};
//=====================================================================================================

} //- end namespace Arageli --------------------------------------------------------
#endif