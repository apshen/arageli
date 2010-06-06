/*****************************************************************************

    algebraic/polynomial_number.cpp

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


#include "polynomial_number.hpp"
#include "../gauss.hpp"

namespace Arageli
{

/// class basis_field

basis_field::basis_field():
    M_m(1,1), LBasisPol_m(0), HBasisPol_m(0)
{
    BasisPol_m = sparse_polynom< big_int>(0);
}

basis_field::basis_field(const s_p_int & pol)
{
    BasisPol_m = sparse_polynom< big_int>(pol);
    LBasisPol_m = CalculateLBasisPol();
    HBasisPol_m = CalculateHBasisPol();
    M_m = CalculateM();
}

basis_field::basis_field(const char *pol)
{
    BasisPol_m = sparse_polynom< big_int>(pol);
    LBasisPol_m = CalculateLBasisPol();
    HBasisPol_m = CalculateHBasisPol();
    M_m = CalculateM();
}

basis_field::basis_field(const basis_field &BP)
{
    BasisPol_m = BP.BasisPol();
    LBasisPol_m = BP.LBasisPol();
    HBasisPol_m = BP.HBasisPol();
    M_m = BP.M();
}

const s_p_int & basis_field::BasisPol() const
{
    return BasisPol_m;
}

rational< big_int> basis_field::M() const
{
    return M_m;
}

big_int basis_field::LBasisPol() const
{
    return LBasisPol_m;
}

big_int basis_field::HBasisPol() const
{
    return HBasisPol_m;
}

big_int basis_field::CalculateLBasisPol()
{
    int csum = 0;
    s_p_int pol = BasisPol();
    for (s_p_int::coef_iterator ci = pol.coefs_begin(), cj = pol.coefs_end(); ci != cj; ++ci)
        csum += abs(*ci);
    return csum;
}

big_int basis_field::CalculateHBasisPol()
{
    s_p_int pol = BasisPol();
    big_int cmax = pol.leading_coef();
    for (s_p_int::coef_iterator ci = pol.coefs_begin(), cj = pol.coefs_end(); ci != cj; ++ci)
        if (cmax < abs(*ci) ) cmax = abs(*ci);
    return cmax;
}

rational< big_int> basis_field::CalculateM()
{
    big_int coef = BasisPol_m.leading_coef();
    return (1 + (rational< big_int>)HBasisPol()/abs(coef));
}

//Metod's
PolynomialNumber::PolynomialNumber():
    L_m(0, 1), H_m(0, 1), F_m(0, 1), U_m(0, 1), SepL_m(0, 1), j_1_B(0), j_2_B(0), j_3_B(0), j_1_C(0), j_2_C(0), j_3_C(0)
{
    BasisPol = NULL;
    Pol_m = "1"; // TODO make this Pol_m = null();
}

PolynomialNumber::PolynomialNumber(basis_field & BP):
    L_m(0, 1), H_m(0, 1), F_m(0, 1), U_m(0, 1), SepL_m(0, 1), j_1_B(0), j_2_B(0), j_3_B(0), j_1_C(0), j_2_C(0), j_3_C(0)
{
    BasisPol = &BP;
    Pol_m = "1"; // TODO make this Pol_m = null();
}

PolynomialNumber::PolynomialNumber(basis_field & BP, rational< big_int> m): //monom< rational< big_int>> m)
    L_m(m), H_m(m), F_m(0, 1), SepL_m(0, 1), j_1_B(0), j_2_B(0), j_3_B(0), j_1_C(0), j_2_C(0), j_3_C(0)
{
    BasisPol = &BP;

    Pol(m);
    U_m = CalculateU();
}

PolynomialNumber::PolynomialNumber(const PolynomialNumber &POL)
{
    if (BasisPol != POL.BasisPol)
    {
        BasisPol = POL.BasisPol;
    }
    Pol(POL.Pol()); 
    F_m = POL.F();
    U_m = POL.U();
    SepL_m = POL.SepL();
    j_1_B = POL.j_1_B;
    j_2_B = POL.j_2_B;
    j_3_B = POL.j_3_B;
    j_1_C = POL.j_1_C;
    j_2_C = POL.j_2_C;
    j_3_C = POL.j_3_C;
}

PolynomialNumber::PolynomialNumber(rational<big_int> x):
    L_m(x), H_m(x), F_m(0, 1), U_m(0, 1), SepL_m(0, 1), j_1_B(0), j_2_B(0), j_3_B(0), j_1_C(0), j_2_C(0), j_3_C(0)
{
    BasisPol = NULL;
    Pol_m = s_p_rat::monom(x, 0);
}

PolynomialNumber::PolynomialNumber(int x):
    L_m(0, 1), H_m(0, 1), F_m(0, 1), U_m(0, 1), SepL_m(0, 1), j_1_B(0), j_2_B(0), j_3_B(0), j_1_C(0), j_2_C(0), j_3_C(0)
{
    BasisPol = NULL;
    Pol_m = s_p_rat::monom(rational< big_int>(x, 1), 0);
}


//access operation to fields of class
//template <typename P>
void PolynomialNumber::Pol(const char* pol)
{
    Pol_m = sparse_polynom< rational< big_int>>(pol);

    if ( !Pol_m.is_null() )
    {
        if ( BasisPol != NULL )
        {
            s_p_rat baspol = BasisPol->BasisPol();
            Pol_m = Pol_m % baspol;
        }

        U_m = CalculateU();
        F_m = CalculateF();
        SepL_m = rational< big_int>(0,1);    //SepL_m = CalculateSepL(); // amplifier
        j_1_B = 0; j_2_B = 0; j_3_B = 0;    //Calculate_j_B(); // amplifier
        j_1_C = 0; j_2_C = 0; j_3_C = 0;    //Calculate_j_C(); //amplifier
    }
    else
    {
        U_m = rational< big_int>(0,1);
        F_m = rational< big_int>(0,1);
        SepL_m = rational< big_int>(0,1);
        j_1_B = 0; j_2_B = 0; j_3_B = 0;
        j_1_C = 0; j_2_C = 0; j_3_C = 0;
    }
}

void PolynomialNumber::Pol(const s_p_rat & pol)
{
    Pol_m = sparse_polynom< rational< big_int>>(pol);

    if ( !Pol_m.is_null() )
    {
        if ( BasisPol != NULL )
        {
            s_p_rat baspol = BasisPol->BasisPol();
            Pol_m = Pol_m % baspol;
        }
        /* */
        U_m = CalculateU();
        F_m = CalculateF();
        SepL_m = rational< big_int>(0,1);    //SepL_m = CalculateSepL(); // amplifier
        j_1_B = 0; j_2_B = 0; j_3_B = 0;    //Calculate_j_B(); // amplifier
        j_1_C = 0; j_2_C = 0; j_3_C = 0;    //Calculate_j_C(); //amplifier
    }
    else
    {
        U_m = rational< big_int>(0,1);
        F_m = rational< big_int>(0,1);
        SepL_m = rational< big_int>(0,1);
        j_1_B = 0; j_2_B = 0; j_3_B = 0;
        j_1_C = 0; j_2_C = 0; j_3_C = 0;
    }
}

void PolynomialNumber::Pol(const rational< big_int> x)
{
    Pol_m = sparse_polynom< rational< big_int>>::monom(x, 0);

    if ( !Pol_m.is_null() )
    {
        if ( BasisPol != NULL )
        {
            s_p_rat baspol = BasisPol->BasisPol();
            Pol_m = Pol_m % baspol;
        }
        /* */
        U_m = CalculateU();
        F_m = CalculateF();
        SepL_m = rational< big_int>(0,1);    //SepL_m = CalculateSepL(); // amplifier
        j_1_B = 0; j_2_B = 0; j_3_B = 0;    //Calculate_j_B(); // amplifier
        j_1_C = 0; j_2_C = 0; j_3_C = 0;    //Calculate_j_C(); //amplifier
    }
    else
    {
        U_m = rational< big_int>(0,1);
        F_m = rational< big_int>(0,1);
        SepL_m = rational< big_int>(0,1);
        j_1_B = 0; j_2_B = 0; j_3_B = 0;
        j_1_C = 0; j_2_C = 0; j_3_C = 0;
    }
}

void PolynomialNumber::Pol(int x)
{
    Pol_m = sparse_polynom< rational< big_int>>::monom(x, 0);

    if ( x != 0 ) //!Pol_m.is_null()
    {
        U_m = CalculateU();
        F_m = CalculateF();
        SepL_m = rational< big_int>(0,1);    //SepL_m = CalculateSepL(); // amplifier
        j_1_B = 0; j_2_B = 0; j_3_B = 0;    //Calculate_j_B(); // amplifier
        j_1_C = 0; j_2_C = 0; j_3_C = 0;    //Calculate_j_C(); //amplifier
    }
    else
    {
        U_m = rational< big_int>(0,1);
        F_m = rational< big_int>(0,1);
        SepL_m = rational< big_int>(0,1);
        j_1_B = 0; j_2_B = 0; j_3_B = 0;
        j_1_C = 0; j_2_C = 0; j_3_C = 0;
    }
}

const s_p_rat & PolynomialNumber::Pol() const
{
    return Pol_m;
}

rational< big_int> PolynomialNumber::F() const
{
    return F_m;
}

rational< big_int> PolynomialNumber::U() const
{
    return U_m;
}

rational< big_int> PolynomialNumber::SepL() const
{
    return SepL_m;
}

big_int PolynomialNumber::Get_j_B(int j) const
{
    switch (j)
    {
    case 1:
        {
            return j_1_B;
            break;
        }
    case 2:
        {
            return j_2_B;
            break;
        }
    case 3:
        {
            return j_3_B;
            break;
        }
    default: return 0; break;
    }
}

big_int PolynomialNumber::Get_j_C(int j) const
{
    switch (j)
    {
    case 1:
        {
            return j_1_C;
            break;
        }
    case 2:
        {
            return j_2_C;
            break;
        }
    case 3:
        {
            return j_3_C;
            break;
        }
    default: return 0; break;
    }
}

PolynomialNumber PolynomialNumber::operator -() const
{
    PolynomialNumber S(*this);
    s_p_rat pol = this->Pol();
    pol = (-pol);
    S.Pol(pol);
    return S;
}

PolynomialNumber abs(const PolynomialNumber &POL) 
{
    PolynomialNumber S(POL);
    if ( S.is_positive(S) ) return S;
    else
    {
        S.Pol(-POL.Pol());
        return S;
    }

}

PolynomialNumber& PolynomialNumber::operator =(const PolynomialNumber &POL)
{
    if ( this->BasisPol != POL.BasisPol && POL.BasisPol != NULL)
    {
        this->BasisPol = POL.BasisPol;
        // could be redefinition for LBasisPol_m, HBasisPol_m, M_m
    }
    if ( this->Pol() != POL.Pol() )
    {
        this->Pol(POL.Pol()); //this->Pol = POL.Pol; // % (s_p_rat)BasisPol;
        // could be redefinition for L_m, H_m, U_m, F_m and other
    }
    return *this;
}

PolynomialNumber PolynomialNumber::operator +(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        PolynomialNumber sumPol(POL);
        sumPol += (*this);
        return sumPol;
    }
}

PolynomialNumber PolynomialNumber::operator -(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        PolynomialNumber subPol(*this);
        subPol -= POL;
        return subPol;
    }
}

PolynomialNumber& PolynomialNumber::operator +=(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        s_p_rat sum;
        sum = Pol();
        if ( POL.BasisPol != NULL )
        {
            *this = POL;
        }
        sum += POL.Pol();
        this->Pol(sum);
        return *this;
    }
    else return *this;
}

PolynomialNumber& PolynomialNumber::operator -=(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        s_p_rat dim = Pol();
        dim -= POL.Pol();
        if ( POL.BasisPol != NULL )
        {
            *this = POL;
        }
        this->Pol(dim);
        return *this;
    }
    else return *this;
}

const PolynomialNumber PolynomialNumber::InversePolNum()
{
    if (BasisPol != NULL )
    {
        if ( gcd(Pol(),BasisPol->BasisPol()) == 1 )
        {
            PolynomialNumber InvPol(*BasisPol);
            s_p_rat pq, U, V, B = BasisPol->BasisPol();
            pq = euclid_bezout(Pol(), B, U, V);
            InvPol.Pol(U);
            return InvPol;
        }
    }
    else
    {
        PolynomialNumber InvPol(*BasisPol);
        if ( Pol().is_const() )
        {
            rational< big_int> num = Pol().leading_coef();
            InvPol.Pol(inverse(num));
        }
        else
        {
            s_p_rat num = inverse(Pol());
            InvPol.Pol(num);
        }

        return InvPol;
    }
}

const PolynomialNumber PolynomialNumber::operator *(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL ) // I
    {
        PolynomialNumber MultPol(*this);
        MultPol *= POL;
        return MultPol;
    }
}

const PolynomialNumber PolynomialNumber::operator /(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )  // I
    {
        PolynomialNumber DivPol(*BasisPol);
        PolynomialNumber InvPol(POL);
        PolynomialNumber Numerator(*this);

        InvPol = InvPol.InversePolNum();
        DivPol = Numerator * InvPol;
        return DivPol;
    }
}

const PolynomialNumber PolynomialNumber::operator %(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        if ( BasisPol != NULL)
        {
            /*PolynomialNumber pol = null(*this);
            return pol;*/
            return null(this); //factory< PolynomialNumber>::
        }
        else
        {
            /*PolynomialNumber pol = null(POL);
            return pol;*/
            return null(POL);
        }
    }
}

PolynomialNumber& PolynomialNumber::operator *=(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL ) // I
    {
        s_p_rat mult = Pol();
        if (POL.BasisPol != NULL)
        {
            *this = POL;
        }
        mult *= POL.Pol();
        this->Pol(mult);
        return *this;
    }
    else return *this;
}

PolynomialNumber& PolynomialNumber::operator /=(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL ) // I
    {
        s_p_rat div = Pol();
        if ( BasisPol != NULL )
        {
            PolynomialNumber S(*BasisPol);
            S.Pol(POL.Pol());
            S = S.InversePolNum();
            div *= S.Pol();
        }
        else
        {
            PolynomialNumber S(POL);
            S = S.InversePolNum();
            div *= S.Pol();
        }
        this->Pol(div);
        return *this;
    }
    else return *this;
}

PolynomialNumber& PolynomialNumber::operator %=(const PolynomialNumber &POL)
{
    *this = null(POL);
    return *this;
}

bool PolynomialNumber::operator ==(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || POL.BasisPol == NULL || BasisPol == NULL )
    {
        if ( Pol() == POL.Pol() ) return true;
        return false;
    }
    else return false;
}

bool PolynomialNumber::operator !=(const PolynomialNumber &POL) const
{
    if ( BasisPol == POL.BasisPol || BasisPol == NULL || POL.BasisPol == NULL )
    {
        if ( Pol() == POL.Pol() ) return false;
        return true;
    }
    return true;
}

int PolynomialNumber::operator <(const PolynomialNumber &POL) const
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    if ( s(POL - *this) > 0 ) return 1;
    else return 0;
}

int PolynomialNumber::operator >(const PolynomialNumber &POL) const
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    if ( s(*this - POL) > 0 ) return 1;
    else return 0;
}

int PolynomialNumber::operator <=(const PolynomialNumber &POL) const
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    if ( s(POL - *this) >= 0 ) return 1;
    else return 0;
}

int PolynomialNumber::operator >=(const PolynomialNumber &POL) const
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    if ( s(*this - POL) >= 0 ) return 1;
    else return false;
}

bool PolynomialNumber::is_positive(const PolynomialNumber &POL) const //
{
    if ( POL.BasisPol != NULL )
    {
        if ( POL.sign() > 0 ) return true; //null<sparse_polynom<rational<>>>("0")
        else return false;
    }
    else
    {
        rational< big_int> coef = POL.Pol().leading_coef();
        if ( Arageli::sign(coef) > 0 ) return true;
        return false;
    }
}

bool PolynomialNumber::is_negative(const PolynomialNumber &POL) const
{
    if ( POL.BasisPol != NULL )
    {
        if ( POL.sign() < 0 ) return true;
        else return false;
    }
    else
    {
        rational< big_int> coef = POL.Pol().leading_coef();
        if ( Arageli::sign(coef) < 0 ) return true;
        return false;
    }
}

bool PolynomialNumber::is_null(const PolynomialNumber &POL) const
{
    if ( POL.BasisPol != NULL )
    {
        if ( POL.sign() == 0 ) return true;
        return false;
    }
    else
    {
        if ( POL.Pol().is_null() ) return true; // == null(POL)
        return false;
    }
}

bool PolynomialNumber::is_null() const
{
    if ( this->BasisPol != NULL )
    {
        if ( this->sign() == 0 ) return true;
        return false;
    }
    else
    {
        if ( this->Pol().is_null() ) return true; // == null(POL)
        return false;
    }
}

bool PolynomialNumber::is_unit() const
{
    return this->Pol().is_unit();
}

/**/

// +1 - equiv BasisPol & Pol;
// 0 - equiv BasisPol, not equiv Pol;
// -1 - other: not equiv BasisPol & Pol
int PolynomialNumber::EquivPolNum(const PolynomialNumber &POL)
{
    if ( BasisPol == POL.BasisPol)
    {
        if ( Pol() == POL.Pol() ) return 1;
        else return 0;
    }
    else return -1;
}

rational< big_int> PolynomialNumber::GetLPol()
{
    rational< big_int> csum = 0;
    s_p_rat pol_iterator(Pol());
    for (s_p_rat::coef_iterator ci = pol_iterator.coefs_begin(), cj = pol_iterator.coefs_end(); ci != cj; ++ci)
        csum += Arageli::abs(*ci);
    return csum;
}

rational< big_int> PolynomialNumber::GetHPol()
{
    rational< big_int> cmax = Pol().leading_coef();
    s_p_rat pol_iterator(Pol());
    for (s_p_rat::coef_iterator ci = pol_iterator.coefs_begin(), cj = pol_iterator.coefs_end(); ci != cj; ++ci)
        if (cmax < Arageli::abs(*ci) ) cmax = Arageli::abs(*ci);
    return cmax;
}

rational< big_int> PolynomialNumber::CalculateU()
{
    if (BasisPol != NULL)
    {
        rational< big_int> UPol = 0;
        int n = BasisPol->BasisPol().degree() - 1;
        int k = Pol().degree();
        rational< big_int> temp = Arageli::power(GetLPol(),n) * Arageli::power(BasisPol->LBasisPol(),k); //GetLPol() -> L(); GetLBasisPol()
        UPol = temp.inverse();
        return UPol;
    }
    else return rational< big_int>(1, 100);
}

rational< big_int> PolynomialNumber::CalculateF()
{
    rational< big_int> F = 0;
    int degree = Pol().degree();
    if ( degree > 1 && BasisPol != NULL)
    {
        for(int i = 2; i <= degree; i++)
        {
            F += (i)*power(BasisPol->M(),i-1); //F += (i)*power(GetM(),i-1);
        }
        F *= GetHPol();
    }
    F += GetHPol();
    return F;
}

rational< big_int> PolynomialNumber::CalculateSepL()
{
    if (BasisPol != NULL)
    {
        //int n = BasisPol->BasisPol().degree();
        int k = Pol().degree();
        int n_k = k + BasisPol->BasisPol().degree();
        rational< big_int> temp1 = rational< big_int>(1, power(n_k, (n_k >>= 1) + 1));  // (1,Arageli::sqrt((big_int)power(n+k,n+k+2)));
        rational< big_int> temp2_1 = rational< big_int>(BasisPol->LBasisPol() * GetLPol().numerator() * pow(2, k), GetLPol().denominator()); //pow(2, k);
        rational< big_int> temp2_2 = Arageli::power(temp2_1, n_k--);
        rational< big_int> temp2 = rational< big_int>(temp2_2.denominator(), temp2_2.numerator()); // 1,power(BasisPol->LBasisPol()*GetLPol()*pow(2,k),n+k-1);
        double temp3 = sqrt(3.0);
        return temp2*temp3*temp1; //sqrt(3)*pow(n+k,(-n-k-2)/2)*power(GetLBasisPol()*GetLPol()*pow(2,k),1-n-k);
    }
    else return rational< big_int>(1, 100);
}

void PolynomialNumber::Calculate_j_B()
{
    if (BasisPol != NULL)
    {
        rational<big_int> temp = rational<big_int>(Arageli::iceil<big_int>(BasisPol->M()), 1);
        rational<big_int> temp1 = temp * (F()/U());

        j_1_B = Arageli::iceil<big_int>(Arageli::log2(temp1.numerator())-Arageli::log2(temp1.denominator()));

        rational<big_int> temp2 = temp * SepL().inverse();

        j_2_B = Arageli::iceil<big_int>(Arageli::log2(temp2.numerator())-Arageli::log2(temp2.denominator()));

        j_3_B = ( j_1_B > j_2_B ? j_2_B : j_1_B );
    }
    else
    {
        j_1_B = 0;
        j_2_B = 0;
        j_3_B = 0;
    }
}

void PolynomialNumber::Calculate_j_C()
{
    if (BasisPol != NULL)
    {
        rational< big_int> temp_1 = F()/U();
        rational< big_int> temp_2 = Arageli::ceil(temp_1);
        rational< big_int> temp_3 = temp_2 * rational< big_int>(225,100) + rational<big_int>(1,1);

        j_1_C = Arageli::iceil<big_int>(Arageli::log2(temp_3.numerator()) - Arageli::log2(temp_3.denominator())) - big_int(1);

        rational< big_int> temp_4 = ceil<big_int>(SepL().inverse()) * rational< big_int>(225,100) + rational<big_int>(1,1);
        big_int temp_5 = Arageli::log2(temp_4.numerator()) - Arageli::log2(temp_4.denominator());

        rational<big_int> temp_6(8, 5);
        rational<big_int> temp_7 = Arageli::log2(temp_6.numerator()) - Arageli::log2(temp_6.denominator());
        rational<big_int> temp_8 = temp_5/temp_7;

        j_2_C = Arageli::iceil<big_int>(temp_8) - big_int(1);

        j_3_C = ( j_1_C > j_2_C ? j_2_C : j_1_C );
    }
    else
    {
        j_1_C = 0;
        j_2_C = 0;
        j_3_C = 0;
    }
}

int PolynomialNumber::sign() const
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    return s(*this);
}


int PolynomialNumber::signPOL_abs(const PolynomialNumber &POL)
{
    sign_binarypolnum_abs_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_u(const PolynomialNumber &POL)
{
    sign_binarypolnum_u_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_sepl(const PolynomialNumber &POL)
{
    sign_binarypolnum_sepl_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_u_sepl(const PolynomialNumber &POL)
{
    sign_binarypolnum_u_sepl_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_j_1_b(const PolynomialNumber &POL)
{
    sign_binarypolnum_j_1_b_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_j_2_b(const PolynomialNumber &POL)
{
    sign_binarypolnum_j_2_b_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signPOL_j_3_b(const PolynomialNumber &POL)
{
    sign_binarypolnum_j_3_b_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_abs(const PolynomialNumber &POL)
{
    sign_chainfraction_abs_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_u(const PolynomialNumber &POL)
{
    sign_chainfraction_u_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_sepl(const PolynomialNumber &POL)
{
    sign_chainfraction_sepl_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_u_sepl(const PolynomialNumber &POL)
{
    sign_chainfraction_u_sepl_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_j_1_c(const PolynomialNumber &POL)
{
    sign_chainfraction_j_1_c_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_j_2_c(const PolynomialNumber &POL)
{
    sign_chainfraction_j_2_c_stop_alg< PolynomialNumber> s;
    return s(POL);
}

int PolynomialNumber::signCFMPOL_j_3_c(const PolynomialNumber &POL)
{
    sign_chainfraction_j_1_c_stop_alg< PolynomialNumber> s;
    return s(POL);
}



//\ Method's //------------------------------------------------------------------------------

/// Shevchenko - Gruzdev signum method's
int PolynomialNumber::signums(const s_p_rat f, const s_p_int g)
{
    s_p_int gmain = s_p_int(g);
    s_p_int gl;
    s_p_rat fmain = s_p_rat(f);
    s_p_rat fl;
    int m = gmain.degree();

    Arageli::matrix<big_int> Frob(m, 0, fromval);
    Arageli::matrix<big_int> D(m,m,0);

    int k = fmain.degree();
    big_int rdet = 0;

    int j = 1;

    while(m%2 == 0 && m > 0)
    {
        Frobeniuss(m, gmain, Frob);
        D = fmain.subs(Frob);
        rdet = det(D);

        if (rdet > 0)
        {
            for(monoms_r mi = fmain.monoms_begin(), mj = fmain.monoms_end(); mi != mj; ++mi)
                if (is_even(mi->degree()))
                    fl.insert(fl.monoms_end(), s_p_rat::monom(mi->coef(), mi->degree()/2));
        }
        else
        {
            for(monoms_r mi = fmain.monoms_begin(), mj = fmain.monoms_end(); mi != mj; ++mi)
                if (is_odd(mi->degree()))
                    fl.insert(fl.monoms_end(), s_p_rat::monom(mi->coef(), (mi->degree()-1)/2));
        }

        k = fl.degree();

        for(monoms mi = gmain.monoms_begin(), mj = gmain.monoms_end(); mi != mj; ++mi)
            gl.insert(gl.monoms_end(), s_p_int::monom(mi->coef(), mi->degree()/2));

        m = gl.degree();

        Frob.resize(m,m);
        Resz(m, Frob);

        D.resize(m,m);
        Resz(m, D);

        gmain = gl;     fmain = fl;

        fl = "0";   gl = "0";

        j++;
    }

    Frobeniuss(m, gmain, Frob);
    D = fmain.subs(Frob);
    rdet = det(D);

    return j * Arageli::sign(rdet);
}

PolynomialNumber& PolynomialNumber::operator --()
{
    sparse_polynom<rational<big_int>> pol = this->Pol();
    pol--;
    this->Pol(pol);
    return *this;
}

PolynomialNumber& PolynomialNumber::operator ++()
{
    sparse_polynom<rational<big_int>> pol = this->Pol();
    pol++;
    this->Pol(pol);
    return *this;
}


} //- end namespace Arageli --------------------------------------------------------