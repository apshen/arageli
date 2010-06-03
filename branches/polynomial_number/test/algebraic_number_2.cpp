// WARNING! This is only compilation test, it isn't integrated in test system.
// Natalia TODO: Please reogranize test as needed and connect to ts.

// algebraic_number_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.hpp"
//#include "arageli/arageli.hpp"

using namespace Arageli;


//- Smith method for matrix< PolynomialNumber> --- does not work!- because std::abs() in smith() does not know PN---
#if 1
//#include <iostream>
//#include <iomanip>


using std::cout;
using std::endl;
using std::size_t;

int pseudo_test (int argc, char* argv[])
{
    basis_field baspol("x^3-2");

    clock_t start, finish;

    PolynomialNumber f11(baspol);
    PolynomialNumber f12(baspol);
    PolynomialNumber f13(baspol);
    PolynomialNumber f21(baspol);
    PolynomialNumber f22(baspol);
    PolynomialNumber f23(baspol);
    PolynomialNumber f31(baspol);
    PolynomialNumber f32(baspol);
    PolynomialNumber f33(baspol);

    f11.Pol("x+1");     f12.Pol("3");   f13.Pol("2");
    f21.Pol("5");       f22.Pol("1");   f23.Pol("2");
    f31.Pol("3");       f32.Pol("8");   f33.Pol("1");

    matrix< PolynomialNumber> S(3, f11, fromval);
/*    S = "((x+1,3,2), (5,1,2), (3,8,1))";
*/

    S(0,1)= f12;  S(0,2)= f13;  S(1,0)= f21;  S(1,1)= f22;  S(1,2)= f23;  S(2,0)= f31;  S(2,1)= f32;  S(2,2)= f33;

    matrix< PolynomialNumber> B(3, f11, fromval);
    matrix< PolynomialNumber> P(3, f12, fromval);
    matrix< PolynomialNumber> Q(3, f12, fromval);

    size_t rk;
    PolynomialNumber d(baspol);

//- output the date -----------------------------------------------------------------------------------------------------
    cout << "\nPolynomialNumber matrix S = \n" << "basis field: " << S(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ ) cout << S(i, j).Pol() << "   ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix B = \n" << "basis field: " << B(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ ) cout << B(i, j).Pol() << "  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix P = \n" << "basis field: " << P(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ ) cout << P(i, j).Pol() << "  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix Q = \n" << "basis field: " << Q(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ ) cout << Q(i, j).Pol() << "  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber d = " << d.BasisPol->BasisPol() << "  d = " << d.Pol() << endl;
//- end outputing ------------------------------------------------------------------------------------------------------

    smith(S, B, P, Q, rk, d);

//- print of results ---------------------------------------------------------------------------------------------------
    cout << "\nPolynomialNumber d.basis_field = " << d.BasisPol->BasisPol() << endl;
    cout << "\nPolynomialNumber d.Pol = " << d.Pol() << endl;

    cout << "\nPolynomialNumber matrix S = \n";
    cout << "basis field: " << S(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ )
            cout << S(i, j).Pol() << "  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix B = \n";
    cout << "basis field: " << B(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ )
            cout << B(i, j).Pol() << "  |  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix P = \n";
    cout << "basis field: " << P(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ )
            cout << P(i, j).Pol() << "  ";
        cout << endl;
    }

    cout << "\nPolynomialNumber matrix Q = \n";
    cout << "basis field: " << Q(0,0).BasisPol->BasisPol() << endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++ )
            cout << Q(i, j).Pol() << "  ";
        cout << endl;
    }
//- end print of results -----------------------------------------------------------------------------------------------


    return 0;
}
#endif