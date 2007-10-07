/*****************************************************************************

    test/test1.cpp

    This file is a part of Arageli library.

    Copyright (C) 2005--2007 Sergey S. Lyalin
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

#include "stdafx.hpp"

#include "test1.hpp"

#define _ARAGELI_TEST_VECTOR_SWAP_BUG

namespace std
{

template <typename T1, typename T2>
inline ostream& operator<< (ostream& out, const pair<T1, T2>& x)
{ return out << "key = " << x.first << ", value = " << x.second; }

}

namespace Arageli
{


/*template <typename P1, typename P2>
void sparse_polynom_tester (std::ostream& report)
{
    // Сейчас, в основном, только тесты на компилируемость.

    binsym_operator_tester<P1, P2>(report);
    binsym_operator_tester<typename P1::monom, typename P2::monom>(report);
    binsym_operator_tester<P1, typename P2::monom>(report);
    binsym_operator_tester<typename P1::coef_type, typename P2::coef_type>(report);
    binsym_operator_tester<P1, typename P2::coef_type>(report);
    binsym_operator_tester<typename P1::monom, typename P2::coef_type>(report);
    binassign_operator_tester<P1, P2>(report);
    binassign_operator_tester<typename P1::monom, typename P2::monom>(report);
    binassign_operator_tester<P1, typename P2::monom>(report);
    binassign_operator_tester<typename P1::coef_type, typename P2::coef_type>(report);
    binassign_operator_tester<P1, typename P2::coef_type>(report);
    binassign_operator_tester<typename P1::monom, typename P2::coef_type>(report);

    binsym_operator_tester<typename P1::monom, P2>(report);
    binsym_operator_tester<typename P1::coef_type, P2>(report);
    binsym_operator_tester<typename P1::coef_type, typename P2::monom>(report);

    binsym_operator_tester<P2, P1>(report);
    binsym_operator_tester<typename P2::monom, typename P1::monom>(report);
    binsym_operator_tester<P2, typename P1::monom>(report);
    binsym_operator_tester<typename P2::coef_type, typename P1::coef_type>(report);
    binsym_operator_tester<P2, typename P1::coef_type>(report);
    binsym_operator_tester<typename P2::monom, typename P1::coef_type>(report);
    binassign_operator_tester<P2, P1>(report);
    binassign_operator_tester<typename P2::monom, typename P1::monom>(report);
    binassign_operator_tester<P2, typename P1::monom>(report);
    binassign_operator_tester<typename P2::coef_type, typename P1::coef_type>(report);
    binassign_operator_tester<P2, typename P1::coef_type>(report);
    binassign_operator_tester<typename P2::monom, typename P1::coef_type>(report);

    binsym_operator_tester<typename P2::monom, P1>(report);
    binsym_operator_tester<typename P2::coef_type, P1>(report);
    binsym_operator_tester<typename P2::coef_type, typename P1::monom>(report);

}


template <typename T1, typename T2>
void binsym_operator_tester (std::ostream& report)
{
    T1 t1;
    T2 t2;

    t1 + t2;
    t1 - t2;
    t1 * t2;
    t1 / t2;
    //t1 % t2;
}


template <typename T1, typename T2>
void binassign_operator_tester (std::ostream& report)
{
    T1 t1;
    T2 t2;

    T1 t12 = t2;

    T1 t13(t2);

    t1 = t2;

    t1 += t2;
    t1 -= t2;
    t1 *= t2;
    t1 /= t2;
    //t1 %= t2;
}
*/


/*
template <typename T> void standard_sparse_polynom_tester_1 (std::ostream& report)
{
    sparse_polynom_tester<sparse_polynom<char>, T>(report);
    sparse_polynom_tester<sparse_polynom<signed char>, T>(report);
    sparse_polynom_tester<sparse_polynom<unsigned char>, T>(report);
    sparse_polynom_tester<sparse_polynom<short>, T>(report);
    sparse_polynom_tester<sparse_polynom<int>, T>(report);
    sparse_polynom_tester<sparse_polynom<long>, T>(report);
    sparse_polynom_tester<sparse_polynom<unsigned short>, T>(report);
    sparse_polynom_tester<sparse_polynom<unsigned int>, T>(report);
    sparse_polynom_tester<sparse_polynom<unsigned long>, T>(report);
    sparse_polynom_tester<sparse_polynom<double>, T>(report);
    sparse_polynom_tester<sparse_polynom<float>, T>(report);
    sparse_polynom_tester<sparse_polynom<long double>, T>(report);
}

void standard_sparse_polynom_tester_2 (std::ostream& report)
{
    //standard_sparse_polynom_tester_1<sparse_polynom<char> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<signed char> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<unsigned char> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<short> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<int> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<long> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<unsigned short> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<unsigned int> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<unsigned long> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<double> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<float> >(report);
    //standard_sparse_polynom_tester_1<sparse_polynom<long double> >(report);
}
*/


void test1_1 ()
{
    big_int a = 1;

    for(int i = 1; i < 1000; ++i)
        a *= i;

    std::cout << "\na = " << a << std::endl;

    rational<big_int> b = a;

    for(int i = 1; i < 1000; ++i)
        b /= i;

    std::cout << "\nb = " << b << std::endl;
}


void test1_2 ()
{
    big_int a;
    a = big_int(1) << big_int(10000);
    a <<= 10000;

    std::cout << "\na = " << a << std::endl;

    big_int b = power(big_int(2), 10000);
    b *= power(big_int(2), 10000);

    std::cout
        << "\nb = " << b << std::endl
        << "\nresults are equal: " << std::boolalpha << (a == b) << std::endl;
}


void test1_3 ()
{
    rational<big_int> b = power(rational<big_int>(2), 10000);
    b *= power(rational<big_int>(2), 10000);

    std::cout
        << "\nb = " << b << std::endl;
}


void test1_4 ()
{
    big_int a(1), b = "2";
    rational<big_int> c = 2, d(2, 7);
    int e = 10, f = 11;

    rational<int> g = a*c*f + 17.0 - d*2;

    std::cout
        << '\n' << a*c
        << '\n' << a*c*f
        << '\n' << a*c*f + 17
        << '\n' << d*2
        << '\n' << -d*2
        << '\n' << a*c*f + 17 - d*2;

    std::cout
        << "\ng = " << rational<big_int>(g) << std::endl
        << "\nresult is correct: " << std::boolalpha << (rational<big_int>(g) == rational<int>(39*7 - 4, 7)) << std::endl;


    std::cout << std::endl;
}


void test1_5 ()
{
    typedef rational<big_int> T;

    matrix<T> a(4), b(4), q(4);
    vector<size_t> basis(4);
    T det;
    std::generate(a.begin(), a.end(), std::rand);
    output_aligned(std::cout << "a == \n", a, "|| ", " ||");

    rref(a, b, q, basis, det);

    output_aligned(std::cout << "\n\na == \n", a, "|| ", " ||");
    output_aligned(std::cout << "\n\nb == \n", b, "|| ", " ||");
    output_aligned(std::cout << "\n\nq == \n", q, "|| ", " ||");
    output_aligned(std::cout << "\n\nbasis == \n", basis, "|| ", " ||");

    std::cout << "\n\ndet == " << det;

    matrix<T> aq = a*q;
    output_aligned(std::cout << "\n\na*q == \n", aq, "|| ", " ||");
    output_aligned(std::cout << "\n\n(a*q)/det == \n", aq/=det, "|| ", " ||");

    std::cout << std::endl;
}


void test1_6 ()
{
    typedef matrix<big_int> MB;
    typedef sparse_polynom<MB> PMB;

    PMB pmb1 = PMB(MB(2, 2, fromval), 2) + PMB(MB(2, 5, fromval));
    const char* nums[] = {"11111111111111111", "22222222", "-3333333333333333333", "4"};
    MB mb1(2);
    std::copy(nums, nums + 4, mb1.begin());

    std::cout
        << "mb1 == " << mb1
        << "\npmb1 == " << pmb1;

    pmb1 += mb1;

    std::cout << "\n\npmb1 += mb1 == " << pmb1;


    std::cout << std::endl;
}


void test1_7 ()
{
    big_int a(std::numeric_limits<int>::min());
    big_int b('1');
    big_int c(123345lu);
    big_int d(true);
    //big_int e(-9223372036854775807);
    big_int f1(1.1);
    big_int f2(-1.1);
    big_int f3(0.1);
    big_int f4(-0.1);
    big_int f5(0.0);
    big_int f6(+1.0);
    big_int f7(-1.0);
    big_int f8(1.0e9);
    big_int f9(1.0e19);
    big_int f10(1.23456789123456789123456789e9);
    big_int f11(1.23456789123456789123456789e19);
    big_int f12(std::numeric_limits<float>::max());
    big_int f13(std::numeric_limits<double>::max());
    big_int f14(std::numeric_limits<long double>::max());

    std::cout << a << ' ' << b << ' ' << c << ' ' << d << ' ' << "-9223372036854775807" << ' ' << "-9223372036854775808";

    std::cout
        << "\n\n"
        << f1 << '\n'
        << f2 << '\n'
        << f3 << '\n'
        << f4 << '\n'
        << f5 << '\n'
        << f6 << '\n'
        << f7 << '\n'
        << f8 << '\n'
        << f9 << '\n'
        << f10 << '\n'
        << f11 << '\n'
        << f12 << '\n'
        << f13 << '\n'
        << f14 << '\n';

    std::cout << '\n' << std::setprecision(300) << std::numeric_limits<double>::max();

    std::cout << std::endl;
}


void test1_8 ()
{
    matrix<rational<int> > mri;
    rational<int> ri;

    //mri = ri*mri;
    mri *= ri;

    sparse_polynom<matrix<int> > pmi;
    matrix<int> mi;

    pmi += mi;


    std::cout << std::endl;
}

void test1_9 ()
{
    big_int a;
    if(a != 0)std::cout << "\ntrue";
    else std::cout << "\nfalse";

    std::cout << std::endl;
}


void test1_10 ()
{
    big_int a(125l);
    int b = a;
    signed char c = a;
    std::cout << a << ' ' << b << ' ' << c;

    a = std::numeric_limits<int>::min();
    b = a;

    std::cout << "\n" << b << ' ' << a;

    a = std::numeric_limits<int>::max();
    b = a;

    std::cout << "\n" << b << ' ' << a;

    std::cout << std::endl;
}


void test1_11 ()
{
    big_int a("99834729384723984273984273498274972394872394729587359834795734958739573957398520934823094820492834098759793473984809284092840923840247884092412398753968701298312");
    float b = a;
    double c = a;
    long double d = a;

    std::cout << a << '\n' << b << '\n' << c << '\n' << d;

    std::cout << std::endl;
}


void test1_12 ()
{
    std::cout << big_int("2") * 12.3 - 2 + 7.0l / big_int(234);

    std::cout << std::endl;
}

void test1_13 ()
{
    big_int a; rational<big_int> b;
    std::cin >> a;
    std::cout << a;
    std::cin >> b;
    std::cout << b;
    std::cin.get();

    std::cout << std::endl;
}


void test1_14 ()
{
    sparse_polynom<int> a;
    std::cout << "Polynom Input Test. Please enter 13 to the next prompt if you want to exit.\n\n";

    do
    {
        std::cout << "Input sparse_polynom: ";
        std::cin >> a;

        std::cout << "State of std::cin is";
        if(std::cin.good())std::cout << " good";
        if(std::cin.eof())std::cout << " eof";
        if(std::cin.bad())std::cout << " bad";
        if(std::cin.fail())std::cout << " fail";

        std::cin.clear();
        std::cin.ignore(1);

        std::cout << "\nPolynom that has been read: " << a << "\n\n";
    }while(a != sparse_polynom<int>(13));

    std::cout << std::endl;
}


void test1_15 ()
{
    typedef matrix<sparse_polynom<matrix<rational<big_int> > > > Mpmrbi;
    typedef vector<int> Vi;

    Vi vi;
    std::cin >> vi;
    std::cout << vi;

    //Mpmrbi a(2, Mpmrbi::value_type("x"));

    //std::cout << a;

    std::cout << std::endl;
}


void test1_16 ()
{
    typedef matrix<sparse_polynom<matrix<rational<big_int> > > > Mpmrbi;

    Mpmrbi a(2, 2, Mpmrbi::value_type("x"));

    std::cout << a;

    std::cout << std::endl;
}


void test1_17 ()
{
    sparse_polynom<big_int> a("x^2-4*x+4"), b("x-2");
    std::cout << "\n\ngcd(" << a << ", " << b << ") == " << gcd(a, b) << "\n\n";

    matrix<sparse_polynom<rational<big_int> > > mprb(2);

    mprb(0, 0) = sparse_polynom<rational<big_int> >(a)/7;
    mprb(1, 0) = b;
    mprb(0, 1) = a*b;
    mprb(1, 1) = a*a;

    output_aligned(std::cout, mprb, "|| ", " ||", "  ");

    std::cout << std::endl;
}


void test1_18 ()
{
    vector<double, false> a(10u, 2.0);
    vector<double, false> b = std::sin(a);
    b = std::sin(a);

    output_aligned(std::cout, a);
    std::cout << std::endl;
}


void test1_19 ()
{
    vector<double, false> a(5);
    for(size_t i = 0; i < a.size(); ++i)
        a[i] = std::rand();

    matrix<double, false> res(100, 5, fromsize);

    for(size_t i = 0; i < res.nrows(); ++i)
    {
        a = std::sin(a);
        for(size_t j = 0; j < a.size(); ++j)
            res(i, j) = a[j];
    }

    output_aligned(std::cout, res);

    std::cout << std::endl;
}


void test1_m ()
{
    typedef sparse_polynom<matrix<rational<big_int> > > P;
    P p1 = P(P::coef_type(2, rational<big_int>(3, 7), fromval), 2) + P(P::coef_type(2, rational<big_int>(1, 2), fromval), 5);
    std::cout << p1 << std::endl;
    p1.leading_coef()(1, 1) -= 5;
    std::cout
        << p1 << "\n\n" << p1*p1 << "\n\n"
        << p1.subs(matrix<sparse_polynom<rational<big_int> > >(2, sparse_polynom<rational<big_int> >(big_int(1), 1), diag)) << '\n'
        << (p1*p1).subs(matrix<sparse_polynom<rational<big_int> > >(2, sparse_polynom<rational<big_int> >(big_int(1), 1), diag)) << '\n';

    std::cout << "\n";

    vector<sparse_polynom<big_int> > vpb1(3u, big_int("123456789123456789"));
    vector<sparse_polynom<int> > vpi1(3u, 123);

    std::cout << '\n' << vpb1 << '\n' << vpi1 << '\n';

    vpb1 += vpi1;

    std::cout << vpb1;

    vpb1 -= vpi1;

    std::cout << vpb1;

    std::cout << std::endl;
}


void test1_20 ()
{
    std::cout << typeid(_Internal::digit).name() << "\n";

    sparse_polynom<rational<big_int> > a, b;
    a = rational<big_int>(24238428);
    a *= a; a *= a;
    b = rational<big_int>(2042343243);
    //sparse_polynom<big_int> a, b;
    //a = big_int(24238428);
    //a *= a; a *= a;
    //b = big_int(23984923749980293);
    //sparse_polynom<long long> a, b;
    //a = long long(24238428);
    //a *= a; a *= a;
    //b = long long(23984923749980293);
    a % b;
    std::cout
        << a << "\n" << b << "\n"
        << a + b << "\n"
        << a * b << "\n"
        << a / b << "\n"
        << a - b << "\n"
        << a % b << "\n";

    std::cout << gcd(a, b);

    std::cout << std::endl;
}


void test1_21 ()
{
    typedef big_int E;
    typedef vector<E> V;
    V v = "(345345345, 204982394723948729847892759573987325,  23432423424234235625)";
    std::cout << v << "\n" << lcm(v);

    std::cout << std::endl;
}


void test1_22 ()
{
    sparse_polynom<rational<big_int> > a;
    while(a != sparse_polynom<rational<big_int> >(13))
    {
        std::cin.clear();
        std::cin >> a;
        if(!std::cin)std::cout << "std::cin is not good\n";
        else std::cout << a << '\n';
    }

    std::cout << std::endl;
}


void test1_23 ()
{
    sparse_polynom<sparse_polynom<rational<> > > v /*= "(x^3 - 12*x)*x^7 + (7)*x^2 - (x)*x"*/;
    //sparse_polynom<rational<> > v /*= "(x^3 - 12*x)*x^7 + (7)*x^2 - (x)*x"*/;
    //input_list(std::cin, v, "", ".", ";", "", "", ",");
    std::cin >> v;
    std::cout << v;

    std::cout << std::endl;
}


void test1_23_5 ()
{
    vector<big_int, true> a;
    vector<big_int, false> b;
    std::swap(a, b);
#ifndef _ARAGELI_TEST_VECTOR_SWAP_BUG
    swap(a, b);
#endif
}


void test1_24 ()
{
    refcntr_proxy<vector<big_int, true> > a;
    refcntr_proxy<vector<big_int, true> , false> b;

    refcntr_proxy<vector<big_int, false> > c;
    refcntr_proxy<vector<big_int, false> , false> d;

    std::swap(a, b);
    std::swap(b, a);
    std::swap(a, a);
    std::swap(b, b);
#ifndef _ARAGELI_TEST_VECTOR_SWAP_BUG
    std::swap(a, c);
    std::swap(b, c);
    std::swap(a, d);
    std::swap(b, d);
    std::swap(c, b);
    std::swap(c, a);
    std::swap(d, a);
    std::swap(d, b);
#endif

    std::cout << std::endl;
}


void test1_25 ()
{
    sparse_polynom<big_int> a("234*x^123-923874*x+23432-x^1");
    sparse_polynom<int> b("1111111*x+22222222*x^2-333333*x^3"), c;

    std::cout << "a = " << a << "\nb = " << b << "\nswapping b <--> a\n";

    std::swap(a, b);
    std::cout << "a = " << a << "\nb = " << b;

    std::cout << "\na = " << a << "\nb = " << b << "\nswapping b <--> c\n";

    std::swap(b, c);
    std::cout << "a = " << a << "\nb = " << b << "\nc = " << c;


    std::cout << std::endl;
}


void test1_26 ()
{
    vector<int> a = "(1, 2, 3)";
    std::cin >> a;
    std::cout << a;

    std::cout << std::endl;
}


void test1_27 ()
{
    rational<int> i;
    rational<big_int> b;
    std::cout << typeid(i + b).name();

    std::cout << std::endl;
}


void test1_28 ()
{
    matrix<sparse_polynom<rational<> > > a, b, p, q;
    //random_int_matrix(a, 2, 2);
    a = "((x, x^3+5), (x^2-x-4, x^4-x^3-4*x^2+5*x-5))";
    output_aligned(std::cout << "a = \n", a);
    size_t rank; sparse_polynom<rational<> > det;
    smith(a, b, p, q, rank, det, ctrl::make_smith_slog(std::cout));

    // normalization pass, actually smith does not this
    for(size_t i = 0; i < b.nrows() && !b(i, i).is_null(); ++i)
    {
        p.div_row(i, b(i, i).leading_coef());
        b(i, i) /= b(i, i).leading_coef();
    }

    output_aligned(std::cout << "\nb = \n", b, "|| ", " ||", "  ");
    output_aligned(std::cout << "\np = \n", p, "|| ", " ||", "  ");
    output_aligned(std::cout << "\nq = \n", q, "|| ", " ||", "  ");
    std::cout << "\nrank = " << rank << "\ndet = " << det;

    matrix<sparse_polynom<rational<> > > paq = p*a*q;
    output_aligned(std::cout << "\np*a*q = \n", paq);
    std::cout << "\nb == p*a*q: " << std::boolalpha << (b == paq);

    std::cout << std::endl;
}




void test1_29 ()
{
    typedef std::map<size_t, int> Lens;

    //{
    //    Lens lens;
    //    for(int i = 0; i < 10000; ++i)
    //        ++lens[big_int::random_in_range(1).length()];


    //    std::copy(lens.begin(), lens.end(), std::ostream_iterator<Lens::value_type>(std::cout, "\n"));

    //    std::cout << "\n";
    //}

    //{
    //    Lens lens;
    //    for(int i = 0; i < 10000; ++i)
    //        ++lens[big_int::random_with_length(50).length()];


    //    std::copy(lens.begin(), lens.end(), std::ostream_iterator<Lens::value_type>(std::cout, "\n"));

    //    std::cout << "\n";
    //}

    {
        Lens lens;
        for(int i = 0; i < 1000000; ++i)
            ++lens[big_int::random_with_length_or_less(50).length()];


        std::copy(lens.begin(), lens.end(), std::ostream_iterator<Lens::value_type>(std::cout, "\n"));
    }
}


void test1_30 ()
{
    Timing tmg;
    big_int res_factorial = factorial_even_odd_multiplication(big_int(20000));
    std::cout << "\nFactorial computation time is " << tmg.time() << " seconds.";

    //big_int res2 = factorial_successive_multiplication(big_int(10000));
    //std::cout << res2;
    //if(res_factorial == res2)
    //    std::cout << "OK";
    //else std::cout << "ERROR";
    //for(big_int i = 0; i < 100; ++i)
    //    if(factorial_by_successive_multiplication(i)
    //        != factorial_by_even_odd_multiplication(i))
    //    {
    //        std::cout
    //            << "\nError: "
    //            << factorial_by_successive_multiplication(i)
    //            << ", " << factorial_by_even_odd_multiplication(i);
    //    }
    //tmg.off();
    //std::cout << "\nFactorial computation time is " << tmg.time() << " seconds.";
    //std::cout << "\nEnd, press Enter...";
    //std::cin.get();
    //std::ofstream file("factorial(100000).txt");
    ////tmg.on();
    //file << res_factorial;
    ////tmg.off();
    ////std::cout << "\nFactorial result outputting time is " << tmg.time() << " seconds.";
    //std::cout << "\nEnd, press Enter...";
    //std::cin.get();

    std::cout << std::endl;
}


void test1_31 ()
{
    vector<big_int> a = "(1,2,3,4,5,6)";

    for(int i = 0; i < 20; ++i)
        std::random_shuffle(a), std::cout << a << '\n';

    std::cout << std::endl;
}


template <typename Pair>
struct Less_second_pair
{
    bool operator() (const Pair& a, const Pair&b) const
    { return a.second < b.second; }
};

template <typename Pair>
struct Outputter_pair
{
    void operator() (const Pair& v) const
    { std::cout << '\n' << v.first << " : " << v.second; }
};


struct Omit_const_first_pair
{
    template <typename T1, typename T2>
    std::pair<T1, T2> operator() (const std::pair<const T1, T2>& a) const
    { return std::pair<T1, T2>(a.first, a.second); }

    template <typename T1, typename T2>
    std::pair<T1, T2> operator() (const std::pair<T1, T2>& a) const
    { return a; }
};


void test1_32 ()
{
    typedef vector<big_int> V;
    typedef std::map<V, big_int> MV;

    V a = "(1, 2, 3, 4, 5)";
    std::cout << "Permutations of a = " << a << std::endl;

    const int
        nech = 10,
        np = factorial(a.size()),
        n = np*nech;

    MV mv;

    for(int i = 0; i < n; ++i)
        ++mv[std::random_shuffle(a), a];

    std::cout
        << "\nTheoretical number of permutations is " << factorial(a.size())
        << "\nEmpirical number of permutations is " << mv.size()
        << "\nAverage value for each chunk is " << nech;

    std::cout << "\n\nDistribution by permutations:\n";
    std::for_each(mv.begin(), mv.end(), Outputter_pair<MV::value_type>());
    typedef vector<std::pair<V, big_int> > VMV;
    VMV vmv;
    std::transform(mv.begin(), mv.end(), std::back_inserter(vmv), Omit_const_first_pair());
    std::sort(vmv, Less_second_pair<VMV::value_type>());
    std::cout << "\n\nSorted by values:\n";
    std::for_each(vmv, Outputter_pair<VMV::value_type>());

    std::cout << std::endl;
}


void test1_33 ()
{
    vector<rational<> > a = "(1, 1/2, 1/3, 2/3, 1/4, 2/4, 3/4, 1/5, 2/5, 3/5, 4/5)";
    std::cout << "\na = " << a;
    std::sort(a);
    std::cout << "\na after sorting = " << a;
    vector<rational<> >::iterator i = std::unique(a);
    std::cout << "\na after unique = " << a;
    a.erase(i, a.end());
    std::cout << "\na after erase = " << a;
    std::cout << "\na + 1 = " << a + 1;
    std::cout << "\na + 1.0 = " << a + 1.0;

    std::cout << std::endl;
}


template <typename T>
struct If_in_interval_0_1
{ bool operator() (const T& x) const { return x > 0 && x < 1; } };


void test1_34 ()
{
    const int n = 10;
    typedef vector<rational<int> > Nums;
    Nums nums(n*n);
    int i = 0;

    for(int num = 1; num <= n; ++num)
        for(int den = 1; den <= n; ++den)
            nums[i++] = Nums::value_type(num, den);

    std::sort(nums);

    nums.erase
    (
        std::partition
        (
            nums.begin(), std::unique(nums),
            If_in_interval_0_1<Nums::value_type>()
        ),
        nums.end()
    );

    vector<double> b = nums;

    output_aligned(std::cout, nums);
    output_aligned(std::cout, b);

    std::cout << std::endl;
}


void test1_35 ()
{
    matrix<big_int> a = "((1, 2, 3), (4, 5, 6))";
    matrix<rational<big_int> > b = "((10, 9), (8, 7), (6, 5), (4, 3))";

    output_aligned(std::cout << "\na = \n", a);
    output_aligned(std::cout << "\nb = \n", b);
    std::cout << "\nSwapping...";
    std::swap(a, b);
    output_aligned(std::cout << "\na = \n", a);
    output_aligned(std::cout << "\nb = \n", b);
    std::cout << "\nTraspose...";
    a.transpose(); b.transpose();
    output_aligned(std::cout << "\na = \n", a);
    output_aligned(std::cout << "\nb = \n", b);

    std::cout << std::endl;
}


void test1_36 ()
{
    matrix<rational<> > a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b = a;

    output_aligned(std::cout << "\na = \n", a);

    a.mult_row(0, 2);
    output_aligned(std::cout << "\na = \n", a);

    a.div_row(2, 3);
    output_aligned(std::cout << "\na = \n", a);

    a.add_rows(0, 1);
    output_aligned(std::cout << "\na = \n", a);

    a.sub_rows(2, 1);
    output_aligned(std::cout << "\na = \n", a);

    a.addmult_rows(1, 0, 10);
    output_aligned(std::cout << "\na = \n", a);

    a.addmult_rows(0, 0, 1);
    output_aligned(std::cout << "\na = \n", a);

    std::swap(a, b);

    std::cout << "\n+++++++++++++++++++++++++\n";

    a.transpose();

    a.mult_col(0, 2);
    output_aligned(std::cout << "\na = \n", a);

    a.div_col(2, 3);
    output_aligned(std::cout << "\na = \n", a);

    a.add_cols(0, 1);
    output_aligned(std::cout << "\na = \n", a);

    a.sub_cols(2, 1);
    output_aligned(std::cout << "\na = \n", a);

    a.addmult_cols(1, 0, 10);
    output_aligned(std::cout << "\na = \n", a);

    a.addmult_cols(0, 0, 1);
    output_aligned(std::cout << "\na = \n", a);

    a.transpose();

    output_aligned(std::cout << "\na = \n", a);

    std::cout << std::endl;
}


//template <typename P> struct Dergee_expansion
//{
//    big_int operator() (const P& x) const
//    { return x.is_null() ? -1 : x.degree(); }
//};


template <typename P> struct Less_degree : std::binary_function<P, P, bool>
{
    bool operator() (const P& a, const P& b) const
    {
        if(a.is_null())return false;
        if(b.is_null())return !a.is_null();
        return a.degree() < b.degree();
    }
};


template <typename P>
matrix<P> canonical_matrix (matrix<P> a)
{
    ARAGELI_ASSERT_0(a.is_square());

    if(a.is_empty())return a;

    if(a.size() == 1)
    {
        P& item = a(0, 0);
        if(!item.is_null())a /= item.leading_coef();
        return a;
    }

    ARAGELI_ASSERT_1(a.size() > 1);

    typedef typename matrix<P>::iterator Iter;
    typedef typename matrix<P>::size_type Size;
    typedef typename P::coef_type Coef;

    Size size = a.ncols();

    output_aligned(std::cerr << "\na from canonical_matrix = \n", a, "|| ", " ||");

    for(Size cur = 0; cur < size - 1; ++cur)
    {
        Iter p = std::min_element
        (
            a.begin() + cur*(size + 1), a.end(),
            Less_degree<P>()
        );

        if(p->is_null())return a;

        Size
            ii = (p - a.begin())/size,
            jj = (p - a.begin())%size;

        std::cerr << "\n(ii, jj) = (" << ii << ", " << jj << ")\n";
        output_aligned(std::cerr << "\na from canonical_matrix = \n", a, "|| ", " ||");

        a.swap_rows(cur, ii);
        a.swap_cols(cur, jj);

        output_aligned(std::cerr << "\na from canonical_matrix = \n", a, "|| ", " ||");

        p = a.begin() + cur*(size + 1);    // it is a(cur, cur)

        ARAGELI_ASSERT_1(!p->is_null());

        Coef lc = p->leading_coef();
        a.div_row(cur, lc);
        ARAGELI_ASSERT_1(is_unit(p->leading_coef()));

        output_aligned(std::cerr << "\na from canonical_matrix = \n", a, "|| ", " ||");

        for(Size j = cur + 1; j < size; ++j)
        {
            std::cerr << "\na(cur, j) % *p  =  a(" << cur << ", " << j << ") % " << *p << "  =  " << a(cur, j) << " % " << *p << "  =  " << a(cur, j) % *p << "\n";

            ARAGELI_ASSERT_1(is_null(a(cur, j) % *p));
            a.addmult_cols(j, cur, -(a(cur, j) / *p));
            ARAGELI_ASSERT_1(is_null(a(cur, j)));
        }

        for(Size i = cur + 1; i < size; ++i)
        {
            std::cerr << "\na(i, cur) % *p  =  a(" << i << ", " << cur << ") % " << *p << "  =  " << a(i, cur) << " % " << *p << "  =  " << a(i, cur) % *p << "\n";

            ARAGELI_ASSERT_1(is_null(a(i, cur) % *p));
            a.addmult_rows(i, cur, -(a(i, cur) / *p));
            ARAGELI_ASSERT_1(is_null(a(i, cur)));
        }
    }

    output_aligned(std::cerr << "\na from canonical_matrix = \n", a, "|| ", " ||");

    Iter p = a.end() - 1;
    if(!p->is_null())*p /= p->leading_coef();

    return a;
}


void test1_37 ()
{
    typedef matrix<sparse_polynom<rational<> > > M;

    M a = "((x^3-x, 2*x^2), (x^2+5*x, 3*x))";
    output_aligned(std::cout << "\na = \n", a);
    std::cout << "\ndet(a) = " << det_int(a) << "\n";
    M ca = canonical_matrix(a);
    output_aligned(std::cout << "\nca = \n", ca);
    std::cout << "\ndet(ca) = " << det_int(ca) << "\n";

    M b = "((x, x^3+5), (x^2-x-4, x^4-x^3-4*x^2+5*x-5))";
    output_aligned(std::cout << "\nb = \n", b);
    std::cout << "\ndet(b) = " << det_int(b) << "\n";
    M cb = canonical_matrix(b);
    output_aligned(std::cout << "\ncb = \n", cb);
    std::cout << "\ndet(cb) = " << det_int(cb) << "\n";

    std::cout << std::endl;
}


void test1_38 ()
{
    matrix<sparse_polynom<rational<> > >
        a = "((x, x^2, x^3), (1, 2*x, 3*x^2), (0, 2, 6*x))";

    output_aligned(std::cout << "\na = \n", a);
    std::cout << "\ndet(a) = " << det_int(a);

    std::cout << std::endl;
}


void test1_39 ()
{
    vector<int> a(10, next_prime(100000)), c;
    big_int b;
    factorize(product(a, b), c);
    std::cout << a << "\n" << b << "\n" << c;
    std::cout << "\nIt is valid: " << (a == c);

    std::cout << std::endl;
}


void test1_40 ()
{
    std::cout << is_prime_division(next_prime(big_int("1234567890")));

    std::cout << std::endl;
}


template <typename T, typename B>
vector<T> power_iterations (const matrix<B>& a, const T& c, size_t numiters)
{
    assert(a.is_square());
    assert(!a.is_empty());

    output_aligned(std::cout << "\na = \n", a);

    size_t n = a.nrows();

    vector<T> p(n);

    for(size_t i = 0; i < n; ++i)
    {
        T nnz = std::count
        (
            a.begin() + i*n, a.begin() + (i+1)*n,    // elements in i-th row
            true
        );

        if(nnz)p[i] = T(1)/nnz;
    }

    output_aligned(std::cout << "\np = \n", p);

    matrix<T> e(n, 1, 1);
    matrix<T> vt(n, 1, T(1)/T(n));
    matrix<T> d(n, 1, fromsize);
    vt.transpose();

    for(size_t i = 0; i < n; ++i)
        if(p[i] == 0)d(i, 0) = 1;

    output_aligned(std::cout << "\ne = \n", e);
    output_aligned(std::cout << "\nvt = \n", vt);
    output_aligned(std::cout << "\nd = \n", d);

    matrix<T> pppt(n);

    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            if(a(i, j))pppt(i, j) = p[i];

    output_aligned(std::cout << "\npp = \n", pppt);

    output_aligned(std::cout << "\nd*vt = \n", d*vt);
    output_aligned(std::cout << "\nppp + d*vt = \n", pppt + d*vt);
    output_aligned(std::cout << "\nc*(ppp + d*vt) = \n", matrix<T>(n, c, diag)*(pppt + d*vt));
    output_aligned(std::cout << "\n(1 - c)*e*vt = \n", matrix<T>(n, 1- c, diag)*e*vt);

    pppt = matrix<T>(n, c, diag)*(pppt + d*vt) + matrix<T>(n, 1 - c, diag)*e*vt;

    output_aligned(std::cout << "\nppp = \n", pppt);
    pppt.transpose();
    output_aligned(std::cout << "\npppt = \n", pppt);

    vector<T> pr = vector<T>(vt.size(), vt.begin());
    output_aligned(std::cout << "\ninitial vector = \n", pr);

    for(size_t i = 0; i < numiters; ++i)
    {
        vector<T> y = pppt*pr;
        std::cout << "\ni = " << i
            << ", delta = " << sum(std::abs(pr - y))
            << "\nfdelta = " << double(sum(std::abs(pr - y)));
        pr = y;
    }

    return pr;
}


template <typename B>
void input_sparse_matrix_from_text_file (std::istream& in, matrix<B>& res)
{
    size_t numrows, numcols, numnzs;
    in >> numrows >> numcols >> numnzs;
    if(!in)return;

    res.assign_fromsize(numrows, numcols);

    size_t nzi = 0;
    for(size_t i = 0; i < numrows; ++i)
    {
        size_t lnlen;
        in >> lnlen;
        if(!in)return;

        assert(nzi + lnlen <= numnzs);

        for(size_t j = 0; j < lnlen; ++j)
        {
            size_t t;
            in >> t;
            if(!in)return;
            res(i, t) = 1;
            ++nzi;
        }
    }

    assert(nzi == numnzs);
}


void test1_41 ()
{
    //matrix<bool> a =
    //    "((1, 1, 0),"
    //    " (0, 1, 1),"
    //    " (1, 0, 1))";

    matrix<bool> a;
    std::ifstream file("../../webgraph/webgraph/lz77_test.in.matrix.txt");
    input_sparse_matrix_from_text_file(file, a);

    vector<rational<> > pr = power_iterations(a, rational<>(85, 100), 100);
    output_aligned(std::cout << "\n\nPageRank = \n", pr);
    output_aligned(std::cout << "\n\nPageRank = \n", vector<double>(pr));

    std::cout << std::endl;
}


void test1_42 ()
{
    sparse_polynom<rational<big_int> > p = "2*x^12-x+23-12*x^2+x^7";
    std::cout << "\np = " << p << "\ndiff(p) = " << diff(p);
    std::cout << "\ngcd(p, diff(p)) = " << gcd(p, diff(p));

    std::cout << std::endl;
}


//void test1_43 ()
//{
//    vector<matrix<double> > v =
//        "(((1, 2), (4, 7)), ((1, 2), (2, 0.4)), ((0, 0), (1, 60)), ((1, 2), (3, 4)), ((0, 0), (-1, 0), (0, 2)))";
//
//    for(size_t i = 0; i < v.size(); ++i)
//    {
//        std::cout << "\nv[" << i << "] =\n";
//        output_aligned(std::cout, v[i]);
//        if(v[i].is_square())std::cout << "\ndet = " << det(v[i]);
//        std::cout << "\nrank = " << rank(v[i]) << std::endl;
//    }
//}


void test1_44 ()
{
    matrix<double>
        m1 = "((1, 2, 3), (4234.222, 5, 6))",
        m2 = "((7, 8, 9000.1))";

    output_aligned_ver_pair(std::cout, m1, m2);

    m1.transpose();
    m2.transpose();

    output_aligned_hor_pair(std::cout << '\n', m1, m2);

    matrix<int> m3 = "((111), (2222), (3333333))";

    output_aligned_corner_triplet_br(std::cout << '\n', m3, m2, m1);
}


//void test1_45 ()
//{
//    matrix<rational<> > m = "((1))";
//    std::cout << det(m) << " " << det_int(matrix<double>(m));
//}

void test1_46 ()
{
    matrix<sparse_polynom<int> > a = "((x, 1, 0), (0, x, 1), (0, 0, x))";
    output_aligned(std::cout, a, "|| ", " ||", " ");
}


void test1_47 ()
{
    typedef rational<sparse_polynom<rational<> > > RF;
    RF
        a1 = "x",
        a2 = "0",
        a3 = "1",
        a4 = "1/3",
        a5 = "x-2*x^4",
        a6 = "x/(x^6-12*x+5)",
        a7 = "(x-2)/(x^2-4*x+4)";

    std::cout
        << "a1 = " << a1 << std::endl
        << "a2 = " << a2 << std::endl
        << "a3 = " << a3 << std::endl
        << "a4 = " << a4 << std::endl
        << "a5 = " << a5 << std::endl
        << "a6 = " << a6 << std::endl
        << "a7 = " << a7 << std::endl;

    std::cout << "a1 + a4 = " << a1 + a4 << std::endl;
    std::cout << "a1 * a4 = " << a1 * a4 << std::endl;
    std::cout << "a1 / a4 = " << a1 / a4 << std::endl;
    std::cout << "a5 / a6 = " << a5 / a6 << std::endl;
    std::cout << "a5 * a6 = " << a5 * a6 << std::endl;
    std::cout << "a6 * x^6-12*x+5 = " << a6 * "x^6-12*x+5" << std::endl;
    //std::cout << "a1 + a4" << a1 + a4 << std::endl;
    //std::cout << "a1 + a4" << a1 + a4 << std::endl;
    //std::cout << "a1 + a4" << a1 + a4 << std::endl;

}


void test1_48 ()
{
    Arageli::big_int a, b, c;
    std::cin >> a >> b;
    c = (a+1)*(b+1);
    std::cout << c << "\n";
}


void test1_49 ()
{
    typedef sparse_polynom<rational<> > P;
    std::cout << gcd(P("1"), P("2")) << '\n';
    std::cout << gcd(P("0"), P("2")) << '\n';
    std::cout << gcd(P("1"), P("0")) << '\n';
    std::cout << gcd(P("x"), P("2")) << '\n';
    std::cout << gcd(P("x^2"), P("2")) << '\n';
    std::cout << gcd(P("x^2"), P("x")) << '\n';
    std::cout << gcd(P("x"), P("x^2")) << '\n';
    std::cout << gcd(P("x^12-x+1"), P("x^12-x+1")) << '\n';
    std::cout << gcd(P("x^2+2*x+1"), P("x+1")) << '\n';
    //std::cut << gcd(P("1"), P("2")) << '\n';
    //std::cut << gcd(P("1"), P("2")) << '\n';
}


void test1_50 ()
{
    typedef sparse_polynom<sparse_polynom<sparse_polynom<rational<> > > > P3v;
    P3v p1 = "(x+(x-(x)))";
    std::cout << p1 << "(x+(x-(x)))^5 = " << std::pow(p1, 5) << '\n';
    P3v q, r;
    divide(std::pow(p1, 5) + P3v("((x))"), std::pow(p1, 3), q, r);
    std::cout << "q = " << q << "\nr = " << r;
}


void test1_51 ()
{
    typedef sparse_polynom<sparse_polynom<rational<> > > P2v;
    P2v p = "((x)-x^2)";    // (x - y^2)
    p *= p;    // (x - y^2)^2 = x^2 - 2*x*y^2 + y^4
    std::cout << p << '\n';
    p = p.subs(P2v("((x))"));
    std::cout << p;
}


//void test1_52 ()
//{
//    typedef sparse_polynom<sparse_polynom<rational<> > > P2v;
//    typedef rational<P2v> RP2v;
//    typedef matrix<RP2v> MRP2v;
//
//    MRP2v a =
//        "((  ((x)) , (((0)))  ),"
//        "(  (((0))), (((x)))  ))";
//
//    //a(0, 0) = RP2v("((x))");
//    //a(1, 1) = RP2v("(((x)))");
//
//    output_aligned(std::cout << "\na = \n", a);
//    std::cout << "\ndet(a) = " << det(a);
//    output_aligned(std::cout << "\ninverse(a) = \n", inverse(a));
//}


void test1_53 ()
{
    sparse_polynom<int> a, b;
    for(;;)
        std::cin >> a >> b, std::cout << cmp(a, b) << '\n';
}


//void test1_54 ()
//{
//    typedef
//        sparse_polynom<
//        sparse_polynom<
//        sparse_polynom<
//        sparse_polynom<
//        sparse_polynom<
//        sparse_polynom<
//        rational<> > > > > > > P6v;
//
//    typedef rational<P6v> RP6v;
//    typedef matrix<RP6v> MRP6v;
//    typedef vector<RP6v> VRP6v;
//
//    MRP6v a =
//        "((   ((x))  ,   (((x)))   ),"
//        "(  ((((x)))), (((((x))))) ))";
//
//    VRP6v b = "( ((((((x)))))), (((((((x))))))) )";
//
//    output_aligned(std::cout << "\na = \n", a);
//    output_aligned(std::cout << "\nb = \n", b);
//    MRP6v inv_a = inverse(a);
//    output_aligned(std::cout << "\ninverse of a = \n", inv_a);
//    VRP6v res1 = inv_a * b;
//    output_aligned(std::cout << "\nres1 = \n", res1);
//    VRP6v res2 = solve_linsys(a, b);
//    output_aligned(std::cout << "\nres2 = \n", res2);
//    std::cout << "\nsolution is equal: " << std::boolalpha << (res1 == res2);
//    std::cout << "\nthe first solution is valid: " << (a*res1 == b);
//    std::cout << "\nthe second solution is valid: " << (a*res2 == b);
//}


void test1_55 ()
{
    for(int t = 0; t < 10; ++t)
    {
        double ct = 0;

        for(int s = 0; s < 10; ++s)
        {
            rational<big_int>
                a = rational<>
                    (
                        big_int::random_with_length_or_less((t+1)*32),
                        big_int::random_with_length_or_less((t+1)*32)
                    ),
                b = rational<>
                    (
                        big_int::random_with_length_or_less((t+1)*32),
                        big_int::random_with_length_or_less((t+1)*32)
                    );

            Timing tm;

            for(int i = 0; i < 10000; ++i)
                a + b;

            ct += tm.time();
        }

        std::cout << "\nt = " << t + 1 << "\t time = " << ct;
    }
}


void test1_56 ()
{
    // ************ WARNING! Please uncomment this. ***************

    //typedef std::complex<rational<big_int> > CRBI;
    //
    //CRBI
    //    a("2039875902357", "29287349872394"),
    //    b("1111", "9824987948721984729234");

    //std::cout
    //    << "\na = " << a << "\nb = " << b
    //    << "\na + b = " << a+b << "\na*b = " << a*b
    //    << "\na/b = " << a/b;

    //sparse_polynom<CRBI> aa = "x^3+(0,1)", bb = "((0,1)*x^5+2*x)";
    //aa.leading_coef() = a;
    //bb -= sparse_polynom<CRBI>::monom(b, 4);

    //std::cout << "\naa = " << aa << "\nbb = " << bb;
    //std::cout << "\naa + bb = " << aa + bb << "\naa*bb = " << aa*bb;

}


void test1_57 ()
{
    sparse_polynom<rational<> > a = "1/7+12*x^23-23/9873*x+x^2";
    std::cout << is_primitive(a);
}


void test1_58 ()
{
    sparse_polynom<rational<> > a = "1/7+12*x^23-23/9873*x+x^2";
    std::cout << is_primitive(a);
}


void test1_59 ()
{
    big_int a = "35278925";
    //big_int a = 10;
    sparse_polynom<rational<> > p;
    for(size_t i = 0; i < a.length(); ++i)
        if(a.bit(i))p += sparse_polynom<rational<> >::monom(1, i);

    //std::cout << "\np = " << p << "\ndiff p = " << diff(sparse_polynom<bool>(p)) <<
    //    "\ngcd = " << gcd(p, sparse_polynom<rational<> >(diff(sparse_polynom<bool>(p))));
}


void test1_60 ()
{
    unsigned long long a = next_prime(487198273);
    unsigned long long b = next_prime(722347161);
    unsigned long long ab = a*b;
    std::cout << "time";
    std::cin.get();
    std::cout << ab << " = " << factorize(ab);
}


void test1_61 ()
{
    // правильные остатки
    std::cout
        << '\n' << (-1) % 4
        << '\n' << (-2) % 4
        << '\n' << (-5) % 4
        << '\n' << big_int(-1) % big_int(4)
        << '\n' << big_int(-2) % big_int(4)
        << '\n' << big_int(-5) % big_int(4);
}


void test1_62 ()
{
    sparse_polynom<int> a = 5;
    //sparse_polynom_input_var(std::cin, a, "x", "");
    std::cout << a << std::endl;
}


void test1_63 ()
{
    std::cout
        << "\nJ(1209478239875929287349875936, 52349124378234343) = "
            << jacobi(big_int("1209478239875929287349875936"), big_int("52349124378234343"))
        << "\nJ(136, 53) = " << jacobi(136, 53)
        << "\nJ(1, 99) = " << jacobi(1, 99);
}

void test1_64 ()
{
    typedef big_int T;

    for(T n = 2; n < 100; n = next_prime(n))
        for(T a = 1; a < n; ++a)
        {
            T r = inverse_mod(a, n);
            std::cout
                << "\n" << a << "^(-1) (mod " << n << ") = "
                << r;

            if(r*a % n != 1)
                std::cerr << "\n!!! ERROR !!!: a = " << a << ", n = " << n;
        }
}


void test1_65 ()
{
    typedef big_int T;
    T end = big_int("100000000000000000000001000");

    for(T i = "100000000000000000000000000"; i < end; ++i)
    {
        std::cout << "\nsqrt(" << i << ") = " << std::sqrt(i);
    }
}

enum simplex_method_result { SMR_OK, SMR_INCOMPATIBLE, SMR_INFINITE };


template <typename T, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& my_output_latex
(
    std::basic_ostream<Ch, ChT>& out, const T& x,
    bool mathmode = false
)
{ return out << x; }


template <typename T, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& my_output_latex
(
    std::basic_ostream<Ch, ChT>& out,
    const rational<T>& x,
    bool mathmode = false,
    size_t msaf = 5
)
{
    typedef std::basic_string<Ch, ChT> S;

    std::basic_ostringstream<Ch, ChT> buf;
    my_output_latex(buf, x.numerator(), false);
    S snumer = buf.str();

    if(x.is_integer())
    {
        if(mathmode)
            out << "$" << snumer << "$";
        else
            out << snumer;

        return out;
    }

    buf.str(S());
    my_output_latex(buf, x.denominator(), false);
    S sdenom = buf.str();

    if(mathmode)out << "$";
    if(snumer.length() + sdenom.length() <= msaf)
        out << "{" << snumer << "}/{" << sdenom << "}";
    else
        out << "\\frac{" << snumer << "}{" << sdenom << '}';
    if(mathmode)out << "$";
    return out;
}


template <typename T, bool REFCNT, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& my_output_latex
(
    std::basic_ostream<Ch, ChT>& out,
    const matrix<T, REFCNT>& x,
    bool mathmode = false,
    size_t msaf = 5
)
{
    typedef std::basic_string<Ch, ChT> S;

    //std::basic_ostringstream<Ch, ChT> buf;
    //buf << x.numerator();

    out << "\\begin{pmatrix}";

    for(size_t i = 0; i < x.nrows(); ++i)
    {
        if(i != 0)out << " \\\\ ";
        for(size_t j = 0; j < x.ncols(); ++j)
        {
            if(j != 0)out << " & ";
            my_output_latex(out, x(i, j), false);
        }
    }

    out << "\\end{pmatrix}";
    return out;
}


template <typename T, bool REFCNT, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& output_simplex_table_latex
(
    std::basic_ostream<Ch, ChT>& out,
    const matrix<T, REFCNT>& x,
    bool mathmode = false,
    size_t pivot_row = std::numeric_limits<size_t>::max(),
    size_t pivot_col = std::numeric_limits<size_t>::max()
)
{
    typedef std::basic_string<Ch, ChT> S;

    //std::basic_ostringstream<Ch, ChT> buf;
    //buf << x.numerator();

    out << " \\left( \\begin{tabular}{c|" << S(x.ncols()-1, 'c') << "c}";

    for(size_t i = 0; i < x.nrows(); ++i)
    {
        if(i != 0)out << " \\\\ ";
        if(i == 1)out << "\\hline ";
        for(size_t j = 0; j < x.ncols(); ++j)
        {
            if(j != 0)out << " & ";
            bool pivoting = (i == pivot_row && j == pivot_col);
            if(pivoting)out << "\\fbox{";
            my_output_latex(out, x(i, j), true);
            if(pivoting)out << "}";
        }
    }

    out << "\\end{tabular} \\right) ";
    return out;
}


template <typename T, bool REFCNT, typename Ch, typename ChT>
std::basic_ostream<Ch, ChT>& output_basis_latex
(
    std::basic_ostream<Ch, ChT>& out,
    const vector<T, REFCNT>& x
)
{
    out << "\n$$\n {\\cal B}=\\left<";
    output_list(out, x, "", "");
    out << "\\right>\n$$";
    return out;
}


template
<
    typename Tc, bool REFCNTc,    // vector c
    typename Ta, bool REFCNTa,    // matrix a
    typename Tb, bool REFCNTb,    // vector b
    typename T,
    typename Tx, bool REFCNTx        // vector x
>
simplex_method_result maximize_simplex_method_l
(
    const vector<Tc, REFCNTc>& c,
    const matrix<Ta, REFCNTa>& a,
    const vector<Tb, REFCNTb>& b,
    T& res,
    vector<Tx, REFCNTx>& x    // may be the same as c
)
{
    ARAGELI_ASSERT_0(c.size() > 0);
    ARAGELI_ASSERT_0(c.size() == a.ncols());
    ARAGELI_ASSERT_0(a.nrows() == b.size());

    std::cout << "Решаем ЗЛП, заданную в канонической форме:";
    std::cout <<
        "\n$$"
        "\\max(cx)$$ $$"
        "\\begin{cases}"
        "Ax=b \\\\ x\\ge 0"
        "\\end{cases}"
        "$$\n";
    std::cout << "где\n$$";
    output_list(std::cout << "c=", c);
    std::cout << "\n$$\n$$\n";
    my_output_latex(std::cout << "A=", a);
    std::cout << "\n$$\n$$\n";
    output_list(std::cout << "b=", b);
    std::cout << "^T\n$$\n";
    std::cout << std::flush;

    // Stage 0.  Building of the first allowable plan.

    typedef matrix<T, false> ST;    // a simplex table
    typedef typename ST::size_type size_type;
    size_type n = a.ncols(), m = a.nrows();
    ST st = a;
    st.insert_col(0, b);
    st.insert_cols(1, m, factory<T>::null());
    st.insert_row(0, factory<T>::null());
    std::fill_n(st.begin() + 1, m, factory<T>::unit());

    for(size_type i = 1; i <= m; ++i)
    {
        st(i, i) = factory<T>::unit();
        st.sub_rows(0, i);
    }

    vector<size_type> basis(m);
    for(size_type i = 1; i <= m; ++i)
        basis[i-1] = i;

    std::cout << "Строим начальную допустимую базу и начальную симплекс-таблицу методом"
        " искусственного базиса.\\par";
    //std::cout << "\n$$\n";
    //output_simplex_table_latex(std::cout, st);
    //std::cout << "\n$$\n";
    //output_basis_latex(std::cout, basis);

    for(;;)
    {
        size_type s = std::find_if
        (
            st.begin() + 1, st.begin() + st.ncols(),
            func::is_negative<T>()
        ) - st.begin();

        if(s == st.ncols())break;

        vector<vector<T, false>, false> norm_rows(m);
        vector<size_type, false> positive_rows(m);
        size_type npr = 0;    // number of positive rows

        for(size_type i = 1; i <= m; ++i)
        {
            const T& cur = st(i, s);
            if(is_positive(cur))
            {
                positive_rows[npr] = i;
                vector<T, false>& cv = norm_rows[npr];
                cv.resize(st.ncols());
                for(size_type j = 0; j < st.ncols(); ++j)
                    cv[j] = st(i, j)/cur;
                ++npr;
            }
        }

        ARAGELI_ASSERT_1(npr > 0);

        size_type r = positive_rows
        [
            min_element (norm_rows.begin(), norm_rows.begin() + npr) -
                norm_rows.begin()
        ];

        std::cout << "\n$$\n";
        output_simplex_table_latex(std::cout, st, false, r, s);
        std::cout << "\n$$\n";
        output_basis_latex(std::cout, basis);

        //std::cout << "Направляющий элемент: $(r = " << r << ", s = " << s << ")$.";

        ARAGELI_ASSERT_1(is_positive(st(r, s)));

        T t = st(r, s);
        st.div_row(r, t);

        for(size_type i = 0; i <= m; ++i)
        {
            if(i == r)continue;
            st.addmult_rows(i, r, -st(i, s));
        }

        basis[r - 1] = s;

    }

    std::cout << "\n$$\n";
    output_simplex_table_latex(std::cout, st, false);
    std::cout << "\n$$\n";
    output_basis_latex(std::cout, basis);

    if(is_negative(st(0, 0)))return SMR_INCOMPATIBLE;
    std::list<size_type> delrows;

    for(size_type r = 1; r <= m; ++r)
    {
        size_type s = basis[r-1];
        ARAGELI_ASSERT_1(s > 0);
        if(s > m)continue;
        ARAGELI_ASSERT_1(is_null(st(r, 0)));

        size_type j = std::find_if
        (
            st.begin() + st.ncols()*r + 1 + m, st.begin() + st.ncols()*(r+1),
            func::sign<T>()
        ) - (st.begin() + st.ncols()*r);

        if(j == st.ncols())
        {
            delrows.push_back(r);
            continue;
        }

        ARAGELI_ASSERT_1(j > m && j < st.ncols());

        T t = st(r, j);
        st.div_row(r, t);
        for(size_type i = 0; i <= m; ++i)
        {
            if(i == r)continue;
            st.addmult_rows(i, r, -st(i, j));
        }

        basis[r-1] = j;
    }

    for
    (
        typename std::list<size_type>::reverse_iterator i = delrows.rbegin();
        i != delrows.rend();
        ++i
    )
    {
        st.erase_row(*i);
        basis.erase(*i-1);
    }

    basis -= m;

    st.erase_cols(1, m);

    m -= delrows.size();


    std::cout << "Симплекс таблица после первой фазы:\n$$\n";
    output_simplex_table_latex(std::cout, st, false);
    std::cout << "\n$$\n";

    for(size_type i = 1; i <= n; ++i)
        st(0, i) = -c[i-1];

    std::cout << "Начинаем вторую фазу симплекс метода.  Допишем вектор $c$:\n$$\n";
    output_simplex_table_latex(std::cout, st, false);
    std::cout << "\n$$\n";

    for(size_type i = 1; i <= m; ++i)
        st.addmult_rows(0, i, -st(0, basis[i - 1]));

    std::cout << "Вторая фаза.\n";
    //output_simplex_table_latex(std::cout, st);
    //std::cout << "\n$$\n";
    //output_basis_latex(std::cout, basis);

    for(;;)
    {
        size_type s = std::find_if
        (
            st.begin() + 1, st.begin() + st.ncols(),
            func::is_negative<T>()
        ) - st.begin();

        if(s == st.ncols())break;

        vector<vector<T, false>, false> norm_rows(m);
        vector<size_type, false> positive_rows(m);
        size_type npr = 0;    // number of positive rows

        for(size_type i = 1; i <= m; ++i)
        {
            const T& cur = st(i, s);
            if(is_positive(cur))
            {
                positive_rows[npr] = i;
                vector<T, false>& cv = norm_rows[npr];
                cv.resize(st.ncols());
                for(size_type j = 0; j < st.ncols(); ++j)
                    cv[j] = st(i, j)/cur;
                ++npr;
            }
        }

        if(npr == 0)return SMR_INFINITE;

        size_type r = positive_rows
        [
            min_element (norm_rows.begin(), norm_rows.begin() + npr) -
                norm_rows.begin()
        ];

        std::cout << "\n$$\n";
        output_simplex_table_latex(std::cout, st, false, r, s);
        std::cout << "\n$$\n";
        output_basis_latex(std::cout, basis);

        ARAGELI_ASSERT_1(is_positive(st(r, s)));

        T t = st(r, s);
        st.div_row(r, t);

        for(size_type i = 0; i <= m; ++i)
        {
            if(i == r)continue;
            st.addmult_rows(i, r, -st(i, s));
        }

        basis[r - 1] = s;

    }

    std::cout << "\n$$\n";
    output_simplex_table_latex(std::cout, st, false);
    std::cout << "\n$$\n";
    output_basis_latex(std::cout, basis);

    res = st(0, 0);

    ARAGELI_ASSERT_0(x.size() >= n);
    x = 0;
    for(size_type i = 0; i < m; ++i)
        x[basis[i]-1] = st(i+1, 0);

    output_list(std::cout << "\n\n{\\bf Ответ:} оптимальный вектор $x=", x);
    std::cout << "^T $, ";
    std::cout << "оптимальное значение:~$" << res << "$.";

    return SMR_OK;
}


void test1_66 ()
{
    typedef rational<> T;
    //vector<T> c = "(0, 2, -1, 1)";
    //matrix<T> a = "((2, 1, 1, 0), (1, 2, 0, 1))";
    //vector<T> b = "(6, 6)";
    vector<T> c;
    matrix<T> a;
    vector<T> b;
    vector<T> x;
    T res;
    const char* status[] = {"OK", "INCOMPATIBLE", "INFINITE"};


    //c = "(-1, -1, 0, -5)";
    //a = "((1, 1, -1, 3), (1, -2, 3, -1), (5, -4, 7, 3))";
    //b = "(1, 1, 5)";
    //x.resize(c.size());
    //
    //std::cout
    //    << "\nStatus = "
    //    << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
    //    << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    //
    //c = "(0, 2, -1, 1)";
    //a = "((2, 1, 1, 0), (1, 2, 0, 1))";
    //b = "(6, 6)";
    //x.resize(c.size());
    //
    //std::cout
    //    << "\nStatus = "
    //    << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
    //    << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    //
    //c = "(2, 1, 0, 0, 0)";
    //a = "((3, 4, 1, 0, 0), (3, 1, 0, 1, 0), (1, 0, 0, 0, 1))";
    //b = "(32, 17, 5)";
    //x.resize(c.size());
    //
    //std::cout
    //    << "\nStatus = "
    //    << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
    //    << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    //
    //c = "(0, 0, 0, 0, 200, 175, -1100, -2)";
    //a =
    //    "((0, 1, 0, 0, -3, -5/4, 7, 1/50),"
    //    " (1, 0, 0, 0, 1/3, 1/6, -1, -1/150),"
    //    " (0, 0, 1, 0, 72/5, -25/4, 175/2, 1/4),"
    //    " (0, 0, 0, 1, 0, 0, 0, 1))";
    //b = "(0, 0, 0, 1)";
    //x.resize(c.size());
    //
    //std::cout
    //    << "\nStatus = "
    //    << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
    //    << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    c = "(4, 4, -1, -4, -4)";
    a =
        "((5, -2, -3, 5, 0),"
        " (-2, -2, 1, 3, -3),"
        " (5, 2, 0, -3, 4))";
    b = "(5, 1, 4)";
    x.resize(c.size());

    std::cout
        << "\nStatus = "
        << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
        << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    c = "(-4, -4, 1, 4, 4)";
    a =
        "((5, -2, -3, 5, 0),"
        " (-2, -2, 1, 3, -3),"
        " (5, 2, 0, -3, 4))";
    b = "(5, 1, 4)";
    x.resize(c.size());

    std::cout
        << "\nStatus = "
        << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
        << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    //c = "(-5, +5, -1, +1, -4, +4, 0, 0, 0, 0, 0)";
    //a =
    //    "((+5, -5, -2, +2, +5, -5, -1, 0,  0,  0,  0),"
    //    " (-2, +2, -2, +2, +2, -2, 0,  -1, 0,  0,  0),"
    //    " (+3, -3, -1, +1, 0,  0,  0,  0,  +1, 0,  0),"
    //    " (-5, +5, -3, +3, +3, -3, 0,  0,  0,  +1, 0),"
    //    " (0,  0,  +3, -3, -4, +4, 0,  0,  0,  0,  +1))";
    //b = "(4, 4, 1, 4, 4)";
    //x.resize(c.size());
    //
    //std::cout
    //    << "\nStatus = "
    //    << status[maximize_simplex_method_l(c, a, b, res, x)] << std::endl
    //    << "++++++++++++++++++++++++++++++++++++++++++++++++++\n";

}


void test1_67 ()
{
    typedef vector<double> V;
    V v;
    std::cin >> v;
    //std::istringstream s;
    //s >> v;
}


void test1_68 ()
{
    typedef monom<double> T;
    T x;

    std::cin >> iomanip::monom_input_list() >> x;
    std::cout << iomanip::monom_output_list() << x << std::endl;
    std::cin >> iomanip::monom_input_var() >> x;
    std::cout << iomanip::monom_output_var(true, "variable", " multiply by ", " in poser of ") << x << '\n';
    std::cout << iomanip::monom_output_aligned() << x;
}


void test1_69 ()
{
    typedef sparse_polynom<double> T;
    T a = "x-2", b = "x^2+5";
    std::cout << (a < b);
}


void test1_70 ()
{
    typedef int T;
    my_output_latex(std::cout << "\n", rational<T>(0));
    my_output_latex(std::cout << "\n", rational<T>(1));
    my_output_latex(std::cout << "\n", rational<T>(1, 13));
    my_output_latex(std::cout << "\n", rational<T>(234, 11));
    my_output_latex(std::cout << "\n", rational<T>(234, 111));
    my_output_latex(std::cout << "\n", rational<T>(23243243, 1141));

    std::cout << "\n";

    matrix<rational<> > m = "((1, 2), (91827/132987687624, 23/2), (1/234, -23223/5))";
    my_output_latex(std::cout << "\n", m);
}


void test1_71 ()
{
    iomanip::vector_input_list vil;
    iomanip::vector_output_aligned voa("|| ", " ||");

    vector<sparse_polynom<rational<> > > v;
    std::cin >> vil >> v;
    std::cout << '\n' << voa << v;
}


void test1_72 ()
{
    matrix<sparse_polynom<int> > mp = "((-x^32-12*x^3+7,12), (x-1,x^2-2*x+1))";
    sparse_polynom<int> p = "-x^32-12*x^3+7";

    std::cout
        << iomanip::matrix_output_aligned() << mp << '\n'
        << iomanip::sparse_polynom_output_aligned() << p;
}


void test1_73 ()
{
    matrix<int> mm = "((1, 2),(3, 4))";
    output_aligned(std::cout << "\n\n", mm);

    {
        matrix<int> m = mm;
        m.insert_col(0, 5);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(1, 5);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(2, 5);
        output_aligned(std::cout << "\n\n", m);
    }

    vector<int> v = "(5, 6)";

    {
        matrix<int> m = mm;
        m.insert_col(0, v);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(1, v);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(2, v);
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<int> m = mm;
        m.insert_col(0, v.begin());
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(1, v.begin());
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.insert_col(2, v.begin());
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<int> m = mm;
        m.erase_col(0);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.erase_col(1);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.erase_cols(0, 2);
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<int> m = mm;
        m.swap_cols(0, 1);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.swap_cols(1, 0);
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<int> m = mm;
        m.mult_col(0, 5);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<int> m = mm;
        m.mult_col(1, 5);
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<rational<int> > m = mm;
        m.div_col(0, 5);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<rational<int> > m = mm;
        m.div_col(1, 5);
        output_aligned(std::cout << "\n\n", m);
    }

    {
        matrix<rational<int> > m = mm;
        m.add_cols(0, 0);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<rational<int> > m = mm;
        m.add_cols(0, 1);
        output_aligned(std::cout << "\n\n", m);
    }
    {
        matrix<rational<int> > m = mm;
        m.add_cols(1, 0);
        output_aligned(std::cout << "\n\n", m);
    }

}

void test1_74 ()
{
    try { big_int(""); } catch(incorrect_string) { std::cout << "All is right!\n"; }
    try { big_int("-"); } catch(incorrect_string) { std::cout << "All is right!\n"; }
    try { big_int("d"); } catch(incorrect_string) { std::cout << "All is right!\n"; }
    try { matrix<big_int>("(234, 234, )"); } catch(incorrect_string) { std::cout << "All is right!\n"; }
    try { sparse_polynom<rational<> >("x^2-1/x"); } catch(incorrect_string) { std::cout << "All is right!\n"; }
}


void test1_75 ()
{
    //unsigned long long a = 92233719531022718101ul;
    //big_int a = "92233719531022718101";
    //std::cout << "\na = " << a << " is prime: " << is_prime(a);
    //std::cout << "\nfactorization is " << partial_factorize_division(a, 1000);
    //a = "9233719531022718101";
    //std::cout << "\na = " << a << " is prime: " << is_prime(a);
    //std::cout << "\nfactorization is " << partial_factorize_division(a, 1000);
    //a = "9223719531022718101";
    //std::cout << "\na = " << a << " is prime: " << is_prime(a);
    //std::cout << "\nfactorization is " << partial_factorize_division(a, 1000);

    //typedef unsigned long long T;

    //T
    //    b1 = next_prime(4100032442u), b2 = next_prime(4102032442u), c = b1*b2;


    //std::cout << b1 << "*" << b2 << " = " << c;
    //
    //vector<T> res;

    //factorize_division(c, res);
    //std::cout << "\nfactorization is " << res;
    //std::cout << "\nfactorization is valid: " << (product(res) == c);

    //a = "9233719531022718101";
    //factorize_division(a, res);
    //std::cout << "\nfactorization is " << res;
    //std::cout << "\nfactorization is valid: " << (product(res) == a);

    //a = "9223719531022718101";
    //factorize_division(a, res);
    //std::cout << "\nfactorization is " << res;
    //std::cout << "\nfactorization is valid: " << (product(res) == a);

    //big_int d = "92233719531022718101";
    //vector<big_int> res2; big_int rest;
    //partial_factorize_division(d, res2, big_int(1000000), rest);
    //std::cout << "\npartial factorization of " << d << " is " << res2 << ", rest is " << rest;

    std::cout << factorize_division(big_int("92233719531022718101"));
    //std::cout << factorize_division(16818466434869209727u);
}


void test1_76 ()
{
    sparse_polynom<int> p;
    p == 0;
}


void test1 (std::ostream& report)
{
#if 1

    test1_1();
    test1_2();
    test1_3();
    test1_4();
    test1_5();
    test1_6();
    test1_7();
    test1_8();
    test1_m();

    test1_9();

    test1_12();
    //test1_13();

    //test1_14();
    //test1_15();
    test1_16();

    test1_17();
    test1_18();
    test1_19();
    test1_20();
    test1_21();
    //test1_22();
    //test1_23();
    test1_24();
    test1_25();
    //test1_26();
    test1_27();
    test1_28();
    //test1_29();
    test1_30();
    test1_31();
    test1_32();
    test1_33();
    test1_34();
    test1_35();
    test1_36();
    //test1_37();
    test1_38();
    test1_39();
    test1_40();
    //test1_41();
    test1_42();
    //test1_43();
    test1_44();

#endif

    //test1_30();
    //test1_43();
    //test1_44();
    //test1_28();
    //test1_46();
    //test1_47();
    //test1_48();
    //test1_49();
    //test1_50();
    //test1_51();
    //test1_52();
    //test1_53();
    //test1_54();
    //test1_55();
    //test1_56();
    //test1_57();
    //test1_58();
    //test1_59();
    //test1_60();
    //test1_61();
    //test1_62();
    //test1_63();
    //test1_64();
    //test1_65();
    //test1_66();
    //test1_68();
    //test1_69();
    //test1_70();
    //test1_71();
    //test1_72();
    //test1_73();
    //test1_74();
    test1_75();

    //P p1 = P(P::coef_type(1), 0) + P(P::coef_type(2, 1), 1) + P(P::coef_type(3, 2), 3) + P::coef_type(7);
    //P p2 = p1*p1 + p1;
    //std::cout << p1 << '\n' << std::pow(p1, P(1)) << std::endl;
    //std::cout
    //    << p1 << '\n' << p2 << '\n'
    //    << gcd(p1, p2) << '\n'
    //    << p1/gcd(p1, p2) << '\n' << p2/gcd(p1, p2) << '\n'
    //    << p1%gcd(p1, p2) << '\n' << p2%gcd(p1, p2) << '\n'
    //    ;
}


// Tests for Arageli::vector  ////////////////////////////////////////////////

void test_vector (std::ostream& report);    // See implementation below.


template <typename CL, typename T, bool REFCNT>
void test_static_declarations (std::ostream& report)
{
    ARAGELI_ASSERT_ALWAYS((equal_types<typename CL::value_type, T>::bvalue));
    ARAGELI_ASSERT_ALWAYS((equal_types<typename CL::const_reference, const T&>::bvalue));
    ARAGELI_ASSERT_ALWAYS(CL::refcounting == REFCNT);
}


template <typename V>
void test_vector_construct (std::ostream&)
{
    V a;
    ARAGELI_ASSERT_ALWAYS(a.size() == 0);

    for(size_t i = 0; i < 100; ++i)
    {
        V b(i);
        ARAGELI_ASSERT_ALWAYS(b.size() == i);
    }

    V c(1), d(2);

    ARAGELI_ASSERT_ALWAYS(c.size() == 1);
    ARAGELI_ASSERT_ALWAYS(d.size() == 2);
    a = c;
    ARAGELI_ASSERT_ALWAYS(a.size() == 1);
    ARAGELI_ASSERT_ALWAYS(a == c);
    c = d;
    ARAGELI_ASSERT_ALWAYS(a.size() == 1);
    ARAGELI_ASSERT_ALWAYS(c.size() == 2);
    ARAGELI_ASSERT_ALWAYS(c == d);

    V e(c), f(a);

    ARAGELI_ASSERT_ALWAYS(e.size() == 2);
    ARAGELI_ASSERT_ALWAYS(f.size() == 1);
    ARAGELI_ASSERT_ALWAYS(e == c);
    ARAGELI_ASSERT_ALWAYS(f == a);
}


template <typename V1, typename V2>
void test_vector_conversion (std::ostream&)
{
    V1 v11(3); V2 v21(2);
    v11 = v21;

    ARAGELI_ASSERT_ALWAYS(v21.size() == 2);
    ARAGELI_ASSERT_ALWAYS(v11.size() == v21.size());

    v21 = v11;

    ARAGELI_ASSERT_ALWAYS(v21.size() == 2);
    ARAGELI_ASSERT_ALWAYS(v21.size() == v11.size());

    V1 v12(v21); V2 v22(v11);
}


template <typename V>
void test_vector_assignment (std::ostream&)
{

}


void test_vector (std::ostream& report)
{
    test_static_declarations
        <vector<int, true>, int, true>
        (report);
    test_static_declarations
        <vector<int, false>, int, false>
        (report);
    test_static_declarations
        <vector<big_int, true>, big_int, true>
        (report);
    test_static_declarations
        <vector<big_int, false>, big_int, false>
        (report);

    test_vector_construct<vector<int, true> >(report);
    test_vector_construct<vector<int, false> >(report);
    test_vector_construct<vector<big_int, true> >(report);
    test_vector_construct<vector<big_int, false> >(report);
}


}

