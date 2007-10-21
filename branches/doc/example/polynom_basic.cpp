#include <iostream>
#include <iomanip>

#include <arageli/arageli.hpp>
#include <arageli/polynom.hpp>

using namespace Arageli;

int main ()
{
    polynom<rational<> >
        f,
        g/*
        ,
                r = "1/7*x^13-5*x^8-7*x^5+10"*/
        ;

//        f = "1/7*x^13-5*x^8-7*x^5+10",
//        g = "-11*x^3+22/17*x^2-x";

/*
    std::cout
        << "f(x) = " << f
        << "\ng(x) = " << g
        << "\nf(x)*g(x) = " << f*g
        << "\nf(x)+g(x) = " << f + g
        << "\nf(x)/g(x) = " << f/g
        << "\nf(x)%g(x) = " << f%g
        << "\nGCD(f(x), g(x)) = " << gcd(f, g)
        << "\ndividing is valid: " << std::boolalpha << ((f/g)*g + f%g == f);
*/

    return 0;
}
