#include <iostream>
#include <iomanip>
#include "arageli/arageli.hpp"

using namespace Arageli;

int main ()
{
    // create matrix with rational elements from string
    matrix<rational<> > A = "((21/3, 3, 4), (3335, 6/5, 75), (81, 9, 10/7))";
    std::cout << "A = \n";
    output_aligned(std::cout, A, "|| ", " ||", "  ");

    // calculate inverse matrix
    matrix<rational<> > AInv(inverse(A));
    std::cout << "\ninversion of A = \n";
    output_aligned(std::cout, AInv, "|| ", " ||", "  ");

    // check inverse operation result
    std::cout
        << "\n\nthe inversion is valid: "
        << std::boolalpha << (A*AInv).is_unit();

    // create 3x3 identity matrix
    matrix<rational<> > E(3, eye);
    std::cout << "\nidentity matrix E = \n";
    output_aligned(std::cout, E, "|| ", " ||", "  ");

    // create 3x3 square matrix
    matrix<sparse_polynom<rational<> > >
        SQ(3, sparse_polynom<rational<> >("3/2*x^2"), fromval);
    std::cout << "\nsquare matrix SQ = \n";
    output_aligned(std::cout, SQ, "|| ", " ||", "  ");

    // create 3x4 fromsize matrix
    matrix<double> NSQ(3, 4, fromsize);
    std::cout << "\nnonsquare matrix NSQ = \n";
    output_aligned(std::cout, NSQ, "|| ", " ||", "  ");

    // create 3x3 diagonal matrix
    matrix<char> D(3, 1, diag);
    std::cout << "\ndiagonal matrix D = \n";
    output_aligned(std::cout, D, "|| ", " ||", "  ");

    return 0;
}
