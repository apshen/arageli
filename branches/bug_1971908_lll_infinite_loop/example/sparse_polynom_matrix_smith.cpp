#include <iostream>
#include <iomanip>
#include <arageli/arageli.hpp>

using namespace Arageli;
using std::cout;
using std::endl;
using std::size_t;
using std::boolalpha;


int main ()
{
    matrix< sparse_polynom<rational<> > > A, B, P, Q;
    size_t rk;
    sparse_polynom<rational<> > d;

    A = "((x+1,x-1), (x-2,x+1))";

    cout << "A = " << endl;
    output_aligned(cout, A, "|| ", " ||", " ");
    cout << endl;

    smith(A, B, P, Q, rk, d);

    cout << "B = " << endl;
    output_aligned(cout, B, "|| ", " ||", " ");
    cout << endl;

    cout << "P = " << endl;
    output_aligned(cout, P, "|| ", " ||", " ");
    cout << endl;

    cout << "Q = " << endl;
    output_aligned(cout, Q, "|| ", " ||", " ");
    cout << endl;

    cout << "det(A) = " << d << endl;
    cout << "det(B) = " << det_int(B) << endl;
    cout << "det(P) = " << det_int(P) << endl;
    cout << "det(Q) = " << det_int(Q) << endl;
    cout << "B == P*A*Q: it's " << boolalpha << (B == P*A*Q) << endl;

    return 0;
}
