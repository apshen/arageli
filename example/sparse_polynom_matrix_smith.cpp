#include <arageli/arageli.hpp>

// Smith's normal diagonal form

using namespace std;
using namespace Arageli;

int main(int argc, char *argv[])
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
