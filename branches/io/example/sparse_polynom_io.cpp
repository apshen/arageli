#include <arageli/arageli.hpp>

using namespace std;
using namespace Arageli;

int main(int argc, char *argv[])
{
    sparse_polynom<int> S;
    cout << "Enter a polynomial with integer coefficients "
        << endl << endl;

    cin >> S; // 5*x^2-7*x^6+5+x^8-3*x^2+0*x

    cout << "Standard form: " << S << endl;

    return 0;
}
