#include <arageli/arageli.hpp>

using namespace std;
using namespace Arageli;

int main(int argc, char *argv[])
{

    vector<int> a, b;
    a = "(4, 6, 16, 8)";
    b = "(2, 3, 8, 6)";

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;

    cout << "GCD for entries in a = " << gcd(a) << endl;
    cout << "LCM for entries in b = " << lcm(b) << endl;

    cout << "GCD for a and b = " << gcd(a, b) << endl;
    cout << "LCM for a and b = " << lcm(a, b) << endl;

    return 0;
}
