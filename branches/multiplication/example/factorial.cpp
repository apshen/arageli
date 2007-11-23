#include <iostream>
#include <iomanip>
#include <arageli/arageli.hpp>

using namespace Arageli;
using std::cout;
using std::endl;
using std::boolalpha;


int main ()
{
    const int nbits_small = 120000;
    const int nbits_big_delta = 10000;
    const int nbits_up_big = 360000;
    const int measure_iters = 500;

    for(int i = nbits_small; i <= nbits_up_big; i += nbits_big_delta)
    {
        big_int small = big_int::random_with_length(nbits_small);
        big_int big = big_int::random_with_length(i);

        timer tm;
        for(int j = 0; j < measure_iters; ++j)
        {
            big_int c = small*big;
            //if (c/small != big)
            //{
            //    std::cout << "\n\nERROR!!!\n\n";
            //}
        }
        tm.stop();
        std::cout << i << '\t' << tm.time() << '\n';
    }
}
