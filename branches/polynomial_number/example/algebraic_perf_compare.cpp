// This example shows performance comparison for two
// classes for representation algebraic numbers.

#include <iostream>

#include <arageli/algebraic.hpp>
#include <arageli/algebraic/polynomial_number.hpp>
#include <arageli/timer.hpp>
#include <arageli/matrix.hpp>
#include <arageli/gauss.hpp>

int main ()
{
    using namespace Arageli;
    using namespace std;

    try
    {

        timer rational_time(false);
        timer algebraic_time(false);
        timer PolynomialNumber_time(false);

        int rational_iters = 0;
        int algebraic_iters = 0;
        int PolynomialNumber_iters = 0;

        // First compare overhead introduced by algebraic classes to
        // calculations with rational numbers.

        cout << "Please wait while performance statistics is being gathered for Arageli::rational<>..." << flush;
        {

            matrix<rational<> > a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

            auto_timer<> atm(rational_time);
            do
            {
                rref(a, b, q);
                bool is_correct = (b == q*a);
                if(!is_correct)
                    throw "[ ERROR ] a != b*q for rational";
                rational_iters++;
            }while(rational_time.precision() > 0.001);
        }
        cout << " done.\n";

        cout << "Please wait while performance statistics is being gathered for Arageli::algebraic<>..." << flush;
        {

            matrix<algebraic<> > a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

            auto_timer<> atm(algebraic_time);
            do
            {
                rref(a, b, q);
                bool is_correct = (b == q*a);
                if(!is_correct)
                    throw "[ ERROR ] a != b*q for algebraic";
                algebraic_iters++;
            }while(algebraic_time.precision() > 0.001);
        }
        cout << " done.\n";

        cout << "Please wait while performance statistics is being gathered for Arageli::PolynomialNumber..." << flush;
        {

            // The following code doesn't work.
            // Natalia TODO
            //matrix<PolynomialNumber> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

            basis_field baspol("x^2-2");

            matrix<PolynomialNumber> a(3, PolynomialNumber(baspol), fromval), b, q;

            // The following code doesn't work without setting of baspol in the previous line.
            // Natalia TODO
            for(int i = 0; i < 9; ++i)
                a(i/3, i%3) = i-1;

            auto_timer<> atm(PolynomialNumber_time);
            do
            {
                rref(a, b, q);
                bool is_correct = (b == q*a);
                if(!is_correct)
                    throw "[ ERROR ] a != b*q for PolynomialNumber";
                PolynomialNumber_iters++;
            }while(PolynomialNumber_time.precision() > 0.001);
        }
        cout << " done.\n";

        cout
            << "\nOverhead introduced by Arageli::algebraic<>: "
            << algebraic_time.time() / algebraic_iters / (rational_time.time() / rational_iters)
            << " times slower wrt to Arageli::rational<>.\n";

        cout
            << "\nOverhead introduced by Arageli::PolynomialNumber: "
            << PolynomialNumber_time.time() / PolynomialNumber_iters / (rational_time.time() / rational_iters)
            << " times slower wrt to Arageli::rational<>.\n";
    }
    catch(const char* msg)
    {
        cerr << msg << '\n';
        return 1;
    }
}
