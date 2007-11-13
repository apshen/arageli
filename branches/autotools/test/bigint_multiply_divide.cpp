/*****************************************************************************

    test/bigint_multiply_divide.cpp

    This file is a part of the Arageli library.

    Copyright (C) 1999--2007 Nikolai Yu. Zolotykh
    REFERENCE ADDITIONAL COPYRIGHTS HERE
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

// Author : Aleksandrov Vladimir
// Tester : Sidnev Alexey
// Creation date : 20.01.2006
// Modification date: 30.01.2006
// Testing date: 12.02.2006
// Description : * /

//1)Убрано деление на 0
//2)Ошибка теста при T=int. При перемножении 2-х больших int'ов
//big_int z=x*y;
//произведение 2-х значений становится
//больше максимально допустимого. В переменную z записывается неправельный результат.

// Sergey S. Lyalin: I think it's not a test, it's a piss poor.

#include "stdafx.hpp"

//#include <time.h>
//#include <strstream>

using namespace Arageli;


template <class T>
bool big_int_multiply_divide (const T A,const T B)
{
    try
    {
        big_int a=A, b=B, d=a*b;

        big_int x=A, y=B;
        big_int z=x*y;

        if (z!=d)
        {
            tout<<"function multiply_divide failed on 1 check with A="<<A<<", B="<<B<<'\n';
            return true;
        }

        if ((a!=(a/b)*b+a%b)||(b!=(b/a)*a+b%a))
        {
            tout<<"function multiply_divide failed on 2 check with A="<<A<<", B="<<B<<'\n';
            return true;
        }
        return false;

        if ((a!=d/b+d%b)||(b!=d/a+d%a))
        {
            tout<<"function multiply_divide failed on 3 check with A="<<A<<", B="<<B<<'\n';
            return true;
        }
        return false;
    }
    catch(Arageli::exception& e)
    {
        tout << "EXCEPTION: ";
        e.output(tout); tout << '\n';
        return true;
    }
}


bool big_int_multiply_divide_test(int param, int count)
{
    bool fail=false;
    RNG gen(param,1<<16);

        //Message 1 fail |= big_int_multiply_divide(0,0);

        char x=3;
        char y=4;
        fail |= big_int_multiply_divide((signed char) (x),(signed char) (y));
        fail |= big_int_multiply_divide((unsigned char) (x),(unsigned char) (y));
        fail |= big_int_multiply_divide(x,y);

        if(fail)
        {
            tout<<"Function big_int_multiply_divide_test failed on first step.\n";
            return fail;
        }


        for (int k=0; k<count; k++)
        {
            int a=gen.Rand()+1;
            int b=gen.Rand()+1;
            fail |= big_int_multiply_divide(a,b);
            fail |= big_int_multiply_divide(double (a),double (b));
            fail |= big_int_multiply_divide(float (a),float (b));
            fail |= big_int_multiply_divide(long (a),long (b));
            fail |= big_int_multiply_divide((long double) (a),(long double) (b));
            fail |= big_int_multiply_divide((signed short) (a),(signed short) (b));
            fail |= big_int_multiply_divide((unsigned short) (a),(unsigned short) (b));
            fail |= big_int_multiply_divide((unsigned) int (a),(unsigned int) (b));
            fail |= big_int_multiply_divide((signed int) (a),(signed int) (b));

            if(fail)
            {
                tout<<"Function big_int_multiply_divide_test failed on "<<k+1<<" step.\n";
                return fail;
            }
        }

        bool a1=1;
        bool a2=1;
        fail |= big_int_multiply_divide(a1,a2);

        if(fail)
        {
            tout<<"Function big_int_multiply_divide_test failed on last step.\n";
            return fail;
        }

    return false;
}

TEST(big_int,multiply_divide,"Test")
{
    bool fail=big_int_multiply_divide_test(100,1000);

    if(fail) return resFAIL;
    return resOK;

}
