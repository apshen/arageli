/*****************************************************************************

    timer.cpp

    This file is a part of Arageli library.

    Copyright (C) 2006--2007 Sergey S. Lyalin
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

/**
    \file timer.cpp
    \brief The timer.hpp file stuff implementation.
*/

#include <limits>
#include "timer.hpp"


namespace Arageli
{


timer::timer (bool turn_on) :
    turn_on_m(false),
    duration(0),
    absprec(0)
{
    first_calibrate();
    if(turn_on)start();
}


void timer::start ()
{
    if(turn_on_m)return;
    start_stamp = std::clock();
    if(start_stamp == std::clock_t(-1))
        throw time_source_isnot_available();
    turn_on_m = true;
}


void timer::stop ()
{
    if(!turn_on_m)return;
    std::clock_t curclock = std::clock();
    if(curclock == std::clock_t(-1))
        throw time_source_isnot_available();
    duration += (curclock - start_stamp);

    ARAGELI_ASSERT_1(is_calibrated());
    absprec += delta;

    turn_on_m = false;
}


std::clock_t timer::clock_time () const
{
    std::clock_t tm = duration;

    if(turn_on_m)
    {
        std::clock_t curclock = std::clock();
        if(curclock == std::clock_t(-1))
            throw time_source_isnot_available();
        tm += (curclock - start_stamp);
    }

    return tm;
}


double timer::precision () const
{
    std::clock_t tm = clock_time();
    std::clock_t curabsprec = absprec;
    if(turn_on_m)
    {
        ARAGELI_ASSERT_1(is_calibrated());
        curabsprec += delta;
    }

    if(tm == 0)
        if(curabsprec == 0)
            return 0;
        else
            return std::numeric_limits<double>::max();
    else
        return double(curabsprec)/double(tm);
}


void timer::calibrate ()
{
    // It is average approximation for duration of
    // the minimum measured time interval.

    const int ncalibs = 10;  // number of calibration runs
    ARAGELI_ASSERT_1(ncalibs >= 1);

    std::clock_t curclock = std::clock();
    if(curclock == std::clock_t(-1))
        throw time_source_isnot_available();

    // Pass the first partial period.
    while(curclock == std::clock());

    std::clock_t startclock = curclock;    // mark start of interval

    // Pass ncalibs whole periods.
    for(int i = 0; i < ncalibs; ++i)
    {
        std::clock_t prevclock = curclock;

        // Wait actively for changing of std::clock returned value.
        do
        {
            curclock = std::clock();
            ARAGELI_ASSERT_1(curclock != -1);
        }while(curclock == prevclock);

        ARAGELI_ASSERT_1(curclock > prevclock);
    }

    delta = (std::clock() - startclock)/ncalibs + 1;    // +1 is for the upper estimate
    sdelta = double(delta)/CLOCKS_PER_SEC;
}


std::clock_t timer::delta = 0;
double timer::sdelta;


std::ostream& operator<< (std::ostream& s, const timer& t)
{
    if(t.is_active())
        throw timer_isnot_stopped();
    s << '(' << t.duration << ',' << t.absprec << ')';
    return s;
}


std::istream& operator>> (std::istream& s, timer& t)
{
    if(_Internal::is_bad_read_literal(s, "("))
        return s;
    std::clock_t duration;
    s >> duration;
    if(!s || _Internal::is_bad_read_literal(s, ","))
        return s;
    std::clock_t absprec;
    s >> absprec;
    if(!s || _Internal::is_bad_read_literal(s, ")"))
        return s;
    t.init(duration, absprec);
    return s;
}


} // namespace Arageli
