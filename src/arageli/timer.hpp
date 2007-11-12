/*****************************************************************************

    timer.hpp

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

/** \file timer.hpp
    \brief Running time measurement.

    This file contains classes for fixing the execution time.
    Use the timer class to measure and temporary store time intervals.
    Use the auto_timer class to exception-safe starting and stopping time
    measurement based on the timer class.
*/


#ifndef _ARAGELI_timer_hpp_
#define _ARAGELI_timer_hpp_

#include "config.hpp"

#include <ctime>
#include <iostream>

#include "exception.hpp"
#include "_utility.hpp"


namespace Arageli
{

/// An exception, it's thrown if the system time source isn't available.
class time_source_isnot_available :
    public exception
{};

/// An exception, it's thrown in functions requiring stopped timer if it isn't so.
class timer_isnot_stopped :
    public exception
{};


/// Measures execution time by marking the begin and the end time stamps.
/** The timer can be switched on and off several times during its life time.
    All time ranges are accumulated. The timer can tell what the current
    relative precision of measured time is.

    The timer uses some system time provider, for example, std::clock,
    or one that is specific for the target platform.

    We do not recommend to use this timer to measure short segments of
    execution without checking a temporal resolution.
    Use resolution function to obtain the temporal resolution.

    All times are expressed in seconds as a double precision floatting
    point number.

    Some functions may throw an exception of time_source_isnot_available
    type if the system time source isn't available. In this case the timer
    object doesn't work and you cannot use it.
*/
class timer
{
public:

    /// Starts time tracking if turn_on == true (it's true by the default).
    timer (bool turn_on = true) :
        turn_on_m(false),
        duration(0),
        absprec(0)
    {
        if(turn_on)start();
    }

    /// Starts new time interval.
    /** If the timer is already activated, the call doesn't have any effect. */
    void start ();

    /// Ends the current time interval.
    /** Accumulates elapsed time with the previous ones.
        If the timer is already stopped, the call doesn't have any effect. */
    void stop ();

    /// Activation flag. If true, the timer counts time.
    /** Functions start and stop can change this property of the timer. */
    bool is_active () const
    {
        return turn_on_m;
    }

    /// The total elapsed time (in seconds).
    /** The returned value includes all previous time intervals fixed by
        start-stop calls and the present time interval (if the timer is
        activated at the moment).
        The returned value is expressed in seconds.
        Note that the timer cannot measure the time without an error.
        The current exactness of the measured time is provided by precision
        function. */
    double time () const
    {
        return double(clock_time())*resolution();
    }

    /// The minimal amount of time that can be measured (in seconds).
    /** Duration of one tick. The returned value is expressed in seconds. */
    static double resolution ()
    {
        return 1.0/CLOCKS_PER_SEC;
    }

    /// Relative precision of measuring value of the total elapsed time.
    /** Determines the maximum relative error of the value returned by
        the time function. If the timer hasn't ever been activated by
        the moment of calling this function and/or has been activated
        but the accumulated time is 0, the functions returns
        std::numeric_limits<double>::max() value. */
    double precision () const;

    /// Reinitializes the timer. Semantics is the same as for the constructor.
    void reset (bool turn_on = true)
    {
        turn_on_m = false;
        duration = 0;
        absprec = 0;
        if(turn_on)
            start();
    }

private:

    friend std::ostream& operator<< (std::ostream& s, const timer& t);
    friend std::istream& operator>> (std::istream& s, timer& t);

    void init (std::clock_t dur, std::clock_t ap)
    {
        duration = dur;
        absprec = ap;
        turn_on_m = false;
    }

    /// The current elapsed time in ticks.
    std::clock_t clock_time () const;

    std::clock_t start_stamp;    ///< Beginning of the current interval in ticks.
    std::clock_t duration;    ///< Total accumulated time in ticks.
    std::clock_t absprec;    ///< Absolute precision of the total accumulated time in ticks.
    bool turn_on_m;    ///< Activation flag.

};


/// Stores the timer state to a stream.
/** The timer should be stopped.

    This function and operator>> allow easily measure one time region
    separated, for example, by shutdowns of the application.
    So, you save the state before the shutdown and load it after restarting,
    and all accumulated time is included.

    Note that the save and load state functions cannot guarantee platform
    independence. You should use it to save and load on the same platform. */
std::ostream& operator<< (std::ostream& s, const timer& t);


/// Loads the timer state previously stored by operator<< from a stream.
/** The previous state of the timer object will be lost.
    See comments on complementary operator<<. */
std::istream& operator>> (std::istream& s, timer& t);


/// Starts and stops a timer object automatically.
/** It is based on the constructor and destructor automation to start and
    stop a timer. We recomend to use it locally in functions to make
    code exception safe. */
template <typename Timer = timer>
class auto_timer
{
    Timer& tm;

public:

    /// Type of the controlled timer object.
    typedef Timer timer_type;

    /// Starts the timer.
    auto_timer (timer_type& tma) :
        tm(tma)
    {
        tm.start();
    }

    /// Stops the timer.
    ~auto_timer ()
    {
        tm.stop();
    }
};


/// Timing of the series of tasks.
/** The series_timing function is intended to measure time for several
    tasks managed by integer parameter.
    WARNING! It isn't implemented yet. The signature may be changed soon. */
template <typename Timer, typename Tasks, typename Times>
Times& series_timing (Tasks& tasks, Times& times, Timer& timer);


} // namespace Arageli


#endif  // #ifndef _Arageli_timer_hpp_
