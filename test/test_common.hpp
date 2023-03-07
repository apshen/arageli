#ifndef _ARAGELI_TEST_test_common_hpp
#define _ARAGELI_TEST_test_common_hpp 

#include <iostream>
#include <ts/ts.hpp>

#define ARAGELI_TS_ALLEXCEPT_CATCH_REGION_BEGIN    \
{    \
    try    \
    {    \
        ARAGELI_EXCEPT_LOCCTRL_REGION_BEGIN


#define ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END    \
        ARAGELI_EXCEPT_LOCCTRL_REGION_END    \
    }    \
    catch(const ::std::exception& e)    \
    {    \
        ::tout    \
            << "\nException from the standard library:"    \
            << "\n\twhat() = " << e.what();    \
        return resEXCEPT;    \
    }    \
    catch(const ::Arageli::exception& e)    \
    {    \
        ::tout << "\nException from Arageli:\n" << e;    \
        return resEXCEPT;    \
    }    \
    catch(...)    \
    {    \
        ::tout << "\nUnknown exception.";    \
        return resEXCEPT;    \
    }    \
} 

#endif
