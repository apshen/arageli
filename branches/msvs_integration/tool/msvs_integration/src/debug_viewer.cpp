/*****************************************************************************

    debug_viewer.cpp

    This file is a part of the Arageli MSVS intergation infrastructure,
    a part of Arageli library.

    Copyright (C) 2007--2008 Andrey M. Kamaev
    Copyright (C) 2007--2008 Sergey S. Lyalin
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
    \file debug_viewer.cpp
    \brief <!--ADD A BRIEF FILE DESCRIPTION HERE-->

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include <windows.h>
#include <cstring>
#include <arageli/arageli.hpp>

using namespace std;

BOOL APIENTRY DllMain
(
    HMODULE hModule,
    DWORD  ul_reason_for_call,
    LPVOID lpReserved
)
{
    return TRUE;
}



typedef struct tagDEBUGHELPER
{
    DWORD dwVersion;

    BOOL (WINAPI *ReadDebuggeeMemory)
    (
        struct tagDEBUGHELPER *pThis,
        DWORD dwAddr,
        DWORD nWant,
        VOID* pWhere,
        DWORD *nGot
    );

    // From here only when dwVersion >= 0x20000.

    DWORDLONG (WINAPI *GetRealAddress) (struct tagDEBUGHELPER *pThis);

    BOOL (WINAPI *ReadDebuggeeMemoryEx)
    (
        struct tagDEBUGHELPER *pThis,
        DWORDLONG qwAddr,
        DWORD nWant,
        VOID* pWhere,
        DWORD *nGot
    );

    int (WINAPI *GetProcessorType) (struct tagDEBUGHELPER *pThis);

} DEBUGHELPER;


typedef HRESULT (WINAPI *CUSTOMVIEWER)
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
);


#define ADDIN_API __declspec(dllexport)

extern "C" ADDIN_API HRESULT WINAPI ArgView_big_int
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
);

extern "C" ADDIN_API HRESULT WINAPI ArgView_big_intd
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
);

extern "C" ADDIN_API HRESULT WINAPI ArgView_big_int_big_struct
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
);


ADDIN_API HRESULT WINAPI ArgView_big_int_big_struct
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
)
{
    DWORD nGot = 0;    
    int ret = 0;
    HRESULT result = E_FAIL;

    void *number;
    struct big_struct
    {
        int sign;               // the sign: 0, 1 or -1
        unsigned long *data;    // the storage for digits
        std::size_t len;        // the number of digits
        int refs;               // the number of points to this number
    } number_val;
    unsigned long *num_data = 0;

    /*
    ret = pHelper->ReadDebuggeeMemoryEx(pHelper, dwAddress, sizeof(void*), &number, &nGot);
    if(nGot !=  sizeof(void*)) goto Finalize;
    if(ret) goto Finalize;
    */

    ret = pHelper->ReadDebuggeeMemoryEx(pHelper, dwAddress, sizeof(big_struct), &number_val, &nGot);
    if(nGot !=  sizeof(big_struct) || ret)
    {
        strcpy(pResult,"<Bad Ptr>");
        result = S_OK;
        goto Finalize;
    }

    if(number_val.len != 0)
    {
        if(number_val.len > 10000000) return E_FAIL;
        try
        {
            num_data = new unsigned long[number_val.len];
        }
        catch(...)
        {
            return E_FAIL;
        }

        ret = pHelper->ReadDebuggeeMemoryEx
        (
            pHelper,
            (DWORDLONG)number_val.data,
            sizeof(unsigned long)*number_val.len,
            num_data,
            &nGot
        );

        if(nGot !=  sizeof(unsigned long)*number_val.len || ret)
        {
            strcpy(pResult,"<Unrecognized number>");
            result = S_OK;
            goto Finalize;
        }
    }
    {
        number_val.data = num_data;
        number = &number_val;
        Arageli::big_int* bn;
        bn = (Arageli::big_int*)(&number);
        std::stringstream num_out;
        std::size_t length = bn->length();
        num_out << '[' << length << " bits] ";
        if(length <= 64)
            num_out << (*bn);
        else
            if(is_positive(*bn))
                num_out << "positive";
            else
                num_out << "negative";
        std::string s = num_out.str();
        const char *res = s.c_str();

        if(strlen(res)+1 > max)
        {
            strcpy(pResult,"<Very Big>");
            result = S_OK;
            goto Finalize;
        }
        strcpy(pResult,res);
        result=S_OK;
        goto Finalize;
    }

Finalize:

    delete num_data;
    return result;
}


ADDIN_API HRESULT WINAPI ArgView_big_int
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
)
{
    DWORD nGot = 0;
    int ret = 0;
    HRESULT result = E_FAIL;

    void *number;
    struct big_struct
    {
        int sign;               // the sign: 0, 1 or -1
        unsigned long *data;    // the storage for digits
        std::size_t len;        // the number of digits
        int refs;               // the number of points to this number
    } number_val;
    unsigned long *num_data=0;

    ret = pHelper->ReadDebuggeeMemoryEx(pHelper, dwAddress, sizeof(void*), &number, &nGot);
    if(nGot !=  sizeof(void*)) goto Finalize;
    if(ret) goto Finalize;

    ret = pHelper->ReadDebuggeeMemoryEx(pHelper, (DWORDLONG)number, sizeof(big_struct), &number_val, &nGot);
    if(nGot !=  sizeof(big_struct) || ret)
    {
        strcpy(pResult,"<Bad Ptr>");
        result = S_OK;
        goto Finalize;
    }

    if(number_val.len != 0)
    {
        if(number_val.len > 10000000) return E_FAIL;
        try
        {
            num_data = new unsigned long[number_val.len];
        }
        catch(...)
        {
            return E_FAIL;
        }

        ret = pHelper->ReadDebuggeeMemoryEx
        (
            pHelper,
            (DWORDLONG)number_val.data,
            sizeof(unsigned long)*number_val.len,
            num_data,
            &nGot
        );

        if(nGot !=  sizeof(unsigned long)*number_val.len || ret)
        {
            strcpy(pResult,"<Unrecognized number>");
            result = S_OK;
            goto Finalize;
        }
    }

    {
        number_val.data = num_data;
        number = &number_val;
        Arageli::big_int* bn;
        bn = (Arageli::big_int*)(&number);

        std::stringstream num_out;
        num_out << (*bn);
        std::string s = num_out.str();
        const char *res = s.c_str();

        if(strlen(res)+1 > max)
        {
            strcpy(pResult,"Very big value");
            result = S_OK;
            goto Finalize;
        }
        strcpy(pResult,res);
        result = S_OK;
        goto Finalize;
    }

Finalize:

    delete num_data;
    return result;
}

ADDIN_API HRESULT WINAPI ArgView_big_intd
(
    DWORD dwAddress,
    DEBUGHELPER *pHelper,
    int nBase,
    BOOL bUniStrings,
    char *pResult,
    size_t max,
    DWORD reserved
)
{
    DWORD nGot;
    HRESULT result = E_FAIL;

    //std::ofstream out("argdebug.log",std::ios_base::out | std::ios_base::app);

    std::stringstream num_out;

    //out << "dwAddress=" << dwAddress
    //    << " pHelper=" << pHelper
    //    << " nBase=" << nBase
    //    << " bUniStrings=" << bUniStrings
    //    << " pResult=" << pResult
    //    << " max=" << max
    //    << " reserved=" << reserved << std::endl;
    
    void *number;
    struct big_struct
    {
        int sign;               // the sign: 0, 1 or -1
        unsigned long *data;    // the storage for digits
        std::size_t len;        // the number of digits
        int refs;               // the number of points to this number
    } number_val;

    unsigned long *num_data=0;
    int ret;

    //ret = pHelper->ReadDebuggeeMemoryEx(pHelper, pHelper->GetRealAddress(pHelper), sizeof(void*), &number, &nGot);
    ret = pHelper->ReadDebuggeeMemoryEx(pHelper, dwAddress, sizeof(void*), &number, &nGot);
    //out <<ret << "number="<< number << std::endl;
    if(nGot !=  sizeof(void*)) goto Finalize;
    if(ret) goto Finalize;
        
    ret =pHelper->ReadDebuggeeMemoryEx(pHelper, (DWORDLONG)number, sizeof(big_struct), &number_val, &nGot);

    //ret = pHelper->ReadDebuggeeMemoryEx(pHelper, pHelper->GetRealAddress(pHelper), sizeof(big_struct), &number_val, &nGot);
    //out <<ret << "sign="<< number_val.sign << std::endl;
    //out <<ret << "data="<< number_val.data << std::endl;
    //out <<ret << "len="<< number_val.len << std::endl;
    //out <<ret << "refs="<< number_val.refs << std::endl;

    if(nGot !=  sizeof(big_struct) || ret)
    {
        strcpy(pResult,"<Bad Ptr>");
        result = S_OK;
        goto Finalize;
    }

    if(number_val.len != 0)
    {
        if(number_val.len > 10000000) return E_FAIL;
        num_data = new unsigned long[number_val.len];
        ret = pHelper->ReadDebuggeeMemoryEx
        (
            pHelper,
            (DWORDLONG)number_val.data,
            sizeof(unsigned long)*number_val.len,
            num_data,
            &nGot
        );

        //out <<ret << "data[0]="<< num_data[0] << std::endl;
        //out <<ret << "data[1]="<< num_data[1] << std::endl;

        if(nGot !=  sizeof(unsigned long)*number_val.len || ret)
        {
            strcpy(pResult,"<Unrecognized number>");
            result = S_OK;
            goto Finalize;
        }
    }

    number_val.data = num_data;
    number = &number_val;
    Arageli::big_int* bn;
    bn = (Arageli::big_int*)(&number);

    //out << "!!!value=" << (*bn) << std::endl;

    {
        num_out << (*bn);
        std::string s = num_out.str();
        const char *res = s.c_str();

        //out << res << std::endl;

        if(strlen(res)+1 > max)
        {
            strcpy(pResult,"Very big value");
            result = S_OK;
            goto Finalize;
        }
        strcpy(pResult,res);

        //out    << std::endl;

        result=S_OK;
        goto Finalize;
    }
    
Finalize:
    delete num_data;
    return result;
}
