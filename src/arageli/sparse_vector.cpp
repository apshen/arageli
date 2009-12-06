/*****************************************************************************

    sparse_vector.cpp

    This file is a part of the Arageli library.

    Copyright (C)2007--2008 Sergey S. Lyalin
    Copyright (C)2008--2009 Valentin K. Kubarev
    University of Nizhni Novgorod, Russia

    The Arageli Library is free software;you can redistribute it and/or
    modify it under the terms of the GNU General Public License version 2
    as published by the Free Software Foundation.

    The Arageli Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY;without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program;if not, write to the Free Software Foundation, Inc.,
    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

    We are also open for dual licensing for the whole library or
    for its particular part. If you are interested to get the library
    in this way, i.e. not under the GNU General Public License,
    please contact Arageli Support Service support.arageli@gmail.com.

*****************************************************************************/

/**
    \file sparse_vector.cpp
    \brief The sparse_vector.hpp file stuff implementation.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "config.hpp"

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SPARSE_VECTOR)

// REFERENCE ADDITIONAL HEADERS HERE

#include "sparse_vector.hpp"


namespace Arageli
{

// Namespace of sparse vector representation
namespace sprsvc_rep
{

/* Class representation of sparse vector functions implementation
    Sparse vector on vector of pairs */
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
sparsevec_vecpair_t<T,REFCNT>::sparsevec_vecpair_t(size_t _n):
        n(_n)
#ifdef ARAGELI_DEBUG_LEVEL_3
    , ordered(true)
#endif

{
}

template <typename T, bool REFCNT>
sparsevec_vecpair_t<T,REFCNT>::~sparsevec_vecpair_t()
{
}

template <typename T, bool REFCNT>
const size_t sparsevec_vecpair_t<T,REFCNT>::real_index(size_t k) const
{
    if(real_size())
    {
        size_t first = 0;
        size_t last  = real_size()-1;
            size_t middle;
        while(first<last)
        {
            middle=(first+last)/2+1;
            if (k==rep[middle].first)
            {
                first=middle;
                last=middle;
            }
            if (k>rep[middle].first)
            {
                first=middle+1;
            }
            else
            {
                last=middle-1;
            }
        }
        return first;
    }
    return 0;
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::ins(size_t i, const T &t)
{
    if(real_size())
    {
        size_t k = real_size();
        if(i<rep[real_size()-1].first)
        {
            k = real_index(i);
            if(rep[k].first<i)
            {
                ++k;
            }
        }
        else if(i==rep[real_size()-1].first)
        {
            rep[real_size()-1].second = t;
            return;
        }
        if(k<real_size()&& rep[k].first==i)
        {
            rep[k].second = t;
        }
        else
        {
            Rep new_rep;
            new_rep.assign(real_size()+1,Elem(0,null<T>()));
            size_t j=0;
            while(j<k)
            {
                new_rep[j].first = rep[j].first;
                new_rep[j].second = rep[j].second;
                ++j;
            }
            new_rep[k].first = i;
            new_rep[k].second = t;
            while(j<real_size())
            {
                new_rep[j+1].first = rep[j].first;
                new_rep[j+1].second = rep[j].second;
                ++j;
            }
            rep.swap(new_rep);
        }
    }
    else
    {
        rep.resize(1);
        rep[0].first = i;
        rep[0].second = t;
    }
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::del(size_t i)
{
    if(real_size())
    {
        size_t k = real_size();
        if(i<rep[real_size()-1].first)
        {
            k = real_index(i);
            if(rep[k].first==i)
            {
                for(size_t i=k;i<real_size()-1;++i)
                {
                    rep[k].first = rep[k+1].first;
                    rep[k].second = rep[k+1].second;
                }
                rep.resize(real_size()-1);
            }
        }
    }
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::raw_set_elem(size_t real_ind, size_t ind, const T & t)
{
    rep[real_ind].first=ind;
    rep[real_ind].second=t;
}

template <typename T, bool REFCNT>
inline void sparsevec_vecpair_t<T,REFCNT>::push_back(size_t ind, const T & t)
{
    ARAGELI_ASSERT_0(ind<n);

#ifdef ARAGELI_DEBUG_LEVEL_3
    if(rep.size()&& rep[rep.size()-1].first>=ind)
    {
        ordered = false;
    }
#endif

    rep.push_back(Elem(ind,t));
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::sorting(Rep &A, int first, int last)
{
    int h, index, countr_index;
    bool condition;
    if (first < last)
    {
        h=-1;
        index=first;countr_index=last;condition=1;
        while (index!=countr_index)
        {
            if ((A[index].first > A[countr_index].first)==condition)
            {
                h=-h;
                std::swap(index,countr_index);
                std::swap(A[index],A[countr_index]);
                condition=!condition;
            }
            countr_index=countr_index+h;
        }
        sorting(A,index+1,last);
        sorting(A,first,index-1);
    }
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::remove_useless()
{
    size_t i=0;
    size_t l=0;
    for(;i<rep.size();++i)
    {
        if(is_null(rep[i].second)|| (i>0 && rep[i-1].first==rep[i].first))
        {
            continue;
        }
        if(i!=l)
        {
            rep[l].first=rep[i].first;
            rep[l].second=rep[i].second;
        }
        ++l;
    }
    rep.resize(l);
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::resize(size_t _n)
{
    if(_n<n && rep.size())
    {
        size_t count=0;
        for(size_t i=0;i<rep.size()&& rep[i].first<_n;++i)
        {
            ++count;
        }
        rep.resize(count);
    }
    n=_n;

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = check_ordered();
#endif
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::reserve_mem(size_t elem_count)
{
    if(elem_count==rep.capacity())
    {
        return;
    }
    else if(elem_count<rep.capacity())
    {
        size_t count = rep.size();
        if(elem_count<count)
        {
            count=elem_count;
        }
        rep.resize(0);
        Rep temp;
        rep.swap(temp);
        rep.resize(count);
        for(size_t i=0;i<count;++i)
        {
            rep[i] = temp[i];
        }
    }
    rep.reserve(elem_count);
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::resize_mem_struct(size_t elem_count)
{
    if(elem_count<rep.size())
    {
        Rep temp;
        rep.swap(temp);
        rep.resize(elem_count);
        for(size_t i=0;i<elem_count;++i)
        {
            rep[i].first = temp[i].first;
            rep[i].second = temp[i].second;
        }
    }
    else if(elem_count!=rep.size())
    {
        rep.resize(elem_count);
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = check_ordered();
#endif
}

template <typename T, bool REFCNT>
void sparsevec_vecpair_t<T,REFCNT>::clear()
{
    reserve_mem(0);

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = true;
#endif
}

template <typename T, bool REFCNT>
const T &    sparsevec_vecpair_t<T,REFCNT>::get(size_t i) const
{
    ARAGELI_ASSERT_0(i>=0 && i<n);

    int j=(int)real_index(i);
    if(j<real_size())
    {
        if(rep[j].first==i)
        {
            return rep[j].second;
        }
    }

    return null<T>();
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/* Class representation of sparse vector functions implementation
    Sparse vector on pair of vectors */
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
sparsevec_pairvec_t<T,REFCNT>::sparsevec_pairvec_t(size_t _n):
        n(_n)
#ifdef ARAGELI_DEBUG_LEVEL_3
    , ordered(true)
#endif
{
}

template <typename T, bool REFCNT>
sparsevec_pairvec_t<T,REFCNT>::~sparsevec_pairvec_t()
{
}

template <typename T, bool REFCNT>
const size_t sparsevec_pairvec_t<T,REFCNT>::real_index(size_t k) const
{
    if(real_size())
    {
        size_t
        middle,
        first=0,
        last=real_size()-1;
        while(first<last)
        {
            middle=(first+last)/2+1;
            if (k==rep.first[middle])
            {
                first=middle;
                last=middle;
            }
            if (k>rep.first[middle])first=middle+1;
            else last=middle-1;
        }
        return first;
    }
    return 0;
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::ins(size_t i, const T &t)
{
    ARAGELI_ASSERT_0(i<n)

    if(real_size())
    {
        size_t k = real_size();
        if(i<rep.first[real_size()-1])
        {
            k = real_index(i);
            if(rep.first[k]<i)
            {
                ++k;
            }
        }
        else if(i==rep.first[real_size()-1])
        {
            rep.second[real_size()-1] = t;
            return;
        }
        if(k<real_size()&& rep.first[k]==i)
        {
            rep.second[k] = t;
        }
        else
        {
            Rep new_rep;
            new_rep.first.assign(real_size()+1,0);
            new_rep.second.assign(real_size()+1,null<T>());
            size_t j=0;
            while(j<k)
            {
                new_rep.first[j] = rep.first[j];
                new_rep.second[j] = rep.second[j];
                ++j;
            }
            new_rep.first[k] = i;
            new_rep.second[k] = t;
            while(j<real_size())
            {
                new_rep.first[j+1] = rep.first[j];
                new_rep.second[j+1] = rep.second[j];
                ++j;
            }
            rep.swap(new_rep);
        }
    }
    else
    {
        rep.first.resize(1);
        rep.second.resize(1);
        rep.first[0] = i;
        rep.second[0] = t;
    }
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::del(size_t i)
{
    if(real_size())
    {
        size_t k = real_size();
        if(i<rep.first[real_size()-1])
        {
            k = real_index(i);
            if(rep.first[k]==i)
            {
                for(size_t i=k;i<real_size()-1;++i)
                {
                    rep.first[k] = rep.first[k+1];
                    rep.second[k] = rep.second[k+1];
                }
                size_t size = real_size()-1;
                rep.first.resize(size);
                rep.second.resize(size);
            }
        }
    }
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::raw_set_elem(size_t real_ind, size_t ind, const T & t)
{
    rep.first[real_ind]=ind;
    rep.second[real_ind]=t;
}

template <typename T, bool REFCNT>
inline void sparsevec_pairvec_t<T,REFCNT>::push_back(size_t ind, const T & t)
{
    ARAGELI_ASSERT_0(ind<n);

#ifdef ARAGELI_DEBUG_LEVEL_3
    if(rep.first.size()&& rep.first[rep.first.size()-1]>=ind)
    {
        ordered = false;
    }
#endif

    rep.first.push_back(ind);
    rep.second.push_back(t);
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::sorting(Rep &A, int first, int last)
{
    int h;
    int index, countr_index;
    bool condition;
    if (first < last)
    {
        h=-1;
        index=first;
        countr_index=last;
        condition=1;
        while (index!=countr_index)
        {
            if ((A.first[index] > A.first[countr_index])==condition)
            {
                h=-h;
                std::swap(index,countr_index);
                std::swap(A.first[index],A.first[countr_index]);
                std::swap(A.second[index],A.second[countr_index]);
                condition=!condition;
            }
            countr_index=countr_index+h;
        }
        sorting(A,index+1,last);
        sorting(A,first,index-1);
    }
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::remove_useless()
{
    size_t i=0;
    size_t l=0;

    for(;i<rep.first.size();++i)
    {
        if(is_null(rep.second[i])|| (i>0 && rep.first[i-1]==rep.first[i]))
        {
            continue;
        }
        if(i!=l)
        {
            rep.first[l]=rep.first[i];
            rep.second[l]=rep.second[i];
        }
        ++l;
    }
    rep.first.resize(l);
    rep.second.resize(l);
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::resize(size_t _n)
{
    if(_n<n && rep.first.size())
    {
        size_t count=0;
        for(size_t i=0;i<rep.first.size()&& rep.first[i]<_n;++i)
        {
            ++count;
        }
        rep.first.resize(count);
        rep.second.resize(count);
    }
    n=_n;

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = check_ordered();
#endif
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::reserve_mem(size_t elem_count)
{
    if(elem_count==rep.first.capacity())
    {
        return;
    }
    else if(elem_count<rep.first.capacity())
    {
        size_t count = rep.first.size();
        if(elem_count<count)
        {
            count=elem_count;
        }
        rep.first.resize(0);
        rep.second.resize(0);
        Rep temp;
        rep.first.swap(temp.first);
        rep.second.swap(temp.second);
        rep.first.resize(count);
        rep.second.resize(count);
        for(size_t i=0;i<count;++i)
        {
            rep.first[i] = temp.first[i];
            rep.second[i] = temp.second[i];
        }
    }
    rep.first.reserve(elem_count);
    rep.second.reserve(elem_count);
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::resize_mem_struct(size_t elem_count)
{
    if(elem_count<rep.first.size())
    {
        Rep temp;
        rep.first.swap(temp.first);
        rep.second.swap(temp.second);
        rep.first.resize(elem_count);
        rep.second.resize(elem_count);
        for(size_t i=0;i<elem_count;++i)
        {
            rep.first[i] = temp.first[i];
            rep.second[i] = temp.second[i];
        }
    }
    else if(elem_count!=rep.first.size())
    {
        rep.first.resize(elem_count);
        rep.second.resize(elem_count);
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = check_ordered();
#endif
}

template <typename T, bool REFCNT>
void sparsevec_pairvec_t<T,REFCNT>::clear()
{
    reserve_mem(0);

#ifdef ARAGELI_DEBUG_LEVEL_3
    ordered = true;
#endif
}

template <typename T, bool REFCNT>
const T &    sparsevec_pairvec_t<T,REFCNT>::get(size_t i) const
{
    ARAGELI_ASSERT_0(i>=0 && i<n);

    int j=(int)real_index(i);
    if(j<real_size())
    {
        if(rep.first[j]==i)
        {
            return rep.second[j];
        }
    }

    return null<T>();
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/* Class representation of sparse vector functions implementation
    Sparse vector on map of pairs */
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
sparsevec_mappair_t<T,REFCNT>::sparsevec_mappair_t(size_t _n):
        n(_n)
{
}

template <typename T, bool REFCNT>
sparsevec_mappair_t<T,REFCNT>::~sparsevec_mappair_t()
{
}

template <typename T, bool REFCNT>
void sparsevec_mappair_t<T,REFCNT>::ins(size_t i, const T &t)
{
    if(rep.size()==0 || i>(--rep.end())->first)
    {
        rep.insert(rep.end(),Elem(i,t));
    }
    else
    {
        rep[i]=t;
    }
}

template <typename T, bool REFCNT>
void sparsevec_mappair_t<T,REFCNT>::del(size_t i)
{
    Rep::iterator it = rep.find(i);
    if(it != rep.end())
    {
        rep.erase(it);
    }
}


template <typename T, bool REFCNT>
void sparsevec_mappair_t<T,REFCNT>::resize(size_t _n)
{
    if(_n<n && rep.size())
    {
        Rep::iterator it = rep.begin();
        Rep::iterator it_end = rep.end();
        while(it!=it_end)
        {
            if(it->first<_n)
            {
                it=rep.erase(it);
            }
        }
    }

    n=_n;
}

template <typename T, bool REFCNT>
void sparsevec_mappair_t<T,REFCNT>::reserve_mem(size_t _n)
{
}

template <typename T, bool REFCNT>
void sparsevec_mappair_t<T,REFCNT>::clear()
{
    rep.clear();
}

template <typename T, bool REFCNT>
const T & sparsevec_mappair_t<T,REFCNT>::get(size_t i) const
{
    if(i>=0 && i<n)
    {
        Rep::const_iterator it = rep.find(i);
        if(it!=rep.end())
        {
            return it->second;
        }
    }
    return null<T>();
}

template <typename T, bool REFCNT>
inline void sparsevec_mappair_t<T,REFCNT>::push_back(size_t ind, const T & t)
{
    ARAGELI_ASSERT_0(ind<n);
    ARAGELI_ASSERT_0(rep.size()==0 || ind>(--rep.end())->first);

    rep.insert(rep.end(), Elem(ind, t));
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}// namespace sprsvc_rep

namespace sprsvc_comp
{
using namespace sprsvc_rep;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Implementation of some template function for sparse vector on vector of pairs
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
void copy_to(sparsevec_vecpair_t<T,REFCNT> &dst, const sparsevec_vecpair_t<T,REFCNT> &src)
{
    dst.n=src.n;
    dst.rep=src.rep;

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.ordered = src.ordered;
#endif
}

template <typename T, bool REFCNT>
void copy_to(sparsevec_vecpair_t<T,REFCNT> &dst, const Arageli::vector<T> &src)
{
    dst.n=src.size();
    size_t count=0;
    for(int i=0;i<dst.n;i++)
    {
        if(!is_null(src[i]))
        {
            ++count;
        }
    }
    dst.rep.resize(count);
    for (size_t i=0,k=0;i<dst.n;++i)
    {
        if(!is_null(src[i]))
        {
            dst.rep[k].first = i;
            dst.rep[k++].second = src[i];
        }
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.ordered = true;
#endif
}

template <typename T, bool REFCNT>
void add(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    res.resize_mem_struct(struct_union(sv1, sv2));
    // numerical part
    typedef sparsevec_vecpair_t<T,REFCNT>::Elem Elem;
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep[j].first < sv2.rep[i].first)
        {
            res.rep[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep[j].first == sv2.rep[i].first)
        {
            res.rep[count++]=Elem(sv1.rep[j].first,sv1.rep[j].second+sv2.rep[i].second);
            ++j;
        }
        else
        {
            res.rep[count++] = Elem(sv2.rep[i].first,sv2.rep[i].second);
        }
    }

    while (j < sv1.real_size())
    {
        res.rep[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
        ++j;
    }
    /*********************************************************************************************/

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = true;
#endif
}

template <typename T, bool REFCNT>
void sub(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    res.resize_mem_struct(struct_union(sv1, sv2));
    // numerical part
    typedef sparsevec_vecpair_t<T,REFCNT>::Elem Elem;
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep[j].first < sv2.rep[i].first)
        {
            res.rep[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep[j].first == sv2.rep[i].first)
        {
            res.rep[count++]=Elem(sv1.rep[j].first,sv1.rep[j].second-sv2.rep[i].second);
            ++j;
        }
        else
        {
            res.rep[count++] = Elem(sv2.rep[i].first,-sv2.rep[i].second);
        }
    }

    while (j < sv1.real_size())
    {
        res.rep[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
        ++j;
    }
    /*********************************************************************************************/

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = true;
#endif
}

template <typename T, bool REFCNT>
void mul(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv, const T &val)
{
    res.n = sv.n;
    res.resize_mem_struct(sv.real_size());
    for (size_t j = 0;j < sv.real_size();++j)
    {
        res.rep[j].first = sv.rep[j].first;
        res.rep[j].second = sv.rep[j].second*val;
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = sv.ordered;
#endif
}

template <typename T, bool REFCNT>
void div(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv, const T &val)
{
    res.n = sv.n;
    res.resize_mem_struct(sv.real_size());
    for (size_t j = 0;j < sv.real_size();++j)
    {
        res.rep[j].first = sv.rep[j].first;
        res.rep[j].second = sv.rep[j].second/val;
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = sv.ordered;
#endif
}

template <typename T, bool REFCNT>
void dot(T &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res = null<T>();
    size_t
        i=0,
        j=0;
    size_t count=0;
    for (;i < sv2.real_size()&& j < sv1.real_size();)
    {
        if(sv1.rep[j].first == sv2.rep[i].first)
        {
            res += sv1.rep[j].second * sv2.rep[i].second;
            ++j;
            ++i;
        }
        else if(sv1.rep[j].first < sv2.rep[i].first)
        {
            ++j;
        }
        else
        {
            ++i;
        }
    }
}

template <typename T, bool REFCNT>
void dot(T &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const Arageli::vector<T> &vec)
{
    ARAGELI_ASSERT_0( sv1.size()==vec.size());
    res = null<T>();
    size_t i=0;

    for (;i < sv1.real_size();++i)
    {
        res += sv1.rep[i].second * vec[sv1.rep[i].first];
    }
}

template <typename T, bool REFCNT>
void add_in(sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    sparsevec_vecpair_t<T,REFCNT>::Rep temp;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    temp.resize(struct_union(sv1, sv2));
    // numerical part
    typedef sparsevec_vecpair_t<T,REFCNT>::Elem Elem;
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep[j].first < sv2.rep[i].first)
        {
            temp[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep[j].first == sv2.rep[i].first)
        {
            temp[count++]=Elem(sv1.rep[j].first,sv1.rep[j].second+sv2.rep[i].second);
            ++j;
        }
        else
        {
            temp[count++] = Elem(sv2.rep[i].first,sv2.rep[i].second);
        }
    }

    while (j < sv1.real_size())
    {
        temp[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
        ++j;
    }
    /*********************************************************************************************/
    temp.swap(sv1.rep);
}

template <typename T, bool REFCNT>
void sub_in(sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    sparsevec_vecpair_t<T,REFCNT>::Rep temp;
    size_t
        i=0,
        j=0,
        count=0;

    /****************************************symbolic part****************************************/
    temp.resize(struct_union(sv1, sv2));
    /****************************************numerical part***************************************/
    typedef sparsevec_vecpair_t<T,REFCNT>::Elem Elem;
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep[j].first < sv2.rep[i].first)
        {
            temp[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep[j].first == sv2.rep[i].first)
        {
            temp[count++]=Elem(sv1.rep[j].first,sv1.rep[j].second-sv2.rep[i].second);
            ++j;
        }
        else
        {
            temp[count++] = Elem(sv2.rep[i].first,-sv2.rep[i].second);
        }
    }

    while (j < sv1.real_size())
    {
        temp[count++] = Elem(sv1.rep[j].first,sv1.rep[j].second);
        ++j;
    }
    /*********************************************************************************************/
    temp.swap(sv1.rep);
}

template <typename T, bool REFCNT>
void mul_in(sparsevec_vecpair_t<T,REFCNT> &sv, const T &val)
{
    for (size_t j = 0;j < sv.real_size();++j)
    {
        sv.rep[j].first     = sv.rep[j].first;
        sv.rep[j].second    = sv.rep[j].second*val;
    }
}

template <typename T, bool REFCNT>
void div_in(sparsevec_vecpair_t<T,REFCNT> &sv, const T &val)
{
    for (size_t j = 0;j < sv.real_size();++j)
    {
        sv.rep[j].first     = sv.rep[j].first;
        sv.rep[j].second    = sv.rep[j].second/val;
    }
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Implementation of some template function for sparse vector on pair of vectors
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
void copy_to(sparsevec_pairvec_t<T,REFCNT> &dst, const sparsevec_pairvec_t<T,REFCNT> &src)
{
    dst.n=src.n;
    dst.rep=src.rep;

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.ordered=src.ordered;
#endif
}

template <typename T, bool REFCNT>
void copy_to(sparsevec_pairvec_t<T,REFCNT> &dst, const Arageli::vector<T> &src)
{
    dst.n=src.size();
    size_t count=0;
    for(int i=0;i<dst.n;i++)
    {
        if(!is_null(src[i]))
        {
            ++count;
        }
    }
    dst.resize_mem_struct(count);
    for (size_t i=0,k=0;i<dst.n;++i)
    {
        if(!is_null(src[i]))
        {
            dst.rep.first[k] = i;
            dst.rep.second[k++] = src[i];
        }
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.ordered=true;
#endif
}

template <typename T, bool REFCNT>
void add(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    res.resize_mem_struct(struct_union(sv1, sv2));
    // numerical part
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep.first[j] < sv2.rep.first[i])
        {
            res.rep.first[count] = sv1.rep.first[j];
            res.rep.second[count++] = sv1.rep.second[j];
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep.first[j] == sv2.rep.first[i])
        {
            res.rep.first[count] = sv1.rep.first[j];
            res.rep.second[count++] = sv1.rep.second[j]+sv2.rep.second[i];
            ++j;
        }
        else
        {
            res.rep.first[count] = sv2.rep.first[i];
            res.rep.second[count++] = sv2.rep.second[i];
        }
    }

    while (j < sv1.real_size())
    {
        res.rep.first[count] = sv1.rep.first[j];
        res.rep.second[count++] = sv1.rep.second[j];
        ++j;
    }
    // +++++++++

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = true;
#endif
}

template <typename T, bool REFCNT>
void sub(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    res.resize_mem_struct(struct_union(sv1, sv2));
    // numerical part
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep.first[j] < sv2.rep.first[i])
        {
            res.rep.first[count] = sv1.rep.first[j];
            res.rep.second[count++] = sv1.rep.second[j];
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep.first[j] == sv2.rep.first[i])
        {
            res.rep.first[count] = sv1.rep.first[j];
            res.rep.second[count++] = sv1.rep.second[j]-sv2.rep.second[i];
            ++j;
        }
        else
        {
            res.rep.first[count] = sv2.rep.first[i];
            res.rep.second[count++] = -sv2.rep.second[i];
        }
    }

    while (j < sv1.real_size())
    {
        res.rep.first[count] = sv1.rep.first[j];
        res.rep.second[count++] = sv1.rep.second[j];
        ++j;
    }
    // +++++++++

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = true;
#endif
}

template <typename T, bool REFCNT>
void mul(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv, const T &val)
{
    res.n = sv.n;
    res.resize_mem_struct(sv.real_size());
    for (size_t j = 0;j < sv.real_size();++j)
    {
        res.rep.first[j] = sv.rep.first[j];
        res.rep.second[j] = sv.rep.second[j]*val;
    };

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = sv.ordered;
#endif
}

template <typename T, bool REFCNT>
void div(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv, const T &val)
{
    res.n = sv.n;
    res.resize_mem_struct(sv.real_size());
    for (size_t j = 0;j < sv.real_size();++j)
    {
        res.rep.first[j] = sv.rep.first[j];
        res.rep.second[j] = sv.rep.second[j]/val;
    };

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.ordered = sv.ordered;
#endif
}

template <typename T, bool REFCNT>
void dot(T &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res = null<T>();
    size_t
        i=0,
        j=0;
    size_t count=0;

    for (;i < sv2.real_size()&& j < sv1.real_size();)
    {
        if(sv1.rep.first[j] == sv2.rep.first[i])
        {
            res += sv1.rep.second[j] * sv2.rep.second[i];
            ++j;
            ++i;
        }
        else if(sv1.rep.first[j] < sv2.rep.first[i])
        {
            ++j;
        }
        else
        {
            ++i;
        }
    }
}

template <typename T, bool REFCNT>
void dot(T &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const Arageli::vector<T> &vec)
{
    ARAGELI_ASSERT_0( sv1.size()==vec.size());

    res = null<T>();
    size_t i=0;
    for (;i < sv1.real_size();++i)
    {
        res += sv1.rep.second[i] * vec[sv1.rep.first[i]];
    }

}

template <typename T, bool REFCNT>
void add_in(sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    sparsevec_pairvec_t<T,REFCNT>::Indeces temp_in;
    sparsevec_pairvec_t<T,REFCNT>::Elems temp_el;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    count = struct_union(sv1, sv2);
    temp_in.resize(count);
    temp_el.resize(count);
    count=0;
    // numerical part
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep.first[j] < sv2.rep.first[i])
        {
            temp_in[count] = sv1.rep.first[j];
            temp_el[count++] = sv1.rep.second[j];
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep.first[j] == sv2.rep.first[i])
        {
            temp_in[count] = sv1.rep.first[j];
            temp_el[count++] = sv1.rep.second[j]+sv2.rep.second[i];
            ++j;
        }
        else
        {
            temp_in[count] = sv2.rep.first[i];
            temp_el[count++] = sv2.rep.second[i];
        }
    }

    while (j < sv1.real_size())
    {
        temp_in[count] = sv1.rep.first[j];
        temp_el[count++] = sv1.rep.second[j];
        ++j;
    }
    // +++++++++++++++
    temp_in.swap(sv1.rep.first);
    temp_el.swap(sv1.rep.second);

}

template <typename T, bool REFCNT>
void sub_in(sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    sparsevec_pairvec_t<T,REFCNT>::Indeces temp_in;
    sparsevec_pairvec_t<T,REFCNT>::Elems temp_el;
    size_t
        i=0,
        j=0,
        count=0;

    // symbolic part
    count = struct_union(sv1, sv2);
    temp_in.resize(count);
    temp_el.resize(count);
    count=0;
    // numerical part
    for (;i < sv2.real_size();++i)
    {
        while (j < sv1.real_size()&& sv1.rep.first[j] < sv2.rep.first[i])
        {
            temp_in[count] = sv1.rep.first[j];
            temp_el[count++] = sv1.rep.second[j];
            ++j;
        }

        if (j < sv1.real_size()&& sv1.rep.first[j] == sv2.rep.first[i])
        {
            temp_in[count] = sv1.rep.first[j];
            temp_el[count++] = sv1.rep.second[j]-sv2.rep.second[i];
            ++j;
        }
        else
        {
            temp_in[count] = sv2.rep.first[i];
            temp_el[count++] = -sv2.rep.second[i];
        }
    }

    while (j < sv1.real_size())
    {
        temp_in[count] = sv1.rep.first[j];
        temp_el[count++] = sv1.rep.second[j];
        ++j;
    }
    // +++++++++++++++
    temp_in.swap(sv1.rep.first);
    temp_el.swap(sv1.rep.second);
}

template <typename T, bool REFCNT>
void mul_in(sparsevec_pairvec_t<T,REFCNT> &sv, const T &val)
{
    for (size_t j = 0;j < sv.real_size();++j)
    {
        sv.rep.first[j]     = sv.rep.first[j];
        sv.rep.second[j]    = sv.rep.second[j]*val;
    }
}

template <typename T, bool REFCNT>
void div_in(sparsevec_pairvec_t<T,REFCNT> &sv, const T &val)
{
    for (size_t j = 0;j < sv.real_size();++j)
    {
        sv.rep.first[j]     = sv.rep.first[j];
        sv.rep.second[j]    = sv.rep.second[j]/val;
    }
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// All sparse vector template function implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename A, typename B>
void copy_to(A &dst, const B &src)
{
    dst.n=src.n;
    dst.clear();
    if(!dst.dynamic())
    {
        dst.reserve_mem(src.real_size());
    }
    B::const_indraw_iterator it_src = src.begin();
    B::const_indraw_iterator it_src_end = src.end();

    for (;it_src != it_src_end;++it_src)
    {
        dst.push_back(it_src.ind(),it_src.el());
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.is_ordered(src.is_ordered());
#endif
}

template <typename T, typename A>
void copy_to(A &dst, const Arageli::vector<T> &src)
{
    dst.n=src.size();
    dst.clear();
    if(!dst.dynamic())
    {
        size_t count=0;
        for(int i=0;i<dst.n;i++)
        {
            if(!is_null(src[i]))
            {
                ++count;
            }
        }
        dst.reserve_mem(count);
    }
    for (size_t i=0;i<dst.n;++i)
    {
        if(!is_null(src[i]))
        {
            dst.push_back(i,src[i]);
        }
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    dst.is_ordered(true);
#endif
}

template <typename T, typename A>
void copy_to(Arageli::vector<T> &dst, const A &src)
{
    A::const_indraw_iterator it_src = src.begin();
    A::const_indraw_iterator it_src_end = src.end();

    dst.assign(src.size(), null<T>());

    for(; it_src!=it_src_end; ++it_src)
    {
        dst[it_src.ind()] = it_src.el();
    }
}

template <typename A, typename B>
bool equal(const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    if(sv1.n==sv2.n)
    {
        A::const_indraw_iterator it_sv1     = sv1.begin();
        A::const_indraw_iterator it_sv1_end = sv1.end();
        B::const_indraw_iterator it_sv2     = sv2.begin();
        B::const_indraw_iterator it_sv2_end = sv2.end();
        while(true)
        {
            while(it_sv1!=it_sv1_end && is_null(it_sv1.el()))
            {
                ++it_sv1;
            }
            while(it_sv2!=it_sv2_end && is_null(it_sv2.el()))
            {
                ++it_sv2;
            }
            if(it_sv2==it_sv2_end)
            {
                if(it_sv1==it_sv1_end)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            if(it_sv1==it_sv1_end)
            {
                return false;
            }
            if (it_sv1.ind()!= it_sv2.ind()|| it_sv1.el()!= it_sv2.el())
            {
                return false;
            }
            ++it_sv1;
            ++it_sv2;

        }
        return true;
    }
    return false;
}

template <typename A, typename B>
bool non_equal(const A &sv1, const B &sv2)
{
    return !equal(sv1, sv2);
}

template <typename A, typename B>
void add(A &res, const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    res.clear();
    if(!res.dynamic())
    {
        //symbolic part
        res.reserve_mem(struct_union(sv1, sv2));
    }
    A::const_indraw_iterator it_sv1     = sv1.begin();
    A::const_indraw_iterator it_sv1_end = sv1.end();
    B::const_indraw_iterator it_sv2     = sv2.begin();
    B::const_indraw_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            res.push_back(it_sv1.ind(), it_sv1.el());
            ++it_sv1;
        }

        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            res.push_back(it_sv1.ind(), it_sv1.el()+ it_sv2.el());
            ++it_sv1;
        }
        else
        {
            res.push_back(it_sv2.ind(), it_sv2.el());
        }
    }

    while (it_sv1 != it_sv1_end)
    {
        res.push_back(it_sv1.ind(), it_sv1.el());
        ++it_sv1;
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.is_ordered(true);
#endif
}

template <typename A, typename B>
void sub(A &res, const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res.n=sv1.n;
    // symbolic part
    res.clear();
    if(!res.dynamic())
    {
        res.reserve_mem(struct_union(sv1, sv2));
    }
    A::const_indraw_iterator it_sv1     = sv1.begin();
    A::const_indraw_iterator it_sv1_end = sv1.end();
    B::const_indraw_iterator it_sv2     = sv2.begin();
    B::const_indraw_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            res.push_back(it_sv1.ind(), it_sv1.el());
            ++it_sv1;
        }

        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            res.push_back(it_sv1.ind(), it_sv1.el()- it_sv2.el());
            ++it_sv1;
        }
        else
        {
            res.push_back(it_sv2.ind(), -it_sv2.el());
        }
    }

    while (it_sv1 != it_sv1_end)
    {
        res.push_back(it_sv1.ind(), it_sv1.el());
        ++it_sv1;
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.is_ordered(true);
#endif
}

template <typename T, typename A>
void mul(A &res, const A &sv, const T &val)
{
    res.n = sv.n;
    res.clear();
    if(!res.dynamic())
    {
        res.reserve_mem(sv.real_size());
    }
    A::const_indraw_iterator it_sv     = sv.begin();
    A::const_indraw_iterator it_sv_end = sv.end();
    for (;it_sv != it_sv_end;++it_sv)
    {
        res.push_back(it_sv.ind(), it_sv.el()* val);
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    res.is_ordered(sv.is_ordered());
#endif
}

template <typename T, typename A>
void div(A &res, const A &sv, const T &val)
{
    res.n = sv.n;
    res.clear();
    if(!res.dynamic())
    {
        res.reserve_mem(sv.real_size());
    }
    A::const_indraw_iterator it_sv     = sv.begin();
    A::const_indraw_iterator it_sv_end = sv.end();
    for (;it_sv != it_sv_end;++it_sv)
    {
        res.push_back(it_sv.ind(), it_sv.el()/ val);
    }


#ifdef ARAGELI_DEBUG_LEVEL_3
    res.is_ordered(sv.is_ordered());
#endif
}

template <typename T, typename A, typename B>
void dot(T &res, const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    res = null<T>();

    A::const_indraw_iterator it_sv1     = sv1.begin();
    A::const_indraw_iterator it_sv1_end = sv1.end();
    B::const_indraw_iterator it_sv2     = sv2.begin();
    B::const_indraw_iterator it_sv2_end = sv2.end();
    for (;it_sv1 != it_sv1_end && it_sv2 != it_sv2_end;)
    {
        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            res += it_sv1.el()* it_sv2.el();
            ++it_sv1;
            ++it_sv2;
        }
        else if(it_sv1.ind()< it_sv2.ind())
        {
            ++it_sv1;
        }
        else
        {
            ++it_sv2;
        }
    }
}

template <typename T, typename A>
void dot(T &res, const A &sv, const Arageli::vector<T> &vec)
{
    ARAGELI_ASSERT_0( sv.size()==vec.size());

    res = null<T>();

    A::const_indraw_iterator it_sv     = sv.begin();
    A::const_indraw_iterator it_sv_end = sv.end();
    for (;it_sv != it_sv_end;++it_sv)
    {
        res += it_sv.el()* vec[it_sv.ind()];
    }
}

template <typename A, typename B>
void add_in(A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    A res;
    res.n=sv1.n;
    res.clear();
    if(!sv1.dynamic())
    {
        //symbolic part
        res.reserve_mem(struct_union(sv1, sv2));
    }
    A::const_indraw_iterator it_sv1     = sv1.begin();
    A::const_indraw_iterator it_sv1_end = sv1.end();
    B::const_indraw_iterator it_sv2     = sv2.begin();
    B::const_indraw_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            res.push_back(it_sv1.ind(),it_sv1.el());
            ++it_sv1;
        }

        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            res.push_back(it_sv1.ind(),it_sv1.el()+ it_sv2.el());
            ++it_sv1;
        }
        else
        {
            res.push_back(it_sv2.ind(),it_sv2.el());
        }
    }

    while (it_sv1 != it_sv1_end)
    {
        res.push_back(it_sv1.ind(),it_sv1.el());
        ++it_sv1;
    }

    sv1.rep.swap(res.rep);
}

template <typename A, typename B>
void sub_in(A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    //symbolic part*
    A res;
    res.n=sv1.n;
    res.clear();
    if(!sv1.dynamic())
    {
        res.reserve_mem(struct_union(sv1, sv2));
    }
    A::const_indraw_iterator it_sv1     = sv1.begin();
    A::const_indraw_iterator it_sv1_end = sv1.end();
    B::const_indraw_iterator it_sv2     = sv2.begin();
    B::const_indraw_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            res.push_back(it_sv1.ind(),it_sv1.el());
            ++it_sv1;
        }

        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            res.push_back(it_sv1.ind(),it_sv1.el()- it_sv2.el());
            ++it_sv1;
        }
        else
        {
            res.push_back(it_sv2.ind(),-it_sv2.el());
        }
    }

    while (it_sv1 != it_sv1_end)
    {
        res.push_back(it_sv1.ind(),it_sv1.el());
        ++it_sv1;
    }

    sv1.rep.swap(res.rep);

}

template <typename T, typename A>
void mul_in(A &sv, const T &val)
{
    A::raw_iterator it_sv     = sv.begin();
    A::raw_iterator it_sv_end = sv.end();
    for (;it_sv != it_sv_end;++it_sv)
    {
        it_sv.el()*= val;
    }
}

template <typename T, typename A>
void div_in(A &sv, const T &val)
{
    A::raw_iterator it_sv     = sv.begin();
    A::raw_iterator it_sv_end = sv.end();
    for (;it_sv != it_sv_end;++it_sv)
    {
        it_sv.el()/= val;
    }
}

template <typename A, typename B>
size_t struct_intersection(const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    size_t count = 0;
    A::const_ind_iterator it_sv1     = sv1.begin();
    A::const_ind_iterator it_sv1_end = sv1.end();
    B::const_ind_iterator it_sv2     = sv2.begin();
    B::const_ind_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            ++it_sv1;
        }
        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            ++count;
            ++it_sv1;
        }
    }
    return count;
}

template <typename A, typename B>
size_t struct_union(const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    size_t count = 0;
    A::const_ind_iterator it_sv1     = sv1.begin();
    A::const_ind_iterator it_sv1_end = sv1.end();
    B::const_ind_iterator it_sv2     = sv2.begin();
    B::const_ind_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            ++count;
            ++it_sv1;
        }

        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            ++count;
            ++it_sv1;
        }
        else
        {
            ++count;
        }
    }

    count += it_sv1.distance(it_sv1_end);

    return count;
}

template <typename A, typename B>
bool struct_disjoint(const A &sv1, const B &sv2)
{
    ARAGELI_ASSERT_0( sv1.size()==sv2.size());
    ARAGELI_ASSERT_2( sv1.is_ordered()&& sv2.is_ordered());

    A::const_ind_iterator it_sv1     = sv1.begin();
    A::const_ind_iterator it_sv1_end = sv1.end();
    B::const_ind_iterator it_sv2     = sv2.begin();
    B::const_ind_iterator it_sv2_end = sv2.end();
    for (;it_sv2 != it_sv2_end;++it_sv2)
    {
        while (it_sv1 != it_sv1_end && it_sv1.ind()< it_sv2.ind())
        {
            ++it_sv1;
        }
        if (it_sv1 != it_sv1_end && it_sv1.ind()== it_sv2.ind())
        {
            return false;
        }
    }



    return true;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}// namespace sprsvc_comp

// Sparse vector function implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename T,typename Param>
sparse_vector<T,Param>::sparse_vector(size_t n=0):
     rep(n)
{
}

template<typename T,typename Param>
template<typename inT, typename inParam>
sparse_vector<T,Param>::sparse_vector(const sparse_vector<inT,inParam> &sp)
{
    sprsvc_comp::copy_to(rep,sp.rep);
}

template<typename T,typename Param>
sparse_vector<T,Param>::sparse_vector(const Arageli::vector<T> &vec)
{
    sprsvc_comp::copy_to(rep,vec);
}

template<typename T,typename Param>
sparse_vector<T,Param>
    ::~sparse_vector()
{
}

template<typename T,typename Param>
void sparse_vector<T,Param>::resize(size_t n)
{
    rep.resize(n);
}

template<typename T,typename Param>
void sparse_vector<T,Param>::reserve_mem(size_t i)
{
    rep.reserve_mem(i);
}

template<typename T,typename Param>
void sparse_vector<T,Param>::resize_mem_struct(size_t elem_count)
{
    rep.resize_mem_struct(elem_count);
}

template<typename T,typename Param>
const T & sparse_vector<T,Param>::operator[](size_t i) const
{
    return rep.get(i);
}

template<typename T,typename Param>
template<typename T1, typename Param1>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator= (const sparse_vector<T1,Param1> &sv)
{
    sprsvc_comp::copy_to(rep, sv.rep);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator= (const sparse_vector &sv)
{
    sprsvc_comp::copy_to(rep, sv.rep);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator= (const Arageli::vector<T> &vec)
{
    sprsvc_comp::copy_to(rep, vec);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator= (const char *ch)
{
    std::istringstream buf(ch);
    // WARNING. It is valid if there are no virtual function.
    static_cast<std::istream&>(buf)>> *this;
    if(!buf && !buf.eof())
        throw incorrect_string(ch);
    return *this;
}

template<typename T,typename Param>
template<typename inT, typename inParam>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator+=(const sparse_vector<inT,inParam> &sv)
{
    sprsvc_comp::add_in(rep, sv.rep);
    return *this;
}

template<typename T,typename Param>
template<typename inT, typename inParam>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator-=(const sparse_vector<inT,inParam> &sv)
{
    sprsvc_comp::sub_in(rep, sv.rep);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator*=(const T &val)
{
    sprsvc_comp::mul_in(rep, val);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param> & sparse_vector<T,Param>::operator/=(const T &val)
{
    sprsvc_comp::div_in(rep, val);
    return *this;
}

template<typename T,typename Param>
sparse_vector<T,Param>::operator vector<T> ()
{
    vector<T> dst;
    copy_to(dst, *this);
    return dst;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool operator==(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    return sprsvc_comp::equal(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool operator!=(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    return sprsvc_comp::non_equal(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_vector<T1,Param1> operator+ (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sparse_vector<T1,Param1> res;
    sprsvc_comp::add(res.rep,sv1.rep,sv2.rep);
    return res;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_vector<T1,Param1> operator- (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sparse_vector<T1,Param1> res;
    sprsvc_comp::sub(res.rep,sv1.rep,sv2.rep);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const big_int &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::mul(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const big_int &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::div(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const float &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::mul(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const float &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::div(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const double &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::mul(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const double &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::div(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const rational<big_int> &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::mul(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const rational<big_int> &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::div(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const big_float &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::mul(res.rep,sv.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const big_float &val)
{
    sparse_vector<T,Param> res;
    sprsvc_comp::div(res.rep,sv.rep,val);
    return res;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<typename T1, typename Param1, typename T2, typename Param2>
inline T1 operator* (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    T1 res;
    sprsvc_comp::dot(res,sv1.rep,sv2.rep);
    return res;
}

template<typename T, typename T1, typename Param1>
inline T1 operator* (const sparse_vector<T1,Param1> &sv, const Arageli::vector<T> &vec)
{
    T1 res;
    sprsvc_comp::dot(res,sv.rep,vec);
    return res;
}

template<typename T,typename Param>
std::istream & operator >> ( std::istream & s, sparse_vector<T,Param> & x)
{
    return input_list(s,x);
}

template<typename T,typename Param>
std::ostream & operator << ( std::ostream & s, const sparse_vector<T,Param> & x)
{
    ARAGELI_ASSERT_2(x.is_ordered());

    return output_list(s,x);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// IO  sparse sparse vector operations
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Simple output of the sparse vector
template <typename Out, typename T, typename Param>
Out& output_list1
(
    Out& out,
    const sparse_vector<T, Param>& x,
    const char* first_bracket,
    const char* second_bracket,
    const char* separator,
    const char* first_pair_bracket,
    const char* second_pair_bracket,
    const char* second_pair_separator
)
{
    out << first_bracket;
    sparse_vector<T,Param>::const_indraw_iterator
        it     = x.begin(),
        it_end = x.end();
    while(it!=it_end)
    {
        out << first_pair_bracket;
        out << it.ind();
        out << second_pair_separator;
        out << it.el();
        out << second_pair_bracket;
        ++it;
        if(it!=it_end)
        {
            out << separator;
        }
    }
    out << second_bracket;

    return out;
}
// Simple output of the sparse vector
template <typename Out, typename T, typename Param>
Out& output_list
(
    Out& out,
    const sparse_vector<T, Param>& x,
    const char* first_bracket,
    const char* second_bracket,
    const char* separator,
    const char* null_symbol,
    const char* index_order_error
)
{
    out << first_bracket;
    size_t i=0;
    sparse_vector<T,Param>::const_indraw_iterator
        it     = x.begin(),
        it_end = x.end();
    while(it!=it_end && it.ind()<x.size())
    {
        if(i<it.ind())
        {
            out << null_symbol;
        }
        else if(i==it.ind())
        {
            out << it.el();
            ++it;
        }
        else
        {
            out << index_order_error;
            out << second_bracket;
            return out;
        }

        ++i;
        if(i<x.size())
        {
            out << separator;
        }
        ARAGELI_ASSERT_1(i<=x.size());
    }
    while(i++ < x.size())
    {
        out << null_symbol;
        if(i<x.size())
        {
            out << separator;
        }
    }

    out << second_bracket;

    return out;
}
/// Simple input of the sparse vector
template <typename In, typename T, typename Param>
In& input_list
(
    In& in,
    sparse_vector<T, Param>& x,
    const char* first_bracket,
    const char* second_bracket,
    const char* separator,
    const char* range
)
{
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(first_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(second_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(separator));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(range));

    if(_Internal::is_bad_read_literal(in, first_bracket))
        return in;

    // (second_bracket == "")is special case
    // in this case vector must be terminated by the some invalid character

    if(*second_bracket && _Internal::read_literal(in, second_bracket))
    {
        // empty vector
        x.clear();
        return in;
    }

    std::list<size_t> buf_ind;// a temporary buffer for index
    std::list<T> buf_val;// a temporary buffer for values

    size_t i=0;
    T t;

    _Internal::turn_off_comma_as_separator<In> _tocas(in);
    if
    (
        separator && separator[0] == ',' ||    // WARNING! Explicit ','
        second_bracket && second_bracket[0] == ','    // WARNING! Explicit ','
    )
        _tocas.activate();

    for(;;)
    {
        in >> t;

        if(_Internal::is_bad_input(in))
            return in;
        else if(_Internal::read_literal(in, separator))
        {// "t,"
            if(!is_null(t))
            {
                buf_ind.push_back(i);
                buf_val.push_back(t);
            }
            ++i;
            continue;
        }
        else if(_Internal::read_literal(in, range))
        {// "t:"
            T t1;
            in >> t1;

            if(_Internal::is_bad_input(in))
                return in;
            else if(_Internal::read_literal(in, separator))
            {// "t:t1,"
                generate_range_helper(t, t1, std::back_inserter(buf_val));
                continue;
            }
            else if(_Internal::read_literal(in, range))
            {// "t:t1:"
                T t2;
                in >> t2;

                if(_Internal::is_bad_input(in))
                    return in;
                else if(_Internal::read_literal(in, separator))
                {// "t:t1:t2,"
                    generate_range_helper(t, t1, t2, std::back_inserter(buf_val));
                    continue;
                }
                else if
                (
                    *second_bracket == 0 ||
                    _Internal::read_literal(in, second_bracket)
                )
                {// "t:t1:t2)"
                    generate_range_helper(t, t1, t2, std::back_inserter(buf_val));
                    break;
                }
            }
            else if(*second_bracket == 0 || _Internal::read_literal(in, second_bracket))
            {// "t:t1)"
                generate_range_helper(t, t1, std::back_inserter(buf_val));
                break;
            }
        }
        else if(*second_bracket == 0 || _Internal::read_literal(in, second_bracket))
        {// "t)"
            if(!is_null(t))
            {
                buf_ind.push_back(i);
                buf_val.push_back(t);
            }
            ++i;
            break;
        }

        in.clear(std::ios_base::badbit);
        return in;
    }

    x.resize(i);
    if(!x.dynamic())
    {
        x.reserve_mem(buf_ind.size());
    }
    std::list<size_t>::iterator it_ind = buf_ind.begin();
    std::list<T>::iterator it_val = buf_val.begin();
    std::list<size_t>::iterator it_ind_end = buf_ind.end();
    std::list<T>::iterator it_val_end = buf_val.end();

    for(;it_ind!=it_ind_end;++it_ind,++it_val)
    {
        x.push_back(*it_ind,*it_val);
    }
    return in;
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, typename Param>
inline bool is_null(const sparse_vector<T, Param>& x)
{
    sparse_vector<T, Param>::const_raw_iterator it;
    sparse_vector<T, Param>::const_raw_iterator it_end;
    it = x.begin();
    it_end = x.end();
    for(;it!=it_end;++it)
    {
        if(!is_null(it.el()))
        {
            return false;
        }
    }

    return true;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline void copy_to(sparse_vector<T1,Param1> &dst, const sparse_vector<T2,Param2> &src)
{
    ARAGELI_ASSERT_0((void*)&dst!=(void*)&src);
    sprsvc_comp::copy_to(dst.rep,src.rep);
}

template <typename T,typename Param1, typename T1>
inline void copy_to(sparse_vector<T1,Param1> &dst, const Arageli::vector<T> &src)
{
    sprsvc_comp::copy_to(dst.rep,src);
}

template <typename T,typename Param1, typename T1>
inline void copy_to(Arageli::vector<T> &dst, const sparse_vector<T1,Param1> &src)
{
    sprsvc_comp::copy_to(dst,src.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline  bool equal(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sprsvc_comp::equal(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool non_equal(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sprsvc_comp::non_equal(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline void add(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    ARAGELI_ASSERT_0((void*)&res!=(void*)&sv1 && (void*)&res!=(void*)&sv2);
    sprsvc_comp::add(res.rep,sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline void sub(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    ARAGELI_ASSERT_0((void*)&res!=(void*)&sv1 && (void*)&res!=(void*)&sv2);
    sprsvc_comp::sub(res.rep,sv1.rep,sv2.rep);
}

template <typename T, typename T1,typename Param1>
inline void mul(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv, const T &val)
{
    ARAGELI_ASSERT_0((void*)&res!=(void*)&sv);
    sprsvc_comp::mul(res.rep,sv.rep,val);
}

template<typename T, typename T1,typename Param1>
inline void div(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv, const T &val)
{
    ARAGELI_ASSERT_0((void*)&res!=(void*)&sv);
    sprsvc_comp::div(res.rep,sv.rep,val);
}

template<typename T, typename T1, typename Param1, typename T2, typename Param2>
inline void dot(T &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sprsvc_comp::dot(res,sv1.rep,sv2.rep);
}

template<typename T, typename Param, typename T1,  typename T2>
inline void dot(T1 &res, const sparse_vector<T,Param> &sv, const Arageli::vector<T2> &vec)
{
    sprsvc_comp::dot(res,sv.rep,vec);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline void add_in(sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sprsvc_comp::add_in(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline void sub_in(sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    sprsvc_comp::sub_in(sv1.rep,sv2.rep);
}

template<typename T, typename T1, typename Param1>
inline void mul_in(sparse_vector<T1,Param1> &sv, const T &val)
{
    sprsvc_comp::mul_in(sv.rep,val);
}

template<typename T, typename T1, typename Param1>
inline void div_in(sparse_vector<T1,Param1> &sv, const T &val)
{
    sprsvc_comp::div_in(sv.rep,val);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline size_t struct_intersection(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    return sprsvc_comp::struct_intersection(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline size_t struct_union(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    return sprsvc_comp::struct_union(sv1.rep,sv2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool struct_disjoint(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2)
{
    return sprsvc_comp::struct_disjoint(sv1.rep,sv2.rep);
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}// namespace Arageli
#else    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)|| ...

#include "sparse_vector.hpp"

namespace Arageli
{

/// IO sparse vector params
const char* sparse_vector_output_list_first_bracket_default = "(";
const char* sparse_vector_output_list_second_bracket_default = ")";
const char* sparse_vector_output_list_separator_default = ", ";
const char* sparse_vector_output_list_null_symbol_default = "0";
const char* sparse_vector_output_list_first_pair_bracket_default = "[";
const char* sparse_vector_output_list_second_pair_bracket_default = "]";
const char* sparse_vector_output_list_pair_separator_default = ", ";
const char* sparse_vector_output_standart_order_error_default = "order index error";
const char* sparse_vector_input_list_first_bracket_default = "(";
const char* sparse_vector_input_list_second_bracket_default = ")";
const char* sparse_vector_input_list_separator_default = ",";
const char* sparse_vector_input_list_range_default = ":";

}// namespace Arageli


#endif    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)|| ...
