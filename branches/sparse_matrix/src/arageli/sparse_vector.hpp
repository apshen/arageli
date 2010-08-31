/*****************************************************************************

    sparse_vector.hpp

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
    \file sparse_vector.hpp
    \brief <!--ADD BRIEF HEADER DESCRIPTION HERE-->

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/
#ifndef _ARAGELI_sparse_vector_hpp_
#define _ARAGELI_sparse_vector_hpp_

#define ARAGELI_DEBUG_LEVEL_3 ARAGELI_DEBUG_LEVEL>=3


#include "config.hpp"
#include <map>

#include <cstddef>
#include <utility>

#include "type_traits.hpp"
#include "vector.hpp"

// REFERENCE ADDITIONAL HEADERS HERE


namespace Arageli
{

struct mappair_t {};

ARAGELI_REGISTER_TAG(mappair_t)

struct vecpair_t {};

ARAGELI_REGISTER_TAG(vecpair_t)

struct pairvec_t {};

ARAGELI_REGISTER_TAG(pairvec_t)

struct listpair_t {};

ARAGELI_REGISTER_TAG(listpair_t)

struct deqpair_t {};

ARAGELI_REGISTER_TAG(deqpair_t)

namespace sprsvc_comp
{
template <typename A, typename B>
void copy_to(A &dst, const B &src);

template <typename T,typename A, typename B>
void copy_to(A &dst, const Arageli::vector<T> &src);

template <typename T>
bool is_notvalid(const T &value);
}

namespace sprsvc_rep
{

namespace iterator_base
{

template <typename TypeIterator, typename ValType>
struct indraw_iterator_base
{
    // pointer
    TypeIterator it;
public:
    indraw_iterator_base()
    {
    }

    indraw_iterator_base(const TypeIterator &i) :
        it(i)
    {
    }

    template <typename TypeIterator, typename ValType>
    indraw_iterator_base(const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename TypeIterator, typename ValType>
    indraw_iterator_base & operator= (const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        it=i.it;
        return *this;
    }

    bool operator==(const indraw_iterator_base &i)
    {
        return it==i.it;
    }

    bool operator!=(const indraw_iterator_base &i)
    {
        return it!=i.it;
    }

    indraw_iterator_base & operator++()
    {
        ++it;
        return *this;
    }

    indraw_iterator_base   operator++(int)
    {
        indraw_iterator_base tmp(it);
        ++it;
        return tmp;
    }

    indraw_iterator_base & operator--()
    {
        --it;
        return *this;
    }
    indraw_iterator_base   operator--(int)
    {
        indraw_iterator_base tmp(it);
        --it;
        return tmp;
    }

    const size_t ind() const
    {
        return it->first;
    }

    inline ValType & el() const
    {
        return it->second;
    }

    void advance(int k)
    {
        std::advance(it,k);
    }

    template <typename TypeIterator, typename ValType>
    int distance(indraw_iterator_base<TypeIterator,ValType> &i)
    {
        return std::distance(it,i.it);
    }
};

template <typename TypeIterator, typename ValType>
struct raw_iterator_base
{
    // pointer
    TypeIterator it;
public:
    raw_iterator_base()
    {
    }

    raw_iterator_base(const TypeIterator &i) :
        it(i)
    {
    }

    template <typename TypeIterator, typename ValType>
    raw_iterator_base(const raw_iterator_base<TypeIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename TypeIterator, typename ValType>
    raw_iterator_base(const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename TypeIterator, typename ValType>
    raw_iterator_base & operator= (const raw_iterator_base<TypeIterator,ValType> &i)
    {
        it=i.it;
        return *this;
    }

    template <typename TypeIterator, typename ValType>
    raw_iterator_base & operator= (const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        it=i.it;
        return *this;
    }

    bool operator==(const raw_iterator_base &i)
    {
        return it==i.it;
    }

    bool operator!=(const raw_iterator_base &i)
    {
        return it!=i.it;
    }

    raw_iterator_base & operator++()
    {
        ++it;
        return *this;
    }

    raw_iterator_base   operator++(int)
    {
        raw_iterator_base tmp(it);
        ++it;
        return tmp;
    }

    raw_iterator_base & operator--()
    {
        --it;
        return *this;
    }
    raw_iterator_base   operator--(int)
    {
        raw_iterator_base tmp(it);
        --it;
        return tmp;
    }

    ValType & el() const
    {
        return it->second;
    }

    void advance(int k)
    {
        std::advance(it,k);
    }

    template <typename TypeIterator, typename ValType>
    int distance(raw_iterator_base<TypeIterator,ValType> &i)
    {
        return std::distance(it,i.it);
    }
};

template <typename TypeIterator>
struct ind_iterator_base
{
    // pointer
    TypeIterator it;
public:
    ind_iterator_base()
    {
    }

    ind_iterator_base(const TypeIterator &i) :
        it(i)
    {
    }

    template <typename TypeIterator, typename ValType>
    ind_iterator_base(const ind_iterator_base<TypeIterator> &i)
    {
        (*this)= i;
    }

    template <typename TypeIterator, typename ValType>
    ind_iterator_base(const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename TypeIterator, typename ValType>
    ind_iterator_base & operator= (const ind_iterator_base<TypeIterator> &i)
    {
        it=i.it;
        return *this;
    }

    template <typename TypeIterator, typename ValType>
    ind_iterator_base & operator= (const indraw_iterator_base<TypeIterator,ValType> &i)
    {
        it=i.it;
        return *this;
    }

    bool operator==(const ind_iterator_base &i)
    {
        return it==i.it;
    }

    bool operator!=(const ind_iterator_base &i)
    {
        return it!=i.it;
    }

    ind_iterator_base & operator++()
    {
        ++it;
        return *this;
    }

    ind_iterator_base   operator++(int)
    {
        ind_iterator_base tmp(it);
        ++it;
        return tmp;
    }

    ind_iterator_base & operator--()
    {
        --it;
        return *this;
    }
    ind_iterator_base   operator--(int)
    {
        ind_iterator_base tmp(it);
        --it;
        return tmp;
    }

    const size_t ind() const
    {
        return it->first;
    }

    void advance(int k)
    {
        std::advance(it,k);
    }

    template <typename TypeIterator>
    int distance(ind_iterator_base<TypeIterator> &i)
    {
        return std::distance(it,i.it);
    }
};

template <typename IndIterator, typename ElemIterator, typename ValType>
struct indraw_complex_iterator_base
{
    // pointer
    IndIterator in_it;
    ElemIterator el_it;

    indraw_complex_iterator_base()
    {
    }

    indraw_complex_iterator_base(const IndIterator &in, ElemIterator &el) :
        in_it(in), el_it(el)
    {
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    indraw_complex_iterator_base(const indraw_complex_iterator_base<IndIterator,ElemIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    indraw_complex_iterator_base & operator= (const indraw_complex_iterator_base<IndIterator,ElemIterator,ValType> &i)
    {
        in_it=i.in_it;
        el_it=i.el_it;
        return *this;
    }

    bool operator==(const indraw_complex_iterator_base &i)
    {
        return in_it==i.in_it && el_it==i.el_it;
    }
    bool operator!=(const indraw_complex_iterator_base &i)
    {
        return in_it!=i.in_it || el_it!=i.el_it;
    }
    indraw_complex_iterator_base & operator++()
    {
        ++in_it;
        ++el_it;
        return *this;
    }

    indraw_complex_iterator_base   operator++(int)
    {
        indraw_complex_iterator_base tmp(in_it,el_it);
        ++(*this);
        return tmp;
    }

    indraw_complex_iterator_base & operator--()
    {
        --in_it;
        --el_it;
        return *this;
    }

    indraw_complex_iterator_base   operator--(int)
    {
        indraw_complex_iterator_base tmp(in_it,el_it);
        --(*this);
        return tmp;
    }

    const size_t ind() const
    {
        return *in_it;
    }

    ValType & el() const
    {
        return *el_it;
    }

    void advance(int k)
    {
        std::advance(in_it,k);
        std::advance(el_it,k);
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    int distance(indraw_complex_iterator_base<IndIterator,ElemIterator,ValType> &i)
    {
        return std::distance(in_it,i.in_it);
    }
};

template <typename ElemIterator, typename ValType>
struct raw_single_iterator_base
{
    // pointer
    ElemIterator el_it;

    raw_single_iterator_base()
    {
    }

    raw_single_iterator_base(const ElemIterator &el) :
        el_it(el)
    {
    }

    template <typename ElemIterator, typename ValType>
    raw_single_iterator_base(const raw_single_iterator_base<ElemIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename ElemIterator, typename ValType>
    raw_single_iterator_base & operator= (const raw_single_iterator_base<ElemIterator,ValType> &i)
    {
        el_it=i.el_it;
        return *this;
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    raw_single_iterator_base(const indraw_complex_iterator_base<IndIterator, ElemIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    raw_single_iterator_base & operator= (const indraw_complex_iterator_base<IndIterator, ElemIterator,ValType> &i)
    {
        el_it=i.el_it;
        return *this;
    }

    bool operator==(const raw_single_iterator_base &i)
    {
        return el_it==i.el_it;
    }

    bool operator!=(const raw_single_iterator_base &i)
    {
        return el_it!=i.el_it;
    }

    raw_single_iterator_base & operator++()
    {
        ++el_it;
        return *this;
    }

    raw_single_iterator_base   operator++(int)
    {
        raw_single_iterator_base tmp(el_it);
        ++el_it;
        return tmp;
    }

    raw_single_iterator_base & operator--()
    {
        --el_it;
        return *this;
    }

    raw_single_iterator_base   operator--(int)
    {
        raw_single_iterator_base tmp(in_it,el_it);
        --el_it;
        return tmp;
    }

    ValType & el() const
    {
        return *el_it;
    }

    void advance(int k)
    {
        std::advance(el_it,k);
    }

    template <typename ElemIterator, typename ValType>
    int distance(raw_single_iterator_base<ElemIterator,ValType> &i)
    {
        return std::distance(el_it,i.el_it);
    }
};

template <typename IndIterator>
struct ind_single_iterator_base
{
    // pointer
    IndIterator in_it;

    ind_single_iterator_base()
    {
    }

    ind_single_iterator_base(const IndIterator &in, const IndIterator &el) :
        in_it(in), in_it(el)
    {
    }

    template <typename IndIterator>
    ind_single_iterator_base(const ind_single_iterator_base<IndIterator> &i)
    {
        (*this)= i;
    }
    template <typename IndIterator>
    ind_single_iterator_base & operator= (const ind_single_iterator_base<IndIterator> &i)
    {
        in_it=i.in_it;
        return *this;
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    ind_single_iterator_base(const indraw_complex_iterator_base<IndIterator,ElemIterator,ValType> &i)
    {
        (*this)= i;
    }

    template <typename IndIterator, typename ElemIterator, typename ValType>
    ind_single_iterator_base & operator= (const indraw_complex_iterator_base<IndIterator,ElemIterator,ValType> &i)
    {
        in_it=i.in_it;
        return *this;
    }

    bool operator==(const ind_single_iterator_base &i)
    {
        return in_it==i.in_it;
    }

    bool operator!=(const ind_single_iterator_base &i)
    {
        return in_it!=i.in_it;
    }

    ind_single_iterator_base & operator++()
    {
        ++in_it;
        return *this;
    }

    ind_single_iterator_base   operator++(int)
    {
        ind_single_iterator_base tmp(in_it);
        ++in_it;
        return tmp;
    }

    ind_single_iterator_base & operator--()
    {
        --in_it;
        return *this;
    }

    ind_single_iterator_base   operator--(int)
    {
        ind_single_iterator_base tmp(in_it,in_it);
        --in_it;
        return tmp;
    }

    const size_t ind() const
    {
        return *in_it;
    }

    void advance(int k)
    {
        std::advance(in_it,k);
    }

    template <typename IndIterator>
    int distance(ind_single_iterator_base<IndIterator> &i)
    {
        return std::distance(in_it,i.in_it);
    }
};

}// namespace iterator


/// Template class representation of sparse vector definition.
/** sparse vector on vector of pairs
    @param T type of vector's items
    @param REFCNT reference counter switch */
template <typename T, bool REFCNT = true>
class sparsevec_vecpair_t
{
public:
    typedef typename    T                               value_type;
    typedef typename    std::pair<size_t,T>             Elem;
    typedef typename    Arageli::vector<Elem, REFCNT>   Rep;
    typedef typename    Rep::iterator                   Rep_iterator;
    typedef typename    Rep::const_iterator             Const_Rep_iterator;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = REFCNT;
public:
    sparsevec_vecpair_t(size_t _n=0);

    ~sparsevec_vecpair_t();

    const size_t size() const
    {
        return n;
    }

    const size_t real_size() const
    {
        return rep.size();
    }

    const size_t capacity() const
    {
        return rep.capacity();
    }

    void resize(size_t _n);

    void reserve_mem(size_t elem_count);

    void resize_mem_struct(size_t elem_count);

    void regularize()
    {
        sorting(rep, 0, real_size()-1);
        remove_useless();


#ifdef ARAGELI_DEBUG_LEVEL_3
        ordered = true;
#endif

    };

#ifdef ARAGELI_DEBUG_LEVEL_3

    bool is_ordered() const
    {
        return ordered;
    }

    void is_ordered(bool ord)
    {
        ordered = ord;
    }

    bool check_ordered() const
    {
        for(size_t i=1;i<rep.size();++i)
        {
            if(rep[i-1].first>=rep[i].first)
            {
                return false;
            }
        }
        return true;
    }

#endif

    void clear();

    void insert(size_t i, const T &t)
    {
        ARAGELI_ASSERT_0(i<n);

        if(!null(t))
        {
            ins(i,t);
        }
        else
        {
            del(i);
        }
    }

    void set_elem(size_t real_ind, size_t ind, const T & t)
    {
        ARAGELI_ASSERT_0(ind<n && real_ind<real_size()&& !sprsvc_comp::is_notvalid(t));

        raw_set_elem(real_ind,ind,t);

#ifdef ARAGELI_DEBUG_LEVEL_3
        ordered = check_ordered();
#endif

    }

    void raw_set_elem(size_t real_ind, size_t ind, const T & t);

    void push_back(size_t ind, const T & t);

    const T & get(size_t i) const;

    bool dynamic() const
    {
        return false;
    }

    void swap(sparsevec_vecpair_t<T,REFCNT> &x)
    {
        rep.swap(x.rep);
    }
protected:

    const size_t real_index(size_t k) const;

public:

    typedef iterator_base::indraw_iterator_base
    <
        Rep_iterator,
        T
    >
    indraw_iterator;

    typedef iterator_base::raw_iterator_base
    <
        Rep_iterator,
        T
    >
    raw_iterator;

    typedef iterator_base::raw_iterator_base
    <
        Const_Rep_iterator,
        const T
    >
    const_raw_iterator;

    typedef iterator_base::ind_iterator_base
    <
        Const_Rep_iterator
    >
    const_ind_iterator;

    typedef iterator_base::indraw_iterator_base
    <
        Const_Rep_iterator,
        const T
    >
    const_indraw_iterator;

    indraw_iterator begin()
    {
        return indraw_iterator(rep.begin());
    }

    indraw_iterator end()
    {
        return indraw_iterator(rep.end());
    }

    const_indraw_iterator begin() const
    {
        return const_indraw_iterator(rep.begin());
    }

    const_indraw_iterator end() const
    {
        return const_indraw_iterator(rep.end());
    }


public:

    indraw_iterator find(size_t i)
    {
        indraw_iterator it = end();
        if(real_size())
        {
            size_t k = real_index(i);
            if(k<real_size())
            {
                if(rep[k].first == i)
                {
                    it.advance((int)k-(int)real_size());
                }
            }
        }

        return it;
    }

    const_indraw_iterator find(size_t i) const
    {
        const_indraw_iterator it = end();
        if(real_size())
        {
            size_t k = real_index(i);
            if(k<real_size())
            {
                if(rep[k].first == i)
                {
                    it.advance((int)k-(int)real_size());
                }
            }
        }

        return it;
    }
public:

    Rep    rep;

#ifdef ARAGELI_DEBUG_LEVEL_3
    bool ordered;
#endif

    size_t n;// vector dimension

private:

    void ins(size_t i, const T &t);

    void del(size_t i);

    void sorting(Rep &A, int first, int last);

    void remove_useless();
};

/// Template class representation of sparse vector definition.
/** sparse vector on pair of vectors
    @param T type of vector's items
    @param REFCNT reference counter switch */
template <typename T, bool REFCNT = true>
class sparsevec_pairvec_t
{
public:
    typedef typename    T                               value_type;
    typedef typename    Arageli::vector<size_t, REFCNT> Indeces;
    typedef typename    Arageli::vector<T, REFCNT>      Elems;
    typedef typename    std::pair<Indeces,Elems>        Rep;
    typedef typename    Indeces::iterator               In_Rep_iterator;
    typedef typename    Elems::iterator                 El_Rep_iterator;
    typedef typename    Indeces::const_iterator         Const_In_Rep_iterator;
    typedef typename    Elems::const_iterator           Const_El_Rep_iterator;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = REFCNT;

    sparsevec_pairvec_t(size_t _n=0);
    ~sparsevec_pairvec_t();

    const size_t size() const
    {
        return n;
    }

    const size_t real_size() const
    {
        return rep.first.size();
    }

    const size_t capacity() const
    {
        return rep.first.capacity();
    }

    void resize(size_t _n);

    void reserve_mem(size_t elem_count);

    void resize_mem_struct(size_t elem_count);

    void regularize()
    {

        sorting(rep, 0, real_size()-1);
        remove_useless();

#ifdef ARAGELI_DEBUG_LEVEL_3
        ordered = true;
#endif

    };

#ifdef ARAGELI_DEBUG_LEVEL_3
    bool is_ordered() const
    {
        return ordered;
    }

    void is_ordered(bool ord)
    {
        ordered = ord;
    }

    bool check_ordered() const
    {
        for(size_t i=1;i<rep.first.size();++i)
        {
            if(rep.first[i-1]>=rep.first[i])
            {
                return false;
            }
        }
        return true;
    }
#endif

    void clear();

    void insert(size_t i, const T &t)
    {
        ARAGELI_ASSERT_0(i<n);

        if(!sprsvc_comp::is_notvalid(t))
        {
            ins(i,t);
        }
        else
        {
            del(i);
        }
    }

    void set_elem(size_t real_ind, size_t ind, const T & t)
    {
        ARAGELI_ASSERT_0(ind<n && real_ind<real_size()&& !sprsvc_comp::is_notvalid(t));

        raw_set_elem(real_ind,ind,t);

#ifdef ARAGELI_DEBUG_LEVEL_3
        ordered = check_ordered();
#endif

    }

    void raw_set_elem(size_t real_ind, size_t ind, const T & t);

    void push_back(size_t ind, const T & t);

    const T & get(size_t i) const;

    bool dynamic() const
    {
        return false;
    }

    void swap(sparsevec_pairvec_t<T,REFCNT> &x)
    {
        rep.swap(x.rep);
    }

protected:

    const size_t real_index(size_t k) const;

public:

    typedef iterator_base::indraw_complex_iterator_base
    <
        Const_In_Rep_iterator,
        El_Rep_iterator,
        T
    >
    indraw_iterator;

    typedef iterator_base::raw_single_iterator_base
    <
        El_Rep_iterator,
        T
    >
    raw_iterator;

    typedef iterator_base::raw_single_iterator_base
    <
        Const_El_Rep_iterator,
        const T
    >
    const_raw_iterator;

    typedef iterator_base::ind_single_iterator_base
    <
        Const_In_Rep_iterator
    >
    const_ind_iterator;

    typedef iterator_base::indraw_complex_iterator_base
    <
        Const_In_Rep_iterator,
        Const_El_Rep_iterator,
        const T
    >
    const_indraw_iterator;

    indraw_iterator begin()
    {
        return indraw_iterator(rep.first.begin(),rep.second.begin());
    }

    indraw_iterator end()
    {
        return indraw_iterator(rep.first.end(),rep.second.end());
    }

    const_indraw_iterator begin() const
    {
        return const_indraw_iterator(rep.first.begin(),rep.second.begin());
    }

    const_indraw_iterator end() const
    {
        return const_indraw_iterator(rep.first.end(),rep.second.end());
    }

public:

    indraw_iterator find(size_t i)
    {
        indraw_iterator it = end();

        if(real_size())
        {
            size_t k = real_index(i);
            if(k<real_size())
            {
                if(rep.first[k] == i)
                {
                    it.advance((int)k-(int)real_size());
                }
            }
        }

        return it;
    }

    const_indraw_iterator find(size_t i) const
    {
        const_indraw_iterator it = end();
        if(real_size())
        {
            size_t k = real_index(i);
            if(k<real_size())
            {
                if(rep.first[k] == i)
                {
                    it.advance((int)k-(int)real_size());
                }
            }
        }

        return it;
    }
public:
    Rep     rep;

#ifdef ARAGELI_DEBUG_LEVEL_3
    bool ordered;
#endif

    size_t  n;// vector dimension
private:

    void ins(size_t i, const T &t);

    void del(size_t i);

    void sorting(Rep &A, int first, int last);

    void remove_useless();
};

/// Template class representation of sparse vector definition.
/** sparse vector on map of pairs
    @param T type of vector's items
    @param REFCNT reference counter switch */
template <typename T, bool REFCNT = true>
class sparsevec_mappair_t
{
public:
    typedef typename    T                       value_type;
    typedef typename    std::pair<size_t, T>    Elem;
    typedef typename    std::map<size_t, T>     Rep;
    typedef typename    Rep::iterator           Rep_iterator;
    typedef typename    Rep::const_iterator     Const_Rep_iterator;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = REFCNT;

    sparsevec_mappair_t(size_t _n=0);

    ~sparsevec_mappair_t();

    const size_t size() const
    {
        return n;
    }

    const size_t real_size() const
    {
        return rep.size();
    }

    const size_t capacity() const
    {
        return rep.size();
    }

    void resize(size_t _n);

    void reserve_mem(size_t elem_count);

    void resize_mem_struct(size_t elem_count)
    {
    }

    void regularize()
    {
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    bool is_ordered() const
    {
        return true;
    }

    void is_ordered(bool _ordered)
    {
    }

    bool check_ordered() const
    {
        return true;
    }
#endif

    void clear();

    void insert(size_t i, const T &t)
    {
        ARAGELI_ASSERT_0(i<n);
        if(!sprsvc_comp::is_notvalid(t))
        {
            ins(i, t);
        }
        // if t==0 that del
    }

    void set_elem(size_t real_ind, size_t ind, const T & t)
    {
        ARAGELI_ASSERT_0(ind<n && !sprsvc_comp::is_notvalid(t));

        raw_set_elem(real_ind,ind,t);
    }

    void push_back(size_t ind, const T & t);

    void raw_set_elem(size_t real_ind, size_t ind, const T & t)
    {
        insert(ind, t);
    }

    const T & get(size_t i) const;

    bool dynamic() const
    {
        return true;
    }

    void swap(sparsevec_mappair_t<T,REFCNT> &x)
    {
        rep.swap(x.rep);
    }

    typedef iterator_base::indraw_iterator_base
    <
        Rep_iterator,
        T
    >
    indraw_iterator;

    typedef iterator_base::indraw_iterator_base
    <
        Const_Rep_iterator,
        const T
    >
    const_indraw_iterator;

    typedef iterator_base::raw_iterator_base
    <
        Rep_iterator,
        T
    >
    raw_iterator;

    typedef iterator_base::raw_iterator_base
    <
        Const_Rep_iterator,
        const T
    >
    const_raw_iterator;


    typedef iterator_base::ind_iterator_base
    <
        Const_Rep_iterator
    >
    const_ind_iterator;

    indraw_iterator begin()
    {
        return indraw_iterator(rep.begin());
    }

    indraw_iterator end()
    {
        return indraw_iterator(rep.end());
    }

    const_indraw_iterator begin() const
    {
        return const_indraw_iterator(rep.begin());
    }

    const_indraw_iterator end() const
    {
        return const_indraw_iterator(rep.end());
    }

    indraw_iterator find(size_t i)
    {
        return indraw_iterator(rep.find(i));
    }

    const_indraw_iterator find(size_t i) const
    {
        return const_indraw_iterator(rep.find(i));
    }

public:
    Rep rep;
    size_t  n;// vector dimension
protected:

    //Rep_iterator last;

    void ins(size_t i, const T &t);

    void del(size_t i);
};

}// namespace sprsvc_rep

/** Namespace of sparse vector representation computation function
    and sparse vector representation operation with memory function */
namespace sprsvc_comp
{
using namespace sprsvc_rep;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Some template function for sparse vector representation on vector of pairs
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
inline void copy_to(sparsevec_vecpair_t<T,REFCNT> &dst, const sparsevec_vecpair_t<T,REFCNT> &src);

template <typename T, bool REFCNT>
inline void copy_to(sparsevec_vecpair_t<T,REFCNT> &dst, const Arageli::vector<T> &src);

//template <typename T, bool REFCNT>
//void add(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2);
//
//template <typename T, bool REFCNT>
//void sub(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2);
//
//template <typename T, bool REFCNT>
//void mul(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv, const T &val);
//
//template <typename T, bool REFCNT>
//void div(sparsevec_vecpair_t<T,REFCNT> &res, const sparsevec_vecpair_t<T,REFCNT> &sv, const T &val);
//
//template <typename T, bool REFCNT>
//void dot(T &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2);
//
//template <typename T, bool REFCNT>
//void dot(T &res, const sparsevec_vecpair_t<T,REFCNT> &sv1, const Arageli::vector<T> &vec);
//
//template <typename T, bool REFCNT>
//void add_in(sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2);
//
//template <typename T, bool REFCNT>
//void sub_in(sparsevec_vecpair_t<T,REFCNT> &sv1, const sparsevec_vecpair_t<T,REFCNT> &sv2);
//
//template <typename T, bool REFCNT>
//void mul_in(sparsevec_vecpair_t<T,REFCNT> &sv, const T &val);
//
//template <typename T, bool REFCNT>
//void div_in(sparsevec_vecpair_t<T,REFCNT> &sv, const T &val);
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Some template function for sparse vector representation on pair of vectors
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, bool REFCNT>
inline void copy_to(sparsevec_pairvec_t<T,REFCNT> &dst, const sparsevec_pairvec_t<T,REFCNT> &src);

template <typename T, bool REFCNT>
inline void copy_to(sparsevec_pairvec_t<T,REFCNT> &dst, const Arageli::vector<T> &src);

template <typename T, bool REFCNT>
inline void add(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2);

template <typename T, bool REFCNT>
inline void sub(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2);

template <typename T, bool REFCNT>
inline void mul(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv, const T &val);

template <typename T, bool REFCNT>
inline void div(sparsevec_pairvec_t<T,REFCNT> &res, const sparsevec_pairvec_t<T,REFCNT> &sv, const T &val);

template <typename T, bool REFCNT>
inline void dot(T &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2);

template <typename T, bool REFCNT>
inline void dot(T &res, const sparsevec_pairvec_t<T,REFCNT> &sv1, const Arageli::vector<T> &vec);

template <typename T, bool REFCNT>
inline void add_in(sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2);

template <typename T, bool REFCNT>
inline void sub_in(sparsevec_pairvec_t<T,REFCNT> &sv1, const sparsevec_pairvec_t<T,REFCNT> &sv2);

template <typename T, bool REFCNT>
inline void mul_in(sparsevec_pairvec_t<T,REFCNT> &sv, const T &val);

template <typename T, bool REFCNT>
inline void div_in(sparsevec_pairvec_t<T,REFCNT> &sv, const T &val);
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// All sparse vector representation template function
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename A, typename B>
inline void copy_to(A &dst, const B &src);

template <typename T, typename A, typename B>
inline void copy_to(A &dst, const Arageli::vector<T> &src);

template <typename T, typename A>
inline void copy_to(Arageli::vector<T> &dst, const A &src);

template <typename A, typename B>
inline bool equal(const A &sv1, const B &sv2);

template <typename A, typename B>
inline bool non_equal(const A &sv1, const B &sv2);

template <typename A, typename B>
inline void add(A &res, const A &sv1, const B &sv2);

template <typename A, typename B>
inline void sub(A &res, const A &sv1, const B &sv2);

template <typename T, typename A>
inline void mul(A &res, const A &sv, const T &val);

template <typename T, typename A>
inline void div(A &res, const A &sv, const T &val);

template <typename T, typename A, typename B>
inline void dot(T &res, const A &sv1, const B &sv2);

template <typename A, typename B>
inline void add_in(A &sv1, const B &sv2);

template <typename A, typename B>
inline void sub_in(A &sv1, const B &sv2);

template <typename T, typename A>
inline void mul_in(A &sv, const T &val);

template <typename T, typename A>
inline void div_in(A &sv, const T &val);

template <typename A, typename B>
inline size_t struct_intersection(const A &sv1, const B &sv2);

template <typename A, typename B>
inline size_t struct_union(const A &sv1, const B &sv2);

template <typename A, typename B>
inline bool struct_disjoint(const A &sv1, const B &sv2);

}// namespace sprsvc_comp

namespace cnfg
{

/// Default configurator for sparse_matrix class.
/** @param T  type of element of the container;
    @param RepType  type of the representation tag;should be one of
        {mappair_t, vecpair_t, pairvec_t, listpair_t, deqpair_t};
    @param REFCNT  reference counter swither for the entire container.
*/
template<typename T, typename RepType = vecpair_t, bool REFCNT = true>
struct sparse_vector_default
{
};

template <typename T, bool REFCNT>
struct sparse_vector_default< T, vecpair_t, REFCNT >
{
    static const bool refcounting = REFCNT;
    typedef vecpair_t representation_tag;
        typedef typename sprsvc_rep::sparsevec_vecpair_t<T, REFCNT> representation;
};

template <typename T, bool REFCNT>
struct sparse_vector_default< T, pairvec_t, REFCNT >
{
    static const bool refcounting = REFCNT;
    typedef vecpair_t representation_tag;
        typedef typename sprsvc_rep::sparsevec_pairvec_t<T, REFCNT> representation;
};

template <typename T, bool REFCNT>
struct sparse_vector_default< T, mappair_t, REFCNT >
{
    static const bool refcounting = REFCNT;
    typedef vecpair_t representation_tag;
        typedef typename sprsvc_rep::sparsevec_mappair_t<T, REFCNT> representation;
};

}

template <typename T, typename RepType, bool REFCNT>
struct type_traits<cnfg::sparse_vector_default<T, RepType, REFCNT> > :
    public type_traits_default<cnfg::sparse_vector_default<T, RepType, REFCNT> >
{
    static const bool is_specialized = true;
    typedef type_category::configurator category_type;
    static const category_type category_value;
};


namespace spct
{

/// Default spectator for sparse_vector class.
struct sparse_vector_idler
{};

}


template <>
struct type_traits<spct::sparse_vector_idler> :
    public type_traits_default<spct::sparse_vector_idler>
{
    static const bool is_specialized = true;
    typedef type_category::spectator category_type;
    static const category_type category_value;
};


namespace _Internal
{


template <typename T, typename Param, typename Category>
struct sparse_vector_param_extractor_helper_1
{
};


/// Gives Param as configurator and set the default spectator.
template <typename T, typename Param>
struct sparse_vector_param_extractor_helper_1<T, Param, type_category::configurator>
{
    /// Choosen configurator.
    typedef Param configurator;

    /// Default spectator.
    typedef spct::sparse_vector_idler spectator;
};


/// Gives Param as spectator and set the default configurator.
template <typename T, typename Param>
struct sparse_vector_param_extractor_helper_1<T, Param, type_category::spectator>
{
    /// Default configurator.
    typedef cnfg::sparse_vector_default<T> configurator;

    /// Choosen spectator.
    typedef Param spectator;
};


/// Uses Param as a tag to choose the representation type.
template <typename T, typename Param>
struct sparse_vector_param_extractor_helper_1<T, Param, type_category::tag>
{
    /// Default configurator with choosen type of representation.
    typedef cnfg::sparse_vector_default<T, Param> configurator;

    /// Choosen spectator.
    typedef spct::sparse_vector_idler spectator;
};


template <typename T, typename Param>
struct sparse_vector_param_extractor
{
    /// Determined configurator.
    typedef typename sparse_vector_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::configurator
        configurator;

    /// Determined spectator.
    typedef typename sparse_vector_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::spectator
        spectator;
};


template <typename T>
struct sparse_vector_param_extractor<T, default_tag>
{
    /// Default configurator.
    typedef cnfg::sparse_vector_default<T> configurator;

    /// Default spectator.
    typedef spct::sparse_vector_idler spectator;
};

template <typename T>
struct sparse_vector_param_extractor<T, vecpair_t>
{
    /// Default configurator.
    typedef cnfg::sparse_vector_default<T, vecpair_t> configurator;

    /// Default spectator.
    typedef spct::sparse_vector_idler spectator;
};

template <typename T>
struct sparse_vector_param_extractor<T, pairvec_t>
{
    /// Default configurator.
    typedef cnfg::sparse_vector_default<T, pairvec_t> configurator;

    /// Default spectator.
    typedef spct::sparse_vector_idler spectator;
};

template <typename T>
struct sparse_vector_param_extractor<T, mappair_t>
{
    /// Default configurator.
    typedef cnfg::sparse_vector_default<T, mappair_t> configurator;

    /// Default spectator.
    typedef spct::sparse_vector_idler spectator;
};

}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Sparse vector with elements of type T.
/** Way of implementation and other configuration parameters are defined
    by the second template parameter (Param). */
template
<
    typename T,
    typename Param = default_tag
>
class sparse_vector
{
    typedef typename _Internal::sparse_vector_param_extractor<T, Param>::configurator configurator;
    typedef typename configurator::representation Rep;

public:
    typedef typename T value_type;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = configurator::refcounting;

    sparse_vector(size_t _n = 0);

    template<typename T1, typename Param1>
    sparse_vector(const sparse_vector<T1,Param1> &sp);

    sparse_vector(const Arageli::vector<T> &vec);

    ~sparse_vector();

    typedef typename Rep::raw_iterator          raw_iterator;
    typedef typename Rep::const_raw_iterator    const_raw_iterator;
    typedef typename Rep::indraw_iterator       indraw_iterator;
    typedef typename Rep::const_indraw_iterator const_indraw_iterator;
    typedef typename Rep::const_ind_iterator    const_ind_iterator;

    indraw_iterator begin()
    {
        return rep.begin();
    }

    indraw_iterator end()
    {
        return rep.end();
    }

    const_indraw_iterator begin() const
    {
        return rep.begin();
    }

    const_indraw_iterator end() const
    {
        return rep.end();
    }

    indraw_iterator find(size_t i)
    {
        return rep.find(i);
    }

    const_indraw_iterator find(size_t i) const
    {
        return rep.find(i);
    }

    void resize(size_t _n);

    void reserve_mem(size_t elem_count);

    void resize_mem_struct(size_t elem_count);

    void regularize()
    {
        rep.regularize();
    }

    void clear()
    {
        rep.clear();
    }

    const size_t size() const
    {
        return rep.size();
    }

    const size_t real_size() const
    {
        return rep.real_size();
    }

    const size_t capacity() const
    {
        return rep.capacity();
    }

    void insert(size_t i, const T &t)
    {
        rep.insert(i, t);
    }

    void set_elem(size_t real_ind, size_t ind, const T & t)
    {
        rep.set_elem(real_ind, ind, t);
    }

    void push_back(size_t ind, const T & t)
    {
        rep.push_back(ind, t);
    }

    const T & operator[](size_t i) const;

    // WARNING!!! do not change template for second vector
    // Use function copy_to instead
    template<typename T1, typename Param1>
    sparse_vector<T,Param> & operator= (const sparse_vector<T1,Param1> &sv);

    sparse_vector<T,Param> & operator= (const sparse_vector &sv);

    sparse_vector<T,Param> & operator= (const Arageli::vector<T> &vec);

    sparse_vector<T,Param> & operator= (const char *ch);

    template<typename T1, typename Param1>
    sparse_vector<T,Param> & operator+=(const sparse_vector<T1,Param1> &sv);

    template<typename T1, typename Param1>
    sparse_vector<T,Param> & operator-=(const sparse_vector<T1,Param1> &sv);

    sparse_vector<T,Param> & operator*=(const T &val);

    sparse_vector<T,Param> & operator/=(const T &val);

    operator vector<T> ();

    bool dynamic() const
    {
        return rep.dynamic();
    }

    void swap(sparse_vector<T,Param> &x)
    {
        rep.swap(x.rep);
    }

protected:

    void raw_set_elem(size_t real_ind, size_t ind, const T & t)
    {
        rep.raw_set_elem(real_ind, ind, t);
    }

#ifdef ARAGELI_DEBUG_LEVEL_3
    bool is_ordered() const
    {
        return rep.is_ordered();
    }
#endif

protected:
    Rep rep;
private:
    //friend section
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    template<typename T1, typename Param1>
    friend class sparse_vector;
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool operator==(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool operator!=(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_vector<T1,Param1> operator+ (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_vector<T1,Param1> operator- (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const big_int &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const big_int &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const float &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const float &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const double &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const double &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const rational<big_int> &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const rational<big_int> &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator* (const sparse_vector<T,Param> &sv, const big_float &val);

    template<typename T, typename Param>
    friend sparse_vector<T,Param> operator/ (const sparse_vector<T,Param> &sv, const big_float &val);
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend T1 operator* (const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T, typename T1, typename Param1>
    friend T1 operator* (const sparse_vector<T1,Param1> &sv1, const Arageli::vector<T> &vec);

    template<typename T,typename Param>
    friend std::ostream & operator <<
    (
        std::ostream & s,
        const sparse_vector<T,Param> & x
    );

    template<typename T,typename Param>
    friend std::istream & operator >>
    (
        std::istream & s,
        sparse_vector<T,Param> & x
    );

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    template<typename T1, typename Param1, typename T2, typename Param2>
    friend void copy_to(sparse_vector<T1,Param1> &dst, const sparse_vector<T2,Param2> &src);

    template <typename T,typename Param1, typename T1>
    friend void copy_to(sparse_vector<T1,Param1> &dst, const Arageli::vector<T> &src);

    template <typename T,typename Param1, typename T1>
    friend void copy_to(Arageli::vector<T> &dst, const sparse_vector<T1,Param1> &src);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool equal(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool non_equal(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend void add(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend void sub(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template <typename T, typename T1,typename Param1>
    friend void mul(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv, const T &val);

    template <typename T, typename T1,typename Param1>
    friend void div(sparse_vector<T1,Param1> &res, const sparse_vector<T1,Param1> &sv, const T &val);

    template <typename T, typename T1, typename Param1, typename T2, typename Param2>
    friend void dot(T &res, const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T, typename Param, typename T1,  typename T2>
    friend void dot(T1 &res, const sparse_vector<T,Param> &sv, const Arageli::vector<T2> &vec);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend void add_in(sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend void sub_in(sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template <typename T, typename T1, typename Param1>
    friend void mul_in(sparse_vector<T1,Param1> &sv, const T &val);

    template <typename T, typename T1, typename Param1>
    friend void div_in(sparse_vector<T1,Param1> &sv, const T &val);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend size_t struct_intersection(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool struct_disjoint(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend size_t struct_union(const sparse_vector<T1,Param1> &sv1, const sparse_vector<T2,Param2> &sv2);
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
};

/// IO  sparse vector operations
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// IO sparse vector params
extern const char* sparse_vector_output_list_first_bracket_default;
extern const char* sparse_vector_output_list_second_bracket_default;
extern const char* sparse_vector_output_list_separator_default;
extern const char* sparse_vector_output_list_null_symbol_default;
extern const char* sparse_vector_output_list_first_pair_bracket_default;
extern const char* sparse_vector_output_list_second_pair_bracket_default;
extern const char* sparse_vector_output_list_pair_separator_default;
extern const char* sparse_vector_output_standart_order_error_default;
extern const char* sparse_vector_input_list_first_bracket_default;
extern const char* sparse_vector_input_list_second_bracket_default;
extern const char* sparse_vector_input_list_separator_default;
extern const char* sparse_vector_input_list_range_default;


/// Simple output of the sparse vector
template <typename Out, typename T, typename Param>
Out& output_list1
(
    Out& out,
    const sparse_vector<T, Param>& x,
    const char* first_bracket = sparse_vector_output_list_first_bracket_default,
    const char* second_bracket = sparse_vector_output_list_second_bracket_default,
    const char* separator = sparse_vector_output_list_separator_default,
    const char* first_pair_bracket = sparse_vector_output_list_first_pair_bracket_default,
    const char* second_pair_bracket = sparse_vector_output_list_second_pair_bracket_default,
    const char* second_pair_separator = sparse_vector_output_list_pair_separator_default
);

/// Simple output of the sparse vector as dense vector
template <typename Out, typename T, typename Param>
Out& output_list
(
    Out& out,
    const sparse_vector<T, Param>& x,
    const char* first_bracket = sparse_vector_output_list_first_bracket_default,
    const char* second_bracket = sparse_vector_output_list_second_bracket_default,
    const char* separator = sparse_vector_output_list_separator_default,
    const char* null_symbol = sparse_vector_output_list_null_symbol_default,
    const char* index_order_error = sparse_vector_output_standart_order_error_default
);

/// Simple input of the sparse vector
template <typename In, typename T, typename Param>
In& input_list
(
    In& in,
    sparse_vector<T, Param>& x,
    const char* first_bracket = sparse_vector_input_list_first_bracket_default,
    const char* second_bracket = sparse_vector_input_list_second_bracket_default,
    const char* separator = sparse_vector_input_list_separator_default,
    const char* range = sparse_vector_input_list_range_default
);


}// namespace Arageli

#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SPARSE_VECTOR
    #include "sparse_vector.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SPARSE_VECTOR
#endif

#endif    // #ifndef _ARAGELI_sparse_vector_hpp_
