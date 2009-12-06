/*****************************************************************************

    sparse_matrix.hpp

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
    \file sparse_matrix.hpp
    \brief Sparse matrix and supporting tools declaration.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#ifndef _ARAGELI_sparse_matrix_hpp_
#define _ARAGELI_sparse_matrix_hpp_

#include <cstddef>
#include <utility>

#include "config.hpp"

#include "type_traits.hpp"

#include "misc.hpp"
#include "sparse_vector.hpp"

// REFERENCE ADDITIONAL HEADERS HERE

namespace Arageli
{

namespace iterator_base
{

struct real
{};
struct idler
{};

template<typename Arg>
class Indexer
{
};

template<>
class Indexer<real>
{
    size_t i;
public:
    Indexer()
    {}
    Indexer(size_t index)
        : i(index)
    {}
    void init(size_t index)
    {
        i = index;
    }
    void operator++()
    {
        ++i;
    }
    void operator--()
    {
        --i;
    }
    void get_index(size_t &index) const
    {
        index = i;
    }
};

template<>
class Indexer<idler>
{
public:
    Indexer()
    {}
    Indexer(size_t index)
    {}
    void init(size_t index)
    {}
    void operator++()
    {}
    void operator--()
    {}
    void get_index(size_t &index) const
    {}
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
class accross_iterator_base
{
protected:
    ElemIterator it;
    VecIterator vec;
#ifdef ARAGELI_DEBUG_LEVEL_3
    VecIterator vec_beg;
#endif
    VecIterator vec_end;
    size_t vec_ind;
    Indexer ix;
public:
    accross_iterator_base()
    {
    }

    accross_iterator_base(ElemIterator &it, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t vec_ind) :
        it(it),
        vec(vec),
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg(vec_beg),
#endif
        vec_end(vec_end),
        vec_ind(vec_ind)
    {
    }

    template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
    accross_iterator_base(const accross_iterator_base <ElemIterator, VecIterator, Indexer, ValType> &i)
    {
        (*this)= i;
    }

    template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
    accross_iterator_base & operator= (const accross_iterator_base <ElemIterator, VecIterator, Indexer, ValType> &i)
    {
        it = i.it;
        vec = i.vec;
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg = i.vec_beg;
#endif
        vec_end = i.vec_end;
        vec_ind = i.vec_ind;
        ix = i.ix;
        return *this;
    }

    bool operator==(const accross_iterator_base &i)
    {
        return vec==i.vec && (vec==vec_end || it==i.it);
    }

    bool operator!=(const accross_iterator_base &i)
    {
        return vec!=i.vec || (vec!=vec_end && it!=i.it);
    }

    void operator++()
    {
        ARAGELI_ASSERT_0(vec!=vec_end);

        do
        {
            ++ix;
            ++vec;
            if(vec!=vec_end)
            {
                it=vec->find(vec_ind);
            }
            else
            {
                break;
            }
        }
        while(it==vec->end());
    }

    void operator--()
    {
        ARAGELI_ASSERT_3(vec!=vec_beg);

        do
        {
            --ix;
            --vec;
            it=vec->find(vec_ind);
            ARAGELI_ASSERT_3(vec!=vec_beg);
        }
        while(it==vec->end());
    }

    accross_iterator_base   operator++(int)
    {
        accross_iterator_base tmp(it);
        ++(*this);
        return tmp;
    }

    accross_iterator_base & operator--(int)
    {
        accross_iterator_base tmp(it);
        --(*this);
        return tmp;
    }

    template<typename A>
    void advance(int k)
    {
        int i=0;
        if(k<0)
        {
            while(i!=k)
            {
                --(*this);
                --ix;
            }
        }
        else
        {
            while(i!=k)
            {
                ++(*this);
                ++ix;
            }
        }
    }
protected:
    size_t ind() const
    {
        size_t i;
        ix.get_index(i);
        return i;
    }

    ValType & el() const
    {
        return it.el();
    }
};

template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
class accross_indraw_iterator_base :
    public accross_iterator_base<ElemIterator, VecIterator, Indexer, ValType>
{
public:
    accross_indraw_iterator_base()
    {
    }

    accross_indraw_iterator_base(ElemIterator &it, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_ind, size_t ind) :
        accross_iterator_base(it, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
            vec_beg,
#endif
            vec_end, line_ind)
    {
        ix = ind;
    }

    size_t ind() const
    {
        return accross_iterator_base::ind();
    }

    ValType & el() const
    {
        return accross_iterator_base::el();
    }
};

template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
class accross_raw_iterator_base
    : public accross_iterator_base<ElemIterator, VecIterator, Indexer, ValType>
{
public:
    accross_raw_iterator_base()
    {}
    accross_raw_iterator_base(ElemIterator &it, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_ind, size_t ind) :
        accross_iterator_base(it, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
            vec_beg,
#endif
            vec_end, line_ind)
    {
    }

    ValType & el() const
    {
        return accross_iterator_base::el();
    }
};

template <typename ElemIterator, typename VecIterator, typename Indexer, typename ValType>
class accross_ind_iterator_base :
    public accross_iterator_base<ElemIterator, VecIterator, Indexer, ValType>
{
public:
    accross_ind_iterator_base()
    {}
    accross_ind_iterator_base(ElemIterator &it, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_ind, size_t ind) :
        accross_iterator_base(it, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
            vec_beg,
#endif
            vec_end, line_ind),
        ix(ind)
    {
    }

    size_t ind() const
    {
        return accross_iterator_base::ind();
    }
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class spmat_iterator_base
{
protected:
    LineIterator it;
    VecIterator vec;
#ifdef ARAGELI_DEBUG_LEVEL_3
    VecIterator vec_beg;
#endif
    VecIterator vec_end;
    Indexer ix;
public:
    spmat_iterator_base()
    {
    }

    spmat_iterator_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end) :
        it(i),
        vec(vec),
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg(vec_beg),
#endif
        vec_end(vec_end)
    {
    }

    template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
    spmat_iterator_base & operator= (const spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType> &i)
    {
        it = i.it;
        vec = i.vec;
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg = i.vec_beg;
#endif
        vec_end = i.vec_end;
        ix = i.ix;
        return *this;
    }

    bool operator==(const spmat_iterator_base &i)
    {
        return vec==i.vec && (vec==vec_end || it==i.it);
    }

    bool operator!=(const spmat_iterator_base &i)
    {
        return vec!=i.vec || (vec!=vec_end && it!=i.it);
    }

    void operator++()
    {
        ARAGELI_ASSERT_0(vec!=vec_end);

        ++it;

        while(it==vec->end())
        {
            ++ix;
            ++vec;
            if(vec!=vec_end)
            {
                it=vec->begin();
            }
            else
            {
                break;
            }
        }
    }

    void operator--()
    {
        if(vec!=vec_end && it!=vec->begin())
        {
            --it;
        }
        else
        {
            do
            {
                ARAGELI_ASSERT_3(vec!=vec_beg)
                --ix;
                --vec;
            }
            while(vec->begin()==vec->end());
            it=vec->end();
            --it;
        }
    }

    spmat_iterator_base   operator++(int)
    {
        spmat_iterator_base tmp(it);
        ++(*this);
        return tmp;
    }

    spmat_iterator_base & operator--(int)
    {
        spmat_iterator_base tmp(it);
        --(*this);
        return tmp;
    }

    template<typename A>
    void advance(int k, A &bx)
    {
        int i=0;
        if(k<0)
        {
            while(i!=k)
            {
                --(*this);
                --ix;
            }
        }
        else
        {
            while(i!=k)
            {
                ++(*this);
                ++ix;
            }
        }
    }
protected:
    ValType & el() const
    {
        return it.el();
    }

    size_t ind_line() const
    {
        size_t i;
        ix.get_index(i);
        return i;
    }

    size_t ind_elem() const
    {
        return it.ind();
    }
};

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class raw_iterator_base
    : public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    raw_iterator_base()
    {}
    raw_iterator_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end)
    {
    }

    ValType & el() const
    {
        return spmat_iterator_base::el();
    }
};

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class indraw_iterator_base :
    public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    indraw_iterator_base()
    {}
    indraw_iterator_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end),
        ix(line_index)
    {
    }

    ValType & el() const
    {
        return spmat_iterator_base::el();
    }

    size_t ind_line() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_elem() const
    {
        return spmat_iterator_base::ind_elem();
    }
};

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class ind_iterator_base :
    public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    ind_iterator_base()
    {}
    ind_iterator_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end),
        ix(line_index)
    {
    }

    size_t ind_line() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_elem() const
    {
        return spmat_iterator_base::ind_elem();
    }
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class indraw_iterator_for_row_base
    : public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    indraw_iterator_for_row_base()
    {}
    indraw_iterator_for_row_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end)
    {
        ix = line_index;
    }

    ValType & el() const
    {
        return spmat_iterator_base::el();
    }

    size_t ind_row() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_col() const
    {
        return spmat_iterator_base::ind_elem();
    }
};

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class ind_iterator_for_row_base
    : public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    ind_iterator_for_row_base()
    {}
    ind_iterator_for_row_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end),
        ix(line_index)
    {
    }

    size_t ind_row() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_col() const
    {
        return spmat_iterator_base::ind_elem();
    }
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class indraw_iterator_for_col_base :
    public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    indraw_iterator_for_col_base()
    {}
    indraw_iterator_for_col_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end)
    {
        ix = line_index;
    }

    ValType & el() const
    {
        return spmat_iterator_base::el();
    }

    size_t ind_col() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_row() const
    {
        return spmat_iterator_base::ind_elem();
    }
};

template <typename LineIterator, typename VecIterator, typename Indexer, typename ValType>
class ind_iterator_for_col_base
    : public spmat_iterator_base<LineIterator, VecIterator, Indexer, ValType>
{
public:
    ind_iterator_for_col_base()
    {}
    ind_iterator_for_col_base(const LineIterator &i, VecIterator vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        VecIterator vec_beg,
#endif
        VecIterator vec_end, size_t line_index) :
        spmat_iterator_base(i, vec,
#ifdef ARAGELI_DEBUG_LEVEL_3
        vec_beg,
#endif
        vec_end),
        ix(line_index)
    {
    }

    size_t ind_col() const
    {
        return spmat_iterator_base::ind_line();
    }

    size_t ind_row() const
    {
        return spmat_iterator_base::ind_elem();
    }
};



}// iterator_base


struct vecpair_col_t {};

ARAGELI_REGISTER_TAG(vecpair_col_t)

struct pairvec_col_t {};

ARAGELI_REGISTER_TAG(pairvec_col_t)

struct mappair_col_t {};

ARAGELI_REGISTER_TAG(mappair_col_t)

struct vecpair_row_t {};

ARAGELI_REGISTER_TAG(vecpair_row_t)

struct pairvec_row_t {};

ARAGELI_REGISTER_TAG(pairvec_row_t)

struct mappair_row_t {};

ARAGELI_REGISTER_TAG(mappair_row_t)

namespace spmt_rep
{

enum spmt_form_t { colform_t, rowform_t };

using namespace iterator_base;

template <typename T, typename ColRep, bool REFCNT = true>
class sparsemat_colform_t
{
public:
    typedef typename T value_type;
    typedef typename sparse_vector<T,ColRep> RepLine;
    typedef typename std::vector<RepLine> Rep;
    typedef typename spmt_form_t Form;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = REFCNT;
public:

    sparsemat_colform_t(size_t nrows = 0, size_t ncols = 0);

    size_t nrows() const
    {
        return n;
    }

    size_t ncols() const
    {
        return m;
    }

    size_t real_size() const;

    size_t capacity() const;

    void resize(size_t nrows, size_t ncols);

    void reserve_line_mem(size_t line_ind, size_t line_size);

    void reserve_mem(size_t &lines_size);

    void reserve_mem(Arageli::vector<size_t> &lines_size);

    void resize_mem_struct_line(size_t line_ind, size_t line_size);

    void resize_mem_struct(size_t &lines_size);

    void resize_mem_struct(Arageli::vector<size_t> &lines_size);

    void resize_struct(size_t m, size_t n);

    void regularize();

    void clear();

    void insert(size_t i, size_t j, const T &t)
    {
        ARAGELI_ASSERT_0(i<nrows()&& j<ncols());

        if(!null(t))
        {
            ins(i,j,t);
        }
        else
        {
            del(i,j);
        }
    }

    template<typename T1, typename Param1>
    void set_row(size_t row_in, const sparse_vector<T1,Param1> &row);

    template<typename T1>
    void set_row(size_t row_in, const Arageli::vector<T1> &row);

    template<typename T1, typename Param1>
    void set_col(size_t col_in, const sparse_vector<T1,Param1> &col);

    template<typename T1>
    void set_col(size_t col_in, const Arageli::vector<T1> &col);

    template<typename T1, typename Param1>
    void set_line(size_t line_in, const sparse_vector<T1,Param1> &line)
    {
        set_col(line_in, line);
    }

    template<typename T1>
    void set_line(size_t line_in, const Arageli::vector<T1> &line)
    {
        set_col(line_in, line);
    }

    void set_elem(size_t col, size_t real_row, size_t row, const T & t)
    {
        rep[col].set_elem(real_row, row, t);
    }

    void raw_set_elem(size_t col, size_t real_row, size_t row, const T & t)
    {
        rep[col].raw_set_elem(real_row, row, t);
    }

    void line_push_back(size_t col, size_t row, const T & t)
    {
        rep[col].push_back(row, t);
    }

    const sparse_vector<T,ColRep> & get_line(size_t i) const
    {
        return rep[i];
    }

    bool is_ordered() const;

    bool rowform() const
    {
        return false;
    }

    bool colform() const
    {
        return true;
    }

    bool dynamic() const
    {
        RepLine x;
        return x.dynamic();
    }

    spmt_form_t rep_form() const
    {
        return colform_t;
    }

    bool is_empty() const
    {
        return m==0 || n==0;
    }

    void swap(sparsemat_colform_t<T,ColRep,REFCNT> &x)
    {
        rep.swap(x.rep);
        std::swap(m,x.m);
        std::swap(n,x.n);
    }

    const T & get(size_t i, size_t j) const
    {
        return rep[j][i];
    }

protected:

    typedef typename Rep::iterator iterator;

    typedef typename Rep::const_iterator const_iterator;

public:

    typedef typename RepLine::indraw_iterator col_indraw_iterator;

    typedef typename RepLine::raw_iterator col_raw_iterator;

    typedef typename accross_indraw_iterator_base<col_indraw_iterator, iterator, Indexer<real>, T> row_indraw_iterator;

    typedef typename accross_raw_iterator_base<col_raw_iterator, iterator, Indexer<idler>, T> row_raw_iterator;

    typedef typename indraw_iterator_for_col_base<col_indraw_iterator,iterator,Indexer<real>, T> indraw_iterator;

    typedef typename raw_iterator_base<col_raw_iterator, iterator, Indexer<idler>, T> raw_iterator;

    typedef typename RepLine::const_indraw_iterator const_col_indraw_iterator;

    typedef typename RepLine::const_raw_iterator const_col_raw_iterator;

    typedef typename RepLine::const_ind_iterator const_col_ind_iterator;

    typedef typename accross_indraw_iterator_base<const_col_indraw_iterator, const_iterator, Indexer<real>, const T> const_row_indraw_iterator;

    typedef typename accross_raw_iterator_base<const_col_raw_iterator, const_iterator,Indexer<idler>, const T> const_row_raw_iterator;

    typedef typename accross_ind_iterator_base<const_col_ind_iterator, const_iterator, Indexer<real>, const T> const_row_ind_iterator;

    typedef typename indraw_iterator_for_col_base<const_col_indraw_iterator, const_iterator, Indexer<real>, const T> const_indraw_iterator;

    typedef typename raw_iterator_base<const_col_raw_iterator, const_iterator, Indexer<idler>, const T> const_raw_iterator;

    typedef typename ind_iterator_for_col_base<const_col_indraw_iterator, const_iterator, Indexer<idler>, const T> const_ind_iterator;

public:
    row_indraw_iterator begin_row(size_t row)
    {
        RepLine::indraw_iterator col_it;
        Rep::iterator rep_it;
        size_t col_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                col_it=rep_it->find(row);
                if(col_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++col_ind;
            }
            while(rep_it!=rep.end());
        }

        return row_indraw_iterator(col_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_it,
#endif
                                    rep.end(), row, col_ind);
    }

    const_row_indraw_iterator begin_row(size_t row) const
    {
        RepLine::const_indraw_iterator col_it;
        Rep::const_iterator rep_it;
        size_t col_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                col_it=rep_it->find(row);
                if(col_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++col_ind;
            }
            while(rep_it!=rep.end());
        }

        return const_row_indraw_iterator(col_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                        rep_it,
#endif
                                        rep.end(), row, col_ind);
    }

    row_indraw_iterator end_row(size_t row)
    {
        RepLine::indraw_iterator col_it;
        Rep::iterator rep_begin;
        Rep::iterator rep_it;

        rep_begin=rep.begin();
        if(rep_begin!=rep.end())
        {
            do
            {
                col_it=rep_begin->find(row);
                if(col_it!=rep_begin->end())
                {
                    break;
                }
                ++rep_begin;
            }
            while(rep_begin!=rep.end());
        }

        return row_indraw_iterator(col_it, rep.end(),
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_begin,
#endif
                                    rep.end(), row, rep.size());
    }

    const_row_indraw_iterator end_row(size_t row) const
    {
        RepLine::const_indraw_iterator col_it;
        Rep::const_iterator rep_begin;
        Rep::const_iterator rep_it;

        rep_begin=rep.begin();
        if(rep_begin!=rep.end())
        {
            do
            {
                col_it=rep_begin->find(row);
                if(col_it!=rep_begin->end())
                {
                    break;
                }
                ++rep_begin;
            }
            while(rep_begin!=rep.end());
        }

        return const_row_indraw_iterator(col_it, rep.end(),
#ifdef ARAGELI_DEBUG_LEVEL_3
                                        rep_begin,
#endif
                                        rep.end(), row, rep.size());
    }

    col_indraw_iterator begin_col(size_t col)
    {
        return rep[col].begin();
    }

    const_col_indraw_iterator begin_col(size_t col) const
    {
        return rep[col].begin();
    }

    col_indraw_iterator end_col(size_t col)
    {
        return rep[col].end();
    }

    const_col_indraw_iterator end_col(size_t col) const
    {
        return rep[col].end();
    }

    indraw_iterator begin()
    {
        RepLine::indraw_iterator col_it;
        Rep::iterator rep_it;
        size_t col_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                col_it=rep_it->begin();
                if(col_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++col_ind;
            }
            while(rep_it!=rep.end());
        }
        return indraw_iterator(col_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
            rep_it,
#endif
            rep.end(), col_ind);
    }

    const_indraw_iterator begin() const
    {
        RepLine::const_indraw_iterator col_it;
        Rep::const_iterator rep_it;
        size_t col_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                col_it=rep_it->begin();
                if(col_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++col_ind;
            }
            while(rep_it!=rep.end());
        }
        return const_indraw_iterator(col_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_it,
#endif
                                    rep.end(), col_ind);
    }

    indraw_iterator end()
    {
        RepLine::indraw_iterator col_it;
        Rep::iterator rep_it;
        Rep::iterator rep_begin=rep.begin();
        rep_it=rep.end();
        return indraw_iterator(col_it, rep_it, rep_begin, rep_it, ncols());
    }

    const_indraw_iterator end() const
    {
        RepLine::const_indraw_iterator col_it;
        Rep::const_iterator rep_it;
        Rep::const_iterator rep_begin=rep.begin();
        rep_it=rep.end();
        return const_indraw_iterator(col_it, rep_it, rep_begin, rep_it, ncols());
    }

public:

    Rep rep;

    size_t m;

    size_t n;

private:

    void ins(size_t i, size_t j, const T &t)
    {
        rep[j].insert(i,t);
    }

    void del(size_t i, size_t j)
    {
        rep[j].insert(i,0);
    }
};


template <typename T, typename RowRep, bool REFCNT = true>
class sparsemat_rowform_t
{
public:
    typedef typename T value_type;
    typedef typename sparse_vector<T,RowRep> RepLine;
    typedef typename std::vector<RepLine> Rep;
    typedef typename spmt_form_t Form;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = REFCNT;
public:

    sparsemat_rowform_t(size_t ncols = 0, size_t nrows = 0);

    size_t ncols() const
    {
        return n;
    }

    size_t nrows() const
    {
        return m;
    }

    size_t real_size() const;

    size_t capacity() const;

    void resize(size_t ncols, size_t nrows);

    void reserve_line_mem(size_t line_ind, size_t line_size);

    void reserve_mem(size_t &lines_size);

    void reserve_mem(Arageli::vector<size_t> &lines_size);

    void resize_mem_struct_line(size_t line_ind, size_t line_size);

    void resize_mem_struct(size_t &lines_size);

    void resize_mem_struct(Arageli::vector<size_t> &lines_size);

    void resize_struct(size_t m, size_t n);

    void regularize();

    void clear();

    void insert(size_t i, size_t j, const T &t)
    {
        ARAGELI_ASSERT_0(i<nrows()&& j<ncols());

        if(!null(t))
        {
            ins(i,j,t);
        }
        else
        {
            del(i,j);
        }
    }

    template<typename T1, typename Param1>
    void set_row(size_t row_in, const sparse_vector<T1,Param1> &row);

    template<typename T1>
    void set_row(size_t row_in, const Arageli::vector<T1> &row);

    template<typename T1, typename Param1>
    void set_col(size_t col_in, const sparse_vector<T1,Param1> &col);

    template<typename T1>
    void set_col(size_t col_in, const Arageli::vector<T1> &col);

    template<typename T1, typename Param1>
    void set_line(size_t line_in, const sparse_vector<T1,Param1> &line)
    {
        set_row(line_in, line);
    }

    template<typename T1>
    void set_line(size_t line_in, const Arageli::vector<T1> &line)
    {
        set_row(line_in, line);
    }

    void set_elem(size_t row, size_t real_col, size_t col, const T & t)
    {
        rep[row].set_elem(real_col, col, t);
    }

    void raw_set_elem(size_t row, size_t real_col, size_t col, const T & t)
    {
        rep[row].raw_set_elem(real_col, col, t);
    }

    void line_push_back(size_t row, size_t col, const T & t)
    {
        rep[row].push_back(col, t);
    }

    const sparse_vector<T,RowRep> & get_line(size_t i) const
    {
        return rep[i];
    }

    bool is_ordered() const;

    bool rowform() const
    {
        return true;
    }

    bool colform() const
    {
        return false;
    }

    bool dynamic() const
    {
        RepLine x;
        return x.dynamic();
    }

    spmt_form_t rep_form() const
    {
        return rowform_t;
    }

    bool is_empty() const
    {
        return m==0 || n==0;
    }

    void swap(sparsemat_rowform_t<T,RowRep,REFCNT> &x)
    {
        rep.swap(x.rep);
        std::swap(m,x.m);
        std::swap(n,x.n);
    }

    const T & get(size_t i, size_t j) const
    {
        return rep[i][j];
    }

protected:

    typedef typename Rep::iterator iterator;
    typedef typename Rep::const_iterator const_iterator;

public:
    typedef typename RepLine::indraw_iterator row_indraw_iterator;

    typedef typename RepLine::raw_iterator row_raw_iterator;

    typedef typename accross_indraw_iterator_base<row_indraw_iterator, iterator, Indexer<real>, T> col_indraw_iterator;

    typedef typename accross_raw_iterator_base<row_raw_iterator, iterator, Indexer<idler>, T> col_raw_iterator;

    typedef typename indraw_iterator_for_row_base<row_indraw_iterator,iterator,Indexer<real>, T> indraw_iterator;

    typedef typename raw_iterator_base<row_raw_iterator, iterator, Indexer<idler>, T> raw_iterator;

    typedef typename RepLine::const_indraw_iterator const_row_indraw_iterator;

    typedef typename RepLine::const_raw_iterator const_row_raw_iterator;

    typedef typename RepLine::const_ind_iterator const_row_ind_iterator;

    typedef typename accross_indraw_iterator_base<const_row_indraw_iterator, const_iterator,Indexer<real>, const T> const_col_indraw_iterator;

    typedef typename accross_raw_iterator_base<const_row_raw_iterator, const_iterator,Indexer<idler>, const T> const_col_raw_iterator;

    typedef typename accross_ind_iterator_base<const_row_ind_iterator, const_iterator, Indexer<real>, const T> const_col_ind_iterator;

    typedef typename indraw_iterator_for_row_base<const_row_indraw_iterator,const_iterator,Indexer<real>,const T> const_indraw_iterator;

    typedef typename raw_iterator_base<const_row_raw_iterator,const_iterator,Indexer<idler>,const T> const_raw_iterator;

    typedef typename ind_iterator_for_row_base<const_row_indraw_iterator,const_iterator, Indexer<idler>, const T> const_ind_iterator;

public:
    col_indraw_iterator begin_col(size_t col)
    {
        RepLine::indraw_iterator row_it;
        Rep::iterator rep_it;
        size_t row_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                row_it=rep_it->find(col);
                if(row_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++row_ind;
            }
            while(rep_it!=rep.end());
        }

        return col_indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_it,
#endif
                                    rep.end(), col, row_ind);
    }

    const_col_indraw_iterator begin_col(size_t col) const
    {
        RepLine::const_indraw_iterator row_it;
        Rep::const_iterator rep_it;
        size_t row_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                row_it=rep_it->find(col);
                if(row_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++row_ind;
            }
            while(rep_it!=rep.end());
        }

        return const_col_indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                        rep_it,
#endif
                                        rep.end(), col, row_ind);
    }

    col_indraw_iterator end_col(size_t col)
    {
        RepLine::indraw_iterator row_it;
        Rep::iterator rep_begin;
        Rep::iterator rep_it;

        rep_begin=rep.begin();
        if(rep_begin!=rep.end())
        {
            do
            {
                row_it=rep_begin->find(col);
                if(row_it!=rep_begin->end())
                {
                    break;
                }
                ++rep_begin;
            }
            while(rep_begin!=rep.end());
        }

        return col_indraw_iterator(row_it, rep.end(),
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_begin,
#endif
                                    rep.end(), col, rep.size());
    }

    const_col_indraw_iterator end_col(size_t col) const
    {
        RepLine::const_indraw_iterator row_it;
        Rep::const_iterator rep_begin;
        Rep::const_iterator rep_it;

        rep_begin=rep.begin();
        if(rep_begin!=rep.end())
        {
            do
            {
                row_it=rep_begin->find(col);
                if(row_it!=rep_begin->end())
                {
                    break;
                }
                ++rep_begin;
            }
            while(rep_begin!=rep.end());
        }

        return const_col_indraw_iterator(row_it, rep.end(),
#ifdef ARAGELI_DEBUG_LEVEL_3
                                        rep_begin,
#endif
                                        rep.end(), col, rep.size());
    }

    row_indraw_iterator begin_row(size_t row)
    {
        return rep[row].begin();
    }

    const_row_indraw_iterator begin_row(size_t row) const
    {
        return rep[row].begin();
    }

    row_indraw_iterator end_row(size_t row)
    {
        return rep[row].end();
    }

    const_row_indraw_iterator end_row(size_t row) const
    {
        return rep[row].end();
    }

    indraw_iterator begin()
    {
        RepLine::indraw_iterator row_it;
        Rep::iterator rep_it;
        size_t row_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                row_it=rep_it->begin();
                if(row_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++row_ind;
            }
            while(rep_it!=rep.end());
        }
        return indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                rep_it,
#endif
                                rep.end(), row_ind);
    }

    const_indraw_iterator begin() const
    {
        RepLine::const_indraw_iterator row_it;
        Rep::const_iterator rep_it;
        size_t row_ind=0;

        rep_it=rep.begin();

        if(rep_it!=rep.end())
        {
            do
            {
                row_it=rep_it->begin();
                if(row_it!=rep_it->end())
                {
                    break;
                }
                ++rep_it;
                ++row_ind;
            }
            while(rep_it!=rep.end());
        }
        return const_indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_it,
#endif
                                    rep.end(), row_ind);
    }

    indraw_iterator end()
    {
        RepLine::indraw_iterator row_it;
        Rep::iterator rep_it;
        Rep::iterator rep_begin=rep.begin();
        rep_it=rep.end();
        return indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                rep_begin,
#endif
                                rep_it, nrows());
    }

    const_indraw_iterator end() const
    {
        RepLine::const_indraw_iterator row_it;
        Rep::const_iterator rep_it;
        Rep::const_iterator rep_begin=rep.begin();
        rep_it=rep.end();
        return const_indraw_iterator(row_it, rep_it,
#ifdef ARAGELI_DEBUG_LEVEL_3
                                    rep_begin,
#endif
                                    rep_it, nrows());
    }

public:

    Rep rep;

    size_t m;

    size_t n;

private:

    void ins(size_t i, size_t j, const T &t)
    {
        rep[i].insert(j,t);
    }

    void del(size_t i, size_t j)
    {
        rep[i].insert(j,0);
    }
};

}

/** Namespace of sparse matrix representation computation function
    and sparse matrix representation operation with memory function */
namespace spmt_comp
{
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// All sparse matrix representation template function
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename A, typename B>
void copy_to(A &dst, const B &src);

template <typename T,typename A>
void copy_to(A &dst, const Arageli::matrix<T> &src);

template <typename T, typename A>
void copy_to(matrix<T> &dst, const A &src);

template <typename A,typename B>
void extract_accross_line(A &dst, const B &src);

template <typename A, typename B, typename SV>
void copy_lines(A &dst, const B &src, size_t line_num);

template <typename A, typename B, typename SV>
void extract_accross_lines(A &dst, const B &src, SV ind);

template <typename A, typename B>
bool equal(const A &sm1, const B &sm2);

template <typename A, typename B>
bool non_equal(const A &sm1, const B &sm2);

template <typename A, typename B>
void turn_over_struct(A &res, const B &sm);

template <typename A, typename B>
void transpose(A &res, const B &sm);

template <typename A, typename B>
void add(A &res, const A &sm1, const B &sm2);

template <typename A, typename B>
void sub(A &res, const A &sm1, const B &sm2);

template <typename T, typename A>
void mul(A &res, const A &sm, const T &val);

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_struct(sparse_vector<T1,Param1> &res, const A &sm, const sparse_vector<T2,Param2> &sv);

template <typename A, typename T1, typename T2>
void mul_struct(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec);

template <typename A, typename T1, typename T2>
void mul_struct_accross(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec);

template <typename A, typename B>
void mul_matrices(A &res, const A &sm1, const B &sm2);

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_right(sparse_vector<T1,Param1> &res, const A &sm, const sparse_vector<T2,Param2> &sv);

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_left(sparse_vector<T1,Param1> &res, const sparse_vector<T2,Param2> &sv, const A &sm);

template <typename A, typename T1, typename T2>
void mul_right(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec);

template <typename A, typename T1, typename T2>
void mul_left(Arageli::vector<T1> &res, const Arageli::vector<T2> &vec, const A &sm);

template <typename T, typename A>
void div(A &res, const A &sm, const T &val);

template <typename A, typename B>
void add_in(A &sm1, const B &sm2);

template <typename A, typename B>
void sub_in(A &sm1, const B &sm2);

template <typename T, typename A>
void mul_in(A &sm, const T &val);

template <typename T, typename A>
void div_in(A &sm, const T &val);

template <typename T, typename A>
void mul_line(A &sm, size_t line, const T &val);

template <typename T, typename A>
void mul_accross_line(A &sm, size_t line, const T &val);
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}

namespace cnfg
{

template<typename T, typename RepType = vecpair_col_t, bool REFCNT = true>
struct sparse_matrix_default
{
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,vecpair_col_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef vecpair_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_colform_t<T,vecpair_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, vecpair_t, false> line_configurator;
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,pairvec_col_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef pairvec_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_colform_t<T,pairvec_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, pairvec_t, false> line_configurator;
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,mappair_col_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef mappair_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_colform_t<T,mappair_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, mappair_t, false> line_configurator;
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,vecpair_row_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef vecpair_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_rowform_t<T,vecpair_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, vecpair_t, false> line_configurator;
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,pairvec_row_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef pairvec_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_rowform_t<T,pairvec_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, pairvec_t, false> line_configurator;
};

template<typename T, bool REFCNT>
struct sparse_matrix_default<T,mappair_row_t,REFCNT>
{
    static const bool refcounting = REFCNT;
    typedef mappair_t line_representation_tag;
    typedef typename spmt_rep::sparsemat_rowform_t<T,mappair_t,REFCNT> representation;

    /// Configurator for one row or column of the matrix.
    /** Here we require sparse_matrix to use sparse_vector as a row (column)
        storage. */
    typedef typename cnfg::sparse_vector_default<T, mappair_t, false> line_configurator;
};



}

template <typename T, typename RepType, bool REFCNT>
struct type_traits<cnfg::sparse_matrix_default<T, RepType, REFCNT> > :
    public type_traits_default<cnfg::sparse_matrix_default<T, RepType, REFCNT> >
{
    static const bool is_specialized = true;
    typedef type_category::configurator category_type;
    static const category_type category_value;
};

namespace spct
{

/// Default spectator for sparse_matrix class.
struct sparse_matrix_idler
{};

}


template <>
struct type_traits<spct::sparse_matrix_idler> :
    public type_traits_default<spct::sparse_matrix_idler>
{
    static const bool is_specialized = true;
    typedef type_category::spectator category_type;
    static const category_type category_value;
};


namespace _Internal
{

template <typename T, typename Param, typename Category>
struct sparse_matrix_param_extractor_helper_1
{
};


/// Gives Param as configurator and set the default spectator.
template <typename T, typename Param>
struct sparse_matrix_param_extractor_helper_1<T, Param, type_category::configurator>
{
    /// Choosen configurator.
    typedef Param configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};


/// Gives Param as spectator and set the default configurator.
template <typename T, typename Param>
struct sparse_matrix_param_extractor_helper_1<T, Param, type_category::spectator>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T> configurator;

    /// Choosen spectator.
    typedef Param spectator;
};


/// Uses Param as a tag to choose the representation type.
template <typename T, typename Param>
struct sparse_matrix_param_extractor_helper_1<T, Param, type_category::tag>
{
    /// Default configurator with choosen type of representation.
    typedef cnfg::sparse_matrix_default<T, Param> configurator;

    /// Choosen spectator.
    typedef spct::sparse_matrix_idler spectator;
};


template <typename T, typename Param>
struct sparse_matrix_param_extractor
{
    /// Determined configurator.
    typedef typename sparse_matrix_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::configurator
        configurator;

    /// Determined spectator.
    typedef typename sparse_matrix_param_extractor_helper_1
    <
        T,
        Param,
        typename type_traits<Param>::type_category
    >::spectator
        spectator;
};


template <typename T>
struct sparse_matrix_param_extractor<T, default_tag>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, vecpair_col_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,vecpair_col_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, pairvec_col_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,pairvec_col_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, mappair_col_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,mappair_col_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, vecpair_row_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,vecpair_row_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, pairvec_row_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,pairvec_row_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

template <typename T>
struct sparse_matrix_param_extractor<T, mappair_row_t>
{
    /// Default configurator.
    typedef cnfg::sparse_matrix_default<T,mappair_row_t> configurator;

    /// Default spectator.
    typedef spct::sparse_matrix_idler spectator;
};

struct spmat_row_form_t
{
};

struct spmat_col_form_t
{
};

template <typename LineParam = vecpair_t, typename Form = spmat_row_form_t >
struct sparse_matrix_param_constructor
{
};

template <>
struct sparse_matrix_param_constructor<vecpair_t, spmat_row_form_t>
{
    typedef vecpair_row_t param;
};

template <>
struct sparse_matrix_param_constructor<pairvec_t, spmat_row_form_t>
{
    typedef pairvec_row_t param;
};

template <>
struct sparse_matrix_param_constructor<mappair_t, spmat_row_form_t>
{
    typedef mappair_row_t param;
};

template <>
struct sparse_matrix_param_constructor<vecpair_t, spmat_col_form_t>
{
    typedef vecpair_col_t param;
};

template <>
struct sparse_matrix_param_constructor<pairvec_t, spmat_col_form_t>
{
    typedef pairvec_col_t param;
};

template <>
struct sparse_matrix_param_constructor<mappair_t, spmat_col_form_t>
{
    typedef mappair_col_t param;
};

}


/// Sparse matrix with elements of type T.
/** Way of implementation and other configuration parameters are defined
    by the second template parameter (Param). */
template
<
    typename T,
    typename Param = default_tag
>
class sparse_matrix
{
    template<typename T1, typename Param1>
    friend class sparse_matrix;

    /// Representation of sparse matrix
    typedef typename _Internal::sparse_matrix_param_extractor<T, Param>::configurator::representation Rep;
    /// Representation form of sparse matrix (g.e. vector of sparse rows, vector of sparse columns)
    typedef typename Rep::Form Form;
public:
    /// Elements of matrix
    typedef typename T value_type;
    /// Representation of sparse vector (sparse row, sparse column)
    typedef typename Rep::RepLine LineType;
    /// Configuration parameter
    typedef typename Param Param;
    /// Line configuration parameter
    typedef typename _Internal::sparse_matrix_param_extractor<T, Param>::configurator::line_representation_tag LineParam;
    /// Reference counting property.
    /**    If it is true then counting is enabled else it is disabled. */
    static const bool refcounting = _Internal::sparse_matrix_param_extractor<T, Param>::configurator::refcounting;
public:
    /// The constructor
    sparse_matrix(size_t row=0, size_t col=0);

    /// The constructor of copy
    template<typename T1, typename Param1>
    sparse_matrix(const sparse_matrix<T1,Param1> &sm);

    /// The type cast constructor
    sparse_matrix(const Arageli::matrix<T> &mt);

    /// Sets col like sparse vector
    template<typename T1, typename Param1>
    void assign_row(size_t i, const sparse_vector<T1,Param1> &sv);

    /// Sets col like sparse vector
    template<typename T1, typename Param1>
    void assign_col(size_t j, const sparse_vector<T1,Param1> &sv);

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /// Iterators section
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular column;
    /** There is information about indices of the elements. If indices aren't important,
        use col_raw_iterator instead. */
    typedef typename Rep::col_indraw_iterator col_indraw_iterator;

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular column;
    /** There is not information about indices of the elements. If indices are important,
        use col_indraw_iterator instead. */
    typedef typename Rep::col_raw_iterator col_raw_iterator;

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular row;
    /** There is information about indices of the elements. If indices aren't important,
        use row_raw_iterator instead. */
    typedef typename Rep::row_indraw_iterator row_indraw_iterator;

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular row;
    /** There is not information about indices of the elements. If indices are important,
        use col_indraw_iterator instead. If only indices are important,use row_ind_iterator instead.*/
    typedef typename Rep::row_raw_iterator row_raw_iterator;

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements;provides information about indices.
    /** There is information about indices of the elements. If indices aren't important,
        use raw_iterator instead. If only indices are important,use ind_iterator instead.*/
    typedef typename Rep::indraw_iterator indraw_iterator;

    /// Non-constant raw iterator for held elements of the matrix only. Iterates among all stored elements for the matrix;
    /** There is not information about indices of the elements. If indices are important,
        use indraw_iterator instead. If only indices are important,use ind_iterator instead.*/
    typedef typename Rep::raw_iterator raw_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular column;
    /** There is no information about indices of the elements. If indices aren't important,
        use col_raw_iterator instead. */
    typedef typename Rep::const_col_indraw_iterator const_col_indraw_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular column;
    /** There is not information about indices of the elements. If indices are important,
        use col_indraw_iterator instead. */
    typedef typename Rep::const_col_raw_iterator const_col_raw_iterator;

    /// Constant index iterator for held indeces of elements of the matrix only. Iterates among all stored indeces of elements for a particular column;
    /** There is not information about values of the elements. If values are important,
        use col_indraw_iterator instead. */
    typedef typename Rep::const_col_ind_iterator const_col_ind_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular row;
    /** There is information about indices of the elements. If indices aren't important,
        use row_raw_iterator instead. */
    typedef typename Rep::const_row_indraw_iterator const_row_indraw_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements for a particular row;
    /** There is not information about indices of the elements. If indices are important,
        use col_indraw_iterator instead. If only indices are important,use row_ind_iterator instead.*/
    typedef typename Rep::const_row_raw_iterator const_row_raw_iterator;

    /// Constant index iterator for held indeces of elements of the matrix only. Iterates among all stored indeces of elements for a particular row;
    /** There is not information about values of the elements. If values are important,
        use col_indraw_iterator instead. If only values of the elements are important,use col_raw_iterator instead*/
    typedef typename Rep::const_row_ind_iterator const_row_ind_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements;provides information about indices.
    /** There is information about indices of the elements. If indices aren't important,
        use raw_iterator instead. If only indices are important,use ind_iterator instead.*/
    typedef typename Rep::const_indraw_iterator const_indraw_iterator;

    /// Constant raw iterator for held elements of the matrix only. Iterates among all stored elements for the matrix;
    /** There is not information about indices of the elements. If indices are important,
        use indraw_iterator instead. If only indices are important,use ind_iterator instead.*/
    typedef typename Rep::const_raw_iterator const_raw_iterator;

    /// Constant index iterator for held indeces of elements of the matrix only. Iterates among all stored indeces of elements for the matrix;
    /** There is not information about values of the elements. If values are important,
        use indraw_iterator instead. If only values of the elements are important,use raw_iterator instead*/
    typedef typename Rep::const_ind_iterator const_ind_iterator;

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /// Iterators initializations section
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Gets row`s iterator for begin row
    row_indraw_iterator begin_row(size_t row)
    {
        return rep.begin_row(row);
    }

    /// Gets constant row`s iterator for begin row
    const_row_indraw_iterator begin_row(size_t row) const
    {
        return rep.begin_row(row);
    }

    /// Gets row`s iterator for end row
    row_indraw_iterator end_row(size_t row)
    {
        return rep.end_row(row);
    }

    /// Gets constant row`s iterator for end row
    const_row_indraw_iterator end_row(size_t row) const
    {
        return rep.end_row(row);
    }

    /// Gets col`s iterator for begin col
    col_indraw_iterator begin_col(size_t col)
    {
        return rep.begin_col(col);
    }

    /// Gets constant col`s iterator for begin col
    const_col_indraw_iterator begin_col(size_t col) const
    {
        return rep.begin_col(col);
    }

    /// Gets constant col`s iterator for end col
    col_indraw_iterator end_col(size_t col)
    {
        return rep.end_col(col);
    }

    /// Gets constant col`s iterator for end col
    const_col_indraw_iterator end_col(size_t col) const
    {
        return rep.end_col(col);
    }

    /// Gets iterator for begin matrix
    indraw_iterator begin()
    {
        return rep.begin();
    }

    /// Gets constant iterator for begin matrix
    const_indraw_iterator begin() const
    {
        return rep.begin();
    }

    /// Gets iterator for end matrix
    indraw_iterator end()
    {
        return rep.end();
    }

    /// Gets constant iterator for end matrix
    const_indraw_iterator end() const
    {
        return rep.end();
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Resizes a matrix
    void resize(size_t rows, size_t cols);

    /** Reserves memory for the chosen line
        (row, column depend on matrix form) */
    void reserve_line_mem(size_t line_ind, size_t elem_count);

    /// Reserves memory for all lines (rows, columns depend on matrix form)
    void reserve_mem(size_t &lines_size);

    /// Reserves memory for the each line (row, column depend on matrix form)
    void reserve_mem(Arageli::vector<size_t> &elem_count);

    /** Resizes internal representation of the chosen line
        (row, column depend on matrix form) */
    void resize_mem_struct_line(size_t line_ind, size_t elem_count);

    /** Resizes internal representation of all lines
        (rows, columns depend on matrix form) */
    void resize_mem_struct(size_t &elem_count);

    /** Resizes internal representation of the each line
        (row, column depend on matrix form) */
    void resize_mem_struct(Arageli::vector<size_t> &elem_count);

    /** Orders internal representation of the each matrix line
        (row, column depend on matrix form) */
    void regularize();

    bool is_ordered() const
    {
        return rep.is_ordered();
    }

    /// Clears matrix
    void clear()
    {
        rep.clear();
    }

    /// Gets rows count
    size_t nrows() const
    {
        return rep.nrows();
    }

    /// Gets columns count
    size_t ncols() const
    {
        return rep.ncols();
    }

    /// Gets non-zero elements count in the matrix
    const size_t real_size() const
    {
        return rep.real_size();
    }

    /// Gets non-zero elements count in the matrix
    const size_t capacity() const
    {
        return rep.capacity();
    }

    /// Inserts elements into the i-th the matrix row and j-th the matrix column
    void insert(size_t i, size_t j, const T &t)
    {
        rep.insert(i, j, t);
    }

    /** Sets element (index in the line, value) into the i-th the matrix row
        and j-th the matrix column according to its real number in a line
        (row, column depend on matrix form) */
    void set_elem(size_t real_ind, size_t ind, const T & t)
    {
        rep.set_elem(real_ind, ind, t);
    }

    /** Puts element (index in the line, value) into the end
        of the i-th matrix line (row, column depend on matrix form) */
    void line_push_back(size_t line, size_t ind, const T & t)
    {
        rep.line_push_back(line, ind, t);
    }

    /// Sets row like sparse vector
    template<typename T1, typename Param1>
    void set_row(size_t row_in, const sparse_vector<T1,Param1> &row)
    {
        rep.set_row(row_in, row);
    }

    /// Sets row like dense vector
    template<typename T1>
    void set_row(size_t row_in, const Arageli::vector<T1> &row)
    {
        rep.set_row(row_in, row);
    }

    /// Sets column like sparse vector
    template<typename T1, typename Param1>
    void set_col(size_t col_in, const sparse_vector<T1,Param1> &col)
    {
        rep.set_col(col_in, col);
    }

    /// Sets column like dense vector
    template<typename T1>
    void set_col(size_t col_in, const Arageli::vector<T1> &col)
    {
        rep.set_col(col_in, col);
    }

    /// Sets line (row, column depend on matrix form) like sparse vector
    template<typename T1, typename Param1>
    void set_line(size_t line_in, const sparse_vector<T1,Param1> &line)
    {
        rep.set_line(line_in, line);
    }

    /// Sets line (row, column depend on matrix form) like dense vector
    template<typename T1>
    void set_line(size_t line_in, const Arageli::vector<T1> &line)
    {
        rep.set_line(line_in, line);
    }

    /// Gets line (row, column depend on matrix form) like sparse vector
    const sparse_vector<T,LineParam> & get_line(size_t line_in) const
    {
        return rep.get_line(line_in);
    }

    /// Do matrix represent as vector of sparse rows?
    bool rowform() const
    {
        return rep.rowform();
    }

    /// Do matrix represent as vector of sparse columns?
    bool colform() const
    {
        return rep.colform();
    }

    /// Gets form of the matrix representation
    Form rep_form() const
    {
        return rep.rep_form();
    }

    /// Is matrix representation dynamic?
    bool dynamic() const
    {
        return rep.dynamic();
    }

    /// Is matrix empty?
    bool is_empty() const
    {
        return rep.is_empty();
    }

    /// Swaps matrices
    void swap(sparse_matrix<T,Param> &x)
    {
        rep.swap(x.rep);
    }

    /// Copies the row
    template<typename V>
    V& copy_row(size_t i, V& res) const;

    /// Copies chosen rows
    template <typename SV, typename M>
    M& copy_rows(const SV& sv, M& res) const;

    /// Copies the column
    template<typename V>
    V& copy_col(size_t i, V& res) const;

    /// Copies chosen columns
    template <typename SV, typename M>
    M& copy_cols(const SV& sv, M& res) const;

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// The operator of type reduction
    const T & operator()(size_t i,size_t j) const;

    // WARNING!!! do not change template for second matrix
    // Use function copy_to instead
    sparse_matrix<T,Param> & operator= (const sparse_matrix<T,Param> &sm);

    sparse_matrix<T,Param> & operator= (const Arageli::matrix<T> &mt);

    sparse_matrix<T,Param> & operator= (const char *ch);

    template<typename T1, typename Param1>
    sparse_matrix<T,Param> & operator+=(const sparse_matrix<T1,Param1> &sm);

    template<typename T1, typename Param1>
    sparse_matrix<T,Param> & operator-=(const sparse_matrix<T1,Param1> &sm);

    sparse_matrix<T,Param> & operator*=(const T &val);

    sparse_matrix<T,Param> & operator/=(const T &val);

    operator matrix<T>();

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Multiplies row i by value x.
    template <typename T1>
    void mult_row (size_t i, const T1& x);

    /// Multiplies col i by value x.
    template <typename T1>
    void mult_col (size_t i, const T1& x);

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool operator==(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend bool operator!=(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_matrix<T1,Param1> operator+ (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_matrix<T1,Param1> operator- (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const big_int &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const big_int &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const float &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const float &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const double &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const double &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const big_float &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const big_float &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const rational<big_int> &val);

    template<typename T, typename Param>
    friend sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const rational<big_int> &val);

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_vector<T2,Param2> operator* (const sparse_matrix<T1,Param1> &sm, const sparse_vector<T2,Param2> &sv);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_vector<T2,Param2> operator* (const sparse_vector<T2,Param2> &sv, const sparse_matrix<T1,Param1> &sm);

    template<typename T1, typename Param1, typename T2, typename Param2>
    friend sparse_matrix<T1,Param1> operator* (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template<typename T, typename T1, typename Param1>
    friend Arageli::vector<T> operator* (const sparse_matrix<T1,Param1> &sm, const Arageli::vector<T> &vec);

    template<typename T, typename T1, typename Param1>
    friend Arageli::vector<T> operator* (const Arageli::vector<T> &vec, const sparse_matrix<T1,Param1> &sm);

    template<typename T,typename Param>
    friend std::ostream & operator <<
    (
        std::ostream & s,
        const sparse_matrix<T,Param> & x
    );

    template<typename T,typename Param>
    friend std::istream & operator >>
    (
        std::istream & s,
        sparse_matrix<T,Param> & x
    );

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

protected:

    Rep rep;

public:
    // Friend functions
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void copy_to(sparse_matrix<T1,Param1> &dst, const sparse_matrix<T2,Param2> &src);

    template <typename T,typename T1,typename Param1>
    friend void copy_to(sparse_matrix<T1,Param1> &dst, const Arageli::matrix<T> &src);

    template <typename T,typename T1,typename Param1>
    friend void copy_to(Arageli::matrix<T> &dst, const sparse_matrix<T1,Param1> &src);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend bool equal(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend bool non_equal(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void transpose(sparse_matrix<T1,Param1> &res, const sparse_matrix<T2,Param2> &sm);

    template <typename T,typename Param>
    friend sparse_matrix<T,Param> transpose(const sparse_matrix<T,Param> &sm);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void add(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void sub(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T,typename T1,typename Param1>
    friend void mul(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm, const T &val);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void mul(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void mul(sparse_vector<T2,Param2> &res, const sparse_matrix<T1,Param1> &sm, const sparse_vector<T2,Param2> &sv);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void mul(sparse_vector<T2,Param2> &res, const sparse_vector<T2,Param2> &sv, const sparse_matrix<T1,Param1> &sm);

    template <typename T1,typename Param1,typename T2,typename T3>
    friend void mul(Arageli::vector<T2> &res, const sparse_matrix<T1,Param1> &sm, const Arageli::vector<T3> &sv);

    template <typename T1,typename Param1,typename T2,typename T3>
    friend void mul(Arageli::vector<T3> &res, const Arageli::vector<T2> &sv, const sparse_matrix<T1,Param1> &sm);

    template <typename T,typename T1,typename Param1>
    friend void div(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm, const T &val);

    template <typename T1,typename Param1,typename T2,typename Param2>
    friend void add_in(sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T1,typename Param1,typename T2,typename Param2>
    void sub_in(sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2);

    template <typename T,typename T1,typename Param1>
    friend void mul_in(sparse_matrix<T1,Param1> &sm, const T &val);

    template <typename T,typename T1,typename Param1>
    friend void div_in(sparse_matrix<T1,Param1> &sm, const T &val);

};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
extern const char* sparse_matrix_output_list_first_bracket_default;
extern const char* sparse_matrix_output_list_second_bracket_default;
extern const char* sparse_matrix_output_list_row_separator_default;
extern const char* sparse_matrix_output_list_first_row_bracket_default;
extern const char* sparse_matrix_output_list_second_row_bracket_default;
extern const char* sparse_matrix_output_list_col_separator_default;
extern const char* sparse_matrix_input_list_first_bracket_default;
extern const char* sparse_matrix_input_list_second_bracket_default;
extern const char* sparse_matrix_input_list_row_separator_default;
extern const char* sparse_matrix_input_list_first_row_bracket_default;
extern const char* sparse_matrix_input_list_second_row_bracket_default;
extern const char* sparse_matrix_input_list_col_separator_default;
extern const char* sparse_matrix_output_aligned_left_col_default;
extern const char* sparse_matrix_output_aligned_right_col_default;
extern const char* sparse_matrix_output_aligned_inter_col_default;


/// Simple output of the matrix.
/** @param out an output stream
    @param x the matrix for outputting */
template <typename T, typename Param>
std::ostream& output_list
(
    std::ostream& out,
    const sparse_matrix<T, Param>& x,
    const char* first_bracket = sparse_matrix_output_list_first_bracket_default,
    const char* second_bracket = sparse_matrix_output_list_second_bracket_default,
    const char* row_separator = sparse_matrix_output_list_row_separator_default,
    const char* first_row_bracket = sparse_matrix_output_list_first_row_bracket_default,
    const char* second_row_bracket = sparse_matrix_output_list_second_row_bracket_default,
    const char* col_separator = sparse_matrix_output_list_col_separator_default
);


/// Simple input of the matrix.
template <typename T, typename Param>
std::istream& input_list
(
    std::istream& out,
    sparse_matrix<T, Param>& x,
    const char* first_bracket = matrix_input_list_first_bracket_default,
    const char* second_bracket = matrix_input_list_second_bracket_default,
    const char* row_separator = matrix_input_list_row_separator_default,
    const char* first_row_bracket = matrix_input_list_first_row_bracket_default,
    const char* second_row_bracket = matrix_input_list_second_row_bracket_default,
    const char* col_separator = matrix_input_list_col_separator_default
);


///// Aligned output of the matrix.
///// It is not realised yet
//template <typename T, typename Param>
//std::ostream& output_aligned
//(
//    std::ostream& out,
//    const sparse_matrix<T, Param>& x,
//    const char* left_col = sparse_matrix_output_aligned_left_col_default,
//    const char* right_col = sparse_matrix_output_aligned_right_col_default,
//    const char* inter_col = sparse_matrix_output_aligned_inter_col_default
//);

}// namespace Arageli


#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_sparse_matrix
    #include "sparse_matrix.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_sparse_matrix
#endif

#endif    // #ifndef _ARAGELI_sparse_matrix_hpp_
