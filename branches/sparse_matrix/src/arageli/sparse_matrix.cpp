/*****************************************************************************

    sparse_matrix.cpp

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
    \file sparse_matrix.cpp
    \brief The sparse_matrix.hpp file stuff implementation.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "config.hpp"

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_sparse_matrix)

// REFERENCE ADDITIONAL HEADERS HERE

#include "sparse_matrix.hpp"


namespace Arageli
{

namespace spmt_rep
{


// Class representation of sparse matrix as vector of columns functions implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, typename ColRep, bool REFCNT>
sparsemat_colform_t<T,ColRep,REFCNT>::sparsemat_colform_t(size_t nrows, size_t ncols)
    :
    n(nrows),
    m(ncols)
{
    rep.assign(ncols,RepLine(nrows));
}

template <typename T, typename ColRep, bool REFCNT>
size_t sparsemat_colform_t<T,ColRep,REFCNT>::real_size() const
{
    size_t real_size=0;
    for(size_t i=0;i<ncols();++i)
    {
        real_size+=rep[i].real_size();
    }
    return real_size;
}

template <typename T, typename ColRep, bool REFCNT>
size_t sparsemat_colform_t<T,ColRep,REFCNT>::capacity() const
{
    size_t capacity=0;
    for(size_t i=0;i<ncols();++i)
    {
        capacity+=rep[i].real_size();
    }
    return capacity;
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::resize(size_t rows, size_t cols)
{
    n = rows;
    m = cols;

    rep.resize(cols);

    for(size_t i=0;i<cols;++i)
    {
        rep[i].resize(rows);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::reserve_line_mem(size_t line_ind, size_t line_size)
{
    rep[line_ind].reserve_mem(line_size);
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::reserve_mem(size_t &lines_size)
{
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].reserve_mem(lines_size);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::reserve_mem(Arageli::vector<size_t> &lines_size)
{
    ARAGELI_ASSERT_0(lines_size.size()==ncols());
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].reserve_mem(lines_size[i]);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::resize_mem_struct_line(size_t line_ind, size_t line_size)
{
    rep[line_ind].resize_mem_struct(line_size);
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::resize_mem_struct(size_t &lines_size)
{
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].resize_mem_struct(lines_size);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::resize_mem_struct(Arageli::vector<size_t> &lines_size)
{
    ARAGELI_ASSERT_0(lines_size.size()==ncols());
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].resize_mem_struct(lines_size[i]);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::resize_struct(size_t _m, size_t _n)
{
    rep.resize(m);
    m=_m;
    n=_n;
    for(size_t i=0;i<m;++i)
    {
        rep[i].resize(n);
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::regularize()
{
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].regularize();
    }
}

template <typename T, typename ColRep, bool REFCNT>
void sparsemat_colform_t<T,ColRep,REFCNT>::clear()
{
    for(size_t i=0;i<ncols();++i)
    {
        rep[i].clear();
    }
}

template <typename T, typename ColRep, bool REFCNT>
template<typename inT, typename inParam>
void sparsemat_colform_t<T,ColRep,REFCNT>::set_row(size_t row_in, const sparse_vector<inT,inParam> &row)
{
    ARAGELI_ASSERT_0(row_in<nrows());

    sparse_vector<T,ColRep>::indraw_iterator it;
    sparse_vector<T,ColRep>::indraw_iterator it_end;

    it = row.begin();
    it_end = row.end();

    for(;it!=it_end;++it)
    {
        rep[it.ind()].insert(row_in,it.el());
    }
}

template <typename T, typename ColRep, bool REFCNT>
template<typename inT>
void sparsemat_colform_t<T,ColRep,REFCNT>::set_row(size_t row_in, const Arageli::vector<inT> &row)
{
    ARAGELI_ASSERT_0(row_in<nrows());

    for(size_t i=0;i<ncols();++i)
    {
        rep[i].insert(row_in,row[i]);
    }
}

template <typename T, typename ColRep, bool REFCNT>
template<typename inT, typename inParam>
void sparsemat_colform_t<T,ColRep,REFCNT>::set_col(size_t col_in, const sparse_vector<inT,inParam> &col)
{
    ARAGELI_ASSERT_0(col_in<ncols());
    rep[col_in] = col;
}

template <typename T, typename ColRep, bool REFCNT>
template<typename inT>
void sparsemat_colform_t<T,ColRep,REFCNT>::set_col(size_t col_in, const Arageli::vector<inT> &col)
{
    ARAGELI_ASSERT_0(col_in<ncols());
    rep[col_in] = col;
}

template <typename T, typename ColRep, bool REFCNT>
bool sparsemat_colform_t<T,ColRep,REFCNT>::is_ordered() const
{
    for(size_t i=0;i<ncols();++i)
    {
        if(!rep[i].is_ordered())
        {
            return false;
        }
    }
    return true;
}

// Class representation of sparse matrix as vector of rows functions implementation
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename T, typename RowRep, bool REFCNT>
sparsemat_rowform_t<T,RowRep,REFCNT>::sparsemat_rowform_t(size_t rows, size_t cols)
    :
    m(rows),
    n(cols)
{
    rep.assign(m,RepLine(n));
}

template <typename T, typename RowRep, bool REFCNT>
size_t sparsemat_rowform_t<T,RowRep,REFCNT>::real_size() const
{
    size_t real_size=0;
    for(size_t i=0;i<nrows();++i)
    {
        real_size+=rep[i].real_size();
    }
    return real_size;
}

template <typename T, typename RowRep, bool REFCNT>
size_t sparsemat_rowform_t<T,RowRep,REFCNT>::capacity() const
{
    size_t capacity=0;
    for(size_t i=0;i<nrows();++i)
    {
        capacity+=rep[i].real_size();
    }
    return capacity;
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::resize(size_t rows, size_t cols)
{
    m = rows;
    n = cols;

    // WARNING vector do not initialize elements!!!
    rep.resize(rows);

    for(size_t i=0;i<rows;++i)
    {
        rep[i].clear();// for element initialization
        rep[i].resize(cols);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::reserve_line_mem(size_t line_ind, size_t line_size)
{
    rep[line_ind].reserve_mem(line_size);
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::reserve_mem(size_t &lines_size)
{
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].reserve_mem(lines_size);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::reserve_mem(Arageli::vector<size_t> &lines_size)
{
    ARAGELI_ASSERT_0(lines_size.size()==nrows());
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].reserve_mem(lines_size[i]);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::resize_mem_struct_line(size_t line_ind, size_t line_size)
{
    rep[line_ind].resize_mem_struct(line_size);
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::resize_mem_struct(size_t &lines_size)
{
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].resize_mem_struct(lines_size);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::resize_mem_struct(Arageli::vector<size_t> &lines_size)
{
    ARAGELI_ASSERT_0(lines_size.size()==nrows());
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].resize_mem_struct(lines_size[i]);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::resize_struct(size_t _m, size_t _n)
{
    // WARNING vector do not initialize elements!!!
    rep.resize(m);
    m=_m;
    n=_n;
    for(size_t i=0;i<m;++i)
    {
        rep[i].clear();// for element initialization
        rep[i].resize(n);
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::regularize()
{
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].regularize();
    }
}

template <typename T, typename RowRep, bool REFCNT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::clear()
{
    for(size_t i=0;i<nrows();++i)
    {
        rep[i].clear();
    }
}

template <typename T, typename RowRep, bool REFCNT>
template<typename inT, typename inParam>
void sparsemat_rowform_t<T,RowRep,REFCNT>::set_col(size_t col_in, const sparse_vector<inT,inParam> &col)
{
    ARAGELI_ASSERT_0(col_in<ncols());

    sparse_vector<T,RowRep>::indraw_iterator it;
    sparse_vector<T,RowRep>::indraw_iterator it_end;

    it = col.begin();
    it_end = col.end();

    for(;it!=it_end;++it)
    {
        rep[it.ind()].insert(col_in,it.el());
    }
}

template <typename T, typename RowRep, bool REFCNT>
template<typename inT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::set_col(size_t col_in, const Arageli::vector<inT> &col)
{
    ARAGELI_ASSERT_0(col_in<ncols());

    for(size_t i=0;i<nrows();++i)
    {
        rep[i].insert(col_in,col[i]);
    }
}

template <typename T, typename RowRep, bool REFCNT>
template<typename inT, typename inParam>
void sparsemat_rowform_t<T,RowRep,REFCNT>::set_row(size_t row_in, const sparse_vector<inT,inParam> &row)
{
    ARAGELI_ASSERT_0(row_in<nrows());
    rep[row_in] = row;
}

template <typename T, typename RowRep, bool REFCNT>
template<typename inT>
void sparsemat_rowform_t<T,RowRep,REFCNT>::set_row(size_t row_in, const Arageli::vector<inT> &row)
{
    ARAGELI_ASSERT_0(row_in<nrows());
    rep[row_in] = row;
}

template <typename T, typename RowRep, bool REFCNT>
bool sparsemat_rowform_t<T,RowRep,REFCNT>::is_ordered() const
{
    for(size_t i=0;i<nrows();++i)
    {
        if(!rep[i].is_ordered())
        {
            return false;
        }
    }
    return true;
}


}// namespace spmt_rep

namespace spmt_comp
{
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// All sparse matrix representation template function
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <typename A, typename B>
void copy_to(A &res, const B &sm)
{
    if(res.rep_form()==sm.rep_form())
    {
        res.rep.resize(sm.m);
        for(size_t i=0;i<sm.m;++i)
        {
            //copy_to(res.rep[i],sm.rep[i]);
            res.rep[i] = sm.rep[i];
        }
        res.m = sm.m;
        res.n = sm.n;
    }
    else
    {
        turn_over_struct(res,sm);
    }
}

template <typename T,typename A>
void copy_to(A &dst, const Arageli::matrix<T> &src)
{
    dst.clear();
    dst.resize(src.nrows(),src.ncols());
    if(dst.rowform())
    {
        A::Rep::iterator it = dst.rep.begin();
        A::Rep::iterator it_end = dst.rep.end();
        for(size_t i=0;i<dst.m;++i, ++it)
        {
            for(size_t j=0;j<dst.n;++j)
            {
                if(!is_null(src(i,j)))
                {
                    it->push_back(j,src(i,j));
                }
            }
        }
    }
    else
    {
        A::Rep::iterator it_end = dst.rep.end();
        A::Rep::iterator it = dst.rep.begin();
        for(size_t j=0;j<dst.m;++j, ++it)
        {
            for(size_t i=0;i<dst.n;++i)
            {
                if(!is_null(src(i,j)))
                {
                    it->push_back(i,src(i,j));
                }
            }
        }
    }
}

template <typename T, typename A>
void copy_to(matrix<T> &dst, const A &src)
{
    dst.assign(src.nrows(), src.ncols(), null<T>());

    A::const_indraw_iterator it = src.begin();
    A::const_indraw_iterator it_end = src.end();

    for( ; it!=it_end; ++it)
    {
        std::cout << "mt(" << it.ind_row() << ", " << it.ind_col() << ") = " << it.el() << std::endl;
        dst(it.ind_row(), it.ind_col()) = it.el();
    }
}

template <typename A,typename B>
void extract_accross_line(A &dst, const B &src, size_t line_num)
{
    dst.clear();
    dst.resize(src.m);

    typedef typename A::value_type value_type;

    if(!dst.dynamic())
    {
        std::vector<size_t> indeces;
        std::vector<const value_type*> pointers;
        B::RepLine::const_indraw_iterator it;
        B::RepLine::const_indraw_iterator it_end;
        for(size_t i=0; i<src.m; ++i)
        {
            it = src.rep[i].begin();
            it_end = src.rep[i].end();

            it = src.rep[i].find(line_num);
            if(it!=it_end)
            {
                indeces.push_back(i);
                pointers.push_back( &( it.el() ) );
            }
        }
        dst.reserve_mem(indeces.size());
        for(size_t i=0; i<indeces.size(); ++i)
        {
            dst.push_back( indeces[i], *(pointers[i]) );
        }
    }
    else
    {
        B::RepLine::const_indraw_iterator it;
        B::RepLine::const_indraw_iterator it_end;
        for(size_t i=0; i<src.m; ++i)
        {
            it = src.rep[i].begin();
            it_end = src.rep[i].end();

            it = src.rep[i].find(line_num);
            if(it!=it_end)
            {
                dst.push_back(i,it.el());
            }
        }
    }
}

template <typename A, typename B, typename SV>
void copy_lines(A &dst, const B &src, SV ind)
{
    ARAGELI_ASSERT_0(src.m>0);
    dst.clear();
    dst.rep.resize(ind.size());
    dst.m = ind.size();
    dst.n = src.n;
    for(size_t i=0; i<ind.size(); ++i)
    {
        dst.rep[i] = src.rep[ind[i]];
    }
}

template <typename A, typename B, typename SV>
void extract_accross_lines(A &dst, const B &src, SV ind)
{
    ARAGELI_ASSERT_0(src.m>0);
    dst.clear();
    dst.rep.resize(ind.size());
    dst.m = ind.size();
    dst.n = src.m;
    for(size_t i=0; i<ind.size(); ++i)
    {
        extract_accross_line(dst.rep[i], src, ind[i]);
    }
}

template <typename A, typename B>
bool equal(const A &sm1, const B &sm2)
{
    if(sm1.rep_form()!=sm2.rep_form())
    {
        A temp;
        turn_over_struct(temp,sm2);

        for(size_t i=0;i<temp.m;++i)
        {
            if(sm1.rep[i]!=temp.rep[i])
            {
                return false;
            }
        }
    }
    for(size_t i=0;i<sm2.m;++i)
    {
        if(sm1.rep[i]!=sm2.rep[i])
        {
            return false;
        }
    }
    return true;
}

template <typename A, typename B>
bool non_equal(const A &sm1, const B &sm2)
{
    return !equal(sm1, sm2);
}

template <typename A, typename B>
void turn_over_struct(A &res, const B &sm)
{
    res.clear();
    if(sm.m)
    {
        res.m = sm.n;
        res.n = sm.m;
        res.resize_struct(sm.n,sm.m);
        if(!sm.dynamic())
        {
            Arageli::vector<size_t> ip;
            B::RepLine::const_ind_iterator it;
            B::RepLine::const_ind_iterator it_end;
            ip.assign(sm.n,0);
            for(size_t i=0;i<sm.rep.size();++i)
            {
                it = sm.rep[i].begin();
                it_end = sm.rep[i].end();
                for(;it!=it_end;++it)
                {
                    ++ip[it.ind()];
                }
            }
            res.reserve_mem(ip);
        }

        B::RepLine::const_indraw_iterator it;
        B::RepLine::const_indraw_iterator it_end;

        for(size_t i=0;i<sm.m;++i)
        {
            it = sm.rep[i].begin();
            it_end = sm.rep[i].end();
            for(;it!=it_end;++it)
            {
                res.rep[it.ind()].push_back(i,it.el());
            }
        }
    }
}

template <typename A, typename B>
void transpose(A &res, const B &sm)
{
    if(res.rep_form()!=sm.rep_form())
    {
        for(size_t i=0;i<sm.rep.size();++i)
        {
            //copy_to(res.rep[i], sm.rep[i]);
            res.rep[i]=sm.rep[i];

        }
        res.m = sm.m;
        res.n = sm.n;
    }
    else
    {
        turn_over_struct(res,sm);
    }
}

template <typename A, typename B>
void add(A &res, const A &sm1, const B &sm2)
{
    ARAGELI_ASSERT_0(sm1.nrows()==sm2.nrows()&& sm1.ncols()==sm2.ncols());
    res.clear();
    res.resize(sm1.nrows(),sm1.ncols());
    if(sm1.rep_form()!=sm2.rep_form())
    {
        A temp;
        turn_over_struct(temp,sm2);
        for(size_t i=0;i<res.m;++i)
        {
            add(res.rep[i],sm1.rep[i],temp.rep[i]);
        }
    }
    else
    {
        for(size_t i=0;i<res.m;++i)
        {
            add(res.rep[i],sm1.rep[i],sm2.rep[i]);
        }
    }

}

template <typename A, typename B>
void sub(A &res, const A &sm1, const B &sm2)
{
    ARAGELI_ASSERT_0(sm1.nrows()==sm2.nrows()&& sm1.ncols()==sm2.ncols());
    res.clear();
    res.resize(sm1.nrows(),sm1.ncols());
    if(sm1.rep_form()!=sm2.rep_form())
    {
        A temp;
        turn_over_struct(temp,sm2);
        for(size_t i=0;i<res.m;++i)
        {
            sub(res.rep[i],sm1.rep[i],temp.rep[i]);
        }
    }
    else
    {
        for(size_t i=0;i<res.m;++i)
        {
            sub(res.rep[i],sm1.rep[i],sm2.rep[i]);
        }
    }
}

template <typename T, typename A>
void mul(A &res, const A &sm, const T &val)
{
    res.clear();
    res.resize(sm.nrows(),sm.ncols());
    for(size_t i=0;i<res.m;++i)
    {
        mul(res.rep[i],sm.rep[i],val);
    }
}

template <typename T, typename A>
void div(A &res, const A &sm, const T &val)
{
    res.clear();
    res.resize(sm.nrows(),sm.ncols());
    for(size_t i=0;i<res.m;++i)
    {
        div(res.rep[i],sm.rep[i],val);
    }
}

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_struct(sparse_vector<T1,Param1> &res, const A &sm, const sparse_vector<T2,Param2> &sv)
{
    ARAGELI_ASSERT_0(sm.n==sv.size());
    res.clear();
    res.resize(sm.m);
    if(!res.dynamic())
    {
        size_t count=0;
        Arageli::vector<size_t> ip(sm.m);
        for(size_t i=0;i<sm.m;++i)
        {
            if(!struct_disjoint(sm.rep[i],sv))
            {
                ip[count++]=i;
            }
        }
        res.reserve_mem(count);
        for(size_t i=0;i<count;++i)
        {
            res.push_back(ip[i],sm.rep[ip[i]]*sv);
        }
    }
    else
    {
        T1 val;
        for(size_t i=0;i<sm.m;++i)
        {
            val = sm.rep[i]*sv;
            if(!is_null(val))
            {
                res.push_back(i,val);
            }
        }
    }
}

template <typename A, typename T1, typename T2>
void mul_struct(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec)
{
    ARAGELI_ASSERT_0(sm.n==vec.size());

    res.resize(sm.m);
	T1 *pRes = &res[0];
	A::Rep::const_iterator itr = sm.rep.begin();
	A::Rep::const_iterator itr_end = sm.rep.end();;
    for( ;itr!=itr_end; ++itr)
    {
        *pRes++ = (*itr)*vec;
    }
}

template <typename A, typename T1, typename T2>
void mul_struct_accross(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec)
{
    ARAGELI_ASSERT_0(sm.m==vec.size());

    res.assign(sm.n, null<T1>());
    A::RepLine::const_indraw_iterator it;
    A::RepLine::const_indraw_iterator it_end;
    for(size_t i=0;i<sm.m; ++i)
    {
		if(!is_null(vec[i]))
		{
			it = sm.rep[i].begin();
			it_end = sm.rep[i].end();
			for(; it!=it_end; ++it)
			{
				res[it.ind()] += it.el()*vec[i];
			}
		}
    }
}

template <typename A, typename B>
void mul_matrices(A &res, const A &sm1, const B &sm2)
{
    ARAGELI_ASSERT_0(sm1.ncols()==sm2.nrows());
    res.clear();
    res.resize(sm1.nrows(),sm2.ncols());
    if(sm1.rowform()&& sm2.rowform())
    {
        B temp;
        transpose(temp,sm2);
        for(size_t i=0;i<temp.m;++i)
        {
            mul_struct(res.rep[i],temp,sm1.rep[i]);
        }
    }
    else if(sm1.colform()&& sm2.rowform())
    {
        A temp1;
        transpose(temp1,sm1);
        B temp2;
        transpose(temp2,sm2);
        A temp3(sm1.nrows(),sm2.ncols());
        for(size_t i=0;i<temp2.m;++i)
        {
            mul_struct(temp3.rep[i],temp2,temp1.rep[i]);
        }
        transpose(res,temp3);
    }
    else if(sm1.rowform()&& sm2.colform())
    {
        B temp;
        temp.resize(sm2.ncols(),sm1.nrows());
        for(size_t i=0;i<sm1.m;++i)
        {
            mul_struct(temp.rep[i],sm2,sm1.rep[i]);
        }
        transpose(res,temp);
    }
    else if(sm1.colform()&& sm2.colform())
    {
        A temp;
        transpose(temp,sm1);
        for(size_t i=0;i<sm2.m;++i)
        {
            mul_struct(res.rep[i],temp,sm2.rep[i]);
        }
    }
}

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_right(sparse_vector<T1,Param1> &res, const A &sm, const sparse_vector<T2,Param2> &sv)
{
    ARAGELI_ASSERT_0(sm.ncols()==sv.size());
    if(sm.rowform())
    {
        mul_struct(res,sm,sv);
    }
    else
    {
        A temp;
        transpose(temp,sm);
        mul_struct(res,temp,sv);
    }
}

template <typename A, typename T1, typename Param1, typename T2, typename Param2>
void mul_left(sparse_vector<T1,Param1> &res, const sparse_vector<T2,Param2> &sv, const A &sm)
{
    ARAGELI_ASSERT_0(sv.size()==sm.nrows());
    if(!sm.rowform())
    {
        mul_struct(res,sm,sv);
    }
    else
    {
        A temp;
        transpose(temp,sm);
        mul_struct(res,temp,sv);
    }
}

template <typename A, typename T1, typename T2>
void mul_right(Arageli::vector<T1> &res, const A &sm, const Arageli::vector<T2> &vec)
{
    ARAGELI_ASSERT_0(sm.ncols()==vec.size());
    if(sm.rowform())
    {
        mul_struct(res,sm,vec);
    }
    else
    {
        mul_struct_accross(res,sm,vec);
    }
}

template <typename A, typename T1, typename T2>
void mul_left(Arageli::vector<T1> &res, const Arageli::vector<T2> &vec, const A &sm)
{
    ARAGELI_ASSERT_0(vec.size()==sm.nrows());
    if(!sm.rowform())
    {
        mul_struct(res,sm,vec);
    }
    else
    {
        mul_struct_accross(res,sm,vec);
    }
}


template <typename A, typename B>
void add_in(A &sm1, const B &sm2)
{
    ARAGELI_ASSERT_0(sm1.nrows()==sm2.nrows()&& sm1.ncols()==sm2.ncols());

    if(sm1.rep_form()!=sm2.rep_form())
    {
        A temp;
        turn_over_struct(temp,sm2);
        for(size_t i=0;i<sm1.m;++i)
        {
            add_in(sm1.rep[i],temp.rep[i]);
        }
    }
    else
    {
        for(size_t i=0;i<sm1.m;++i)
        {
            add_in(sm1.rep[i],sm2.rep[i]);
        }
    }
}

template <typename A, typename B>
void sub_in(A &sm1, const B &sm2)
{
    ARAGELI_ASSERT_0(sm1.nrows()==sm2.nrows()&& sm1.ncols()==sm2.ncols());

    if(sm1.rep_form()!=sm2.rep_form())
    {
        A temp;
        turn_over_struct(temp,sm2);
        for(size_t i=0;i<sm1.m;++i)
        {
            sub_in(sm1.rep[i],temp.rep[i]);
        }
    }
    else
    {
        for(size_t i=0;i<sm1.m;++i)
        {
            sub_in(sm1.rep[i],sm2.rep[i]);
        }
    }
}

template <typename T, typename A>
void mul_in(A &sm, const T &val)
{
    for(size_t i=0;i<sm.m;++i)
    {
        mul_in(sm.rep[i],val);
    }
}

template <typename T, typename A>
void div_in(A &sm, const T &val)
{
    for(size_t i=0;i<sm.m;++i)
    {
        div_in(sm.rep[i],val);
    }
}

template <typename T, typename A>
void mul_line(A &sm, size_t line, const T &val)
{
    ARAGELI_ASSERT_0(line<sm.m);
    sm.rep[line] *= val;
}

template <typename T, typename A>
void mul_accross_line(A &sm, size_t line, const T &val)
{
    ARAGELI_ASSERT_0(line<sm.n);
    A::RepLine::indraw_iterator it;
    A::RepLine::indraw_iterator it_end;
    for(size_t i=0; i<sm.m; ++i)
    {
        it = sm.rep[i].begin();
        it_end = sm.rep[i].end();

        it = sm.rep[i].find(line);
        if(it!=it_end)
        {
            it.el() *= val;
        }
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}// namespace spmt_comp

template <typename T,typename Param>
sparse_matrix<T,Param>::sparse_matrix(size_t rows=0, size_t cols=0)
:
    rep(rows,cols)
{
}

template <typename T,typename Param>
template<typename inT, typename inParam>
sparse_matrix<T,Param>::sparse_matrix(const sparse_matrix<inT,inParam> &sm)
{
    copy_to(rep,sm.rep);
}

template <typename T,typename Param>
sparse_matrix<T,Param>::sparse_matrix(const Arageli::matrix<T> &mt)
{
    copy_to(rep,mt);
}

template <typename T,typename Param>
template<typename T1, typename Param1>
void sparse_matrix<T,Param>::assign_row(size_t i, const sparse_vector<T1,Param1> &sv)
{
    ARAGELI_ASSERT_0(i<nrows());
    ARAGELI_ASSERT_0(sv.size()==ncols());
    if(rowform())
    {
        rep.rep[i] = sv;
    }
    else
    {
        sparse_vector<T1,Param1>::const_indraw_iterator it = sv.begin();
        sparse_vector<T1,Param1>::const_indraw_iterator it_end = sv.end();
        for(; it!=it_end; ++it)
        {
            rep.rep[it.ind()].insert(i, it.el());
        }
    }
}

template <typename T,typename Param>
template<typename T1, typename Param1>
void sparse_matrix<T,Param>::assign_col(size_t j, const sparse_vector<T1,Param1> &sv)
{
    ARAGELI_ASSERT_0(j<ncols());
    ARAGELI_ASSERT_0(sv.size()==nrows());
    if(colform())
    {
        rep.rep[j] = sv;
    }
    else
    {
        sparse_vector<T1,Param1>::const_indraw_iterator it = sv.begin();
        sparse_vector<T1,Param1>::const_indraw_iterator it_end = sv.end();
        for(; it!=it_end; ++it)
        {
            rep.rep[it.ind()].insert(j, it.el());
        }
    }
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::resize(size_t rows, size_t cols)
{
    rep.resize(rows,cols);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::reserve_line_mem(size_t line_ind, size_t line_size)
{
    rep.reserve_line_mem(line_ind,line_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::reserve_mem(size_t &lines_size)
{
    rep.reserve_mem(lines_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::reserve_mem(Arageli::vector<size_t> &lines_size)
{
    rep.reserve_mem(lines_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::resize_mem_struct_line(size_t line_ind, size_t line_size)
{
    rep.resize_mem_struct_line(line_ind,line_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::resize_mem_struct(size_t &lines_size)
{
    rep.resize_mem_struct(lines_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::resize_mem_struct(Arageli::vector<size_t> &lines_size)
{
    rep.resize_mem_struct(lines_size);
}

template <typename T,typename Param>
void sparse_matrix<T,Param>::regularize()
{
    rep.regularize();
}

template <typename T,typename Param>
template<typename V>
V& sparse_matrix<T,Param>::copy_row(size_t i, V& res) const
{
    if(rowform())
    {
        res = rep.get_line(i);
    }
    else
    {
        spmt_comp::extract_accross_line(res, rep, i);
    }

    return res;
}

template <typename T,typename Param>
template<typename V>
V& sparse_matrix<T,Param>::copy_col(size_t i, V& res) const
{
    if(colform())
    {
        res = rep.get_line(i);
    }
    else
    {
        spmt_comp::extract_accross_line(res, rep, i);
    }

    return res;
}

template <typename T,typename Param>
template <typename SV, typename M>
M& sparse_matrix<T,Param>::copy_rows(const SV& sv, M& res) const
{
    if(res.rowform())
    {
        if(rowform())
        {
            spmt_comp::copy_lines(res.rep, rep, sv);
        }
        else
        {
            spmt_comp::extract_accross_lines(res.rep, rep, sv);
        }
    }
    else
    {
        M temp;
        if(rowform())
        {
            spmt_comp::copy_lines(temp.rep, rep, sv);
        }
        else
        {
            spmt_comp::extract_accross_lines(temp.rep, rep, sv);
        }
        spmt_comp::turn_over_struct(res.rep, temp.rep);
    }

    return res;
}

template <typename T,typename Param>
template <typename SV, typename M>
M& sparse_matrix<T,Param>::copy_cols(const SV& sv, M& res) const
{
    if(res.colform())
    {
        if(colform())
        {
            spmt_comp::copy_lines(res.rep, rep, sv);
        }
        else
        {
            spmt_comp::extract_accross_lines(res.rep, rep, sv);
        }
    }
    else
    {
        M temp;
        if(colform())
        {
            spmt_comp::copy_lines(temp.rep, rep, sv);
        }
        else
        {
            spmt_comp::extract_accross_lines(temp.rep, rep, sv);
        }
        spmt_comp::turn_over_struct(res.rep, temp.rep);
    }

    return res;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T,typename Param>
const T & sparse_matrix<T,Param>::operator()(size_t i,size_t j) const
{
    return rep.get(i,j);
}

template <typename T,typename Param>
//template<typename T1, typename Param1>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator= (const sparse_matrix<T,Param> &sm)
{
    spmt_comp::copy_to(rep,sm.rep);
    return *this;
}

template <typename T,typename Param>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator= (const Arageli::matrix<T> &mt)
{
    spmt_comp::copy_to(rep,mt);
    return *this;
}

template <typename T,typename Param>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator= (const char *ch)
{
    std::istringstream buf(ch);
    // WARNING. It is valid if there are no virtual function.
    static_cast<std::istream&>(buf)>> *this;
    if(!buf && !buf.eof())
        throw incorrect_string(ch);
    return *this;
}

template <typename T,typename Param>
template<typename inT, typename inParam>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator+=(const sparse_matrix<inT,inParam> &sm)
{
    spmt_comp::add_in(rep,sm.rep);
    return *this;
}

template <typename T,typename Param>
template<typename inT, typename inParam>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator-=(const sparse_matrix<inT,inParam> &sm)
{
    spmt_comp::sub_in(rep,sm.rep);
    return *this;
}

template <typename T,typename Param>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator*=(const T &val)
{
    spmt_comp::mul_in(rep,val);
    return *this;
}

template <typename T,typename Param>
sparse_matrix<T,Param> & sparse_matrix<T,Param>::operator/=(const T &val)
{
    spmt_comp::div_in(rep,val);
    return *this;
}

template <typename T,typename Param>
sparse_matrix<T,Param>::operator matrix<T>()
{
    matrix<T> mt;
    copy_to(mt, *this);
    return mt;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// Multiplies row i by value x.
template <typename T,typename Param>
template <typename T1>
void sparse_matrix<T,Param>::mult_row (size_t i, const T1& x)
{
    if(rowform())
    {
        spmt_comp::mul_line(rep, i, x);
    }
    else
    {
        spmt_comp::mul_accross_line(rep, i, x);
    }
}

/// Multiplies col i by value x.
template <typename T,typename Param>
template <typename T1>
void sparse_matrix<T,Param>::mult_col (size_t i, const T1& x)
{
    if(colform())
    {
        spmt_comp::mul_line(rep, i, x);
    }
    else
    {
        spmt_comp::mul_accross_line(rep, i, x);
    }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool operator==(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    return spmt_comp::equal(sm1.rep,sm2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline bool operator!=(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    return spmt_comp::non_equal(sm1.rep,sm2.rep);
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_matrix<T1,Param1> operator+ (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    sparse_matrix<T1,Param1> res;
    spmt_comp::add(res.rep,sm1.rep,sm2.rep);
    return res;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_matrix<T1,Param1> operator- (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    sparse_matrix<T1,Param1> res;
    spmt_comp::sub(res.rep,sm1.rep,sm2.rep);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const big_int &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::mul(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const big_int &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::div(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const float &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::mul(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const float &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::div(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const double &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::mul(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const double &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::div(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const big_float &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::mul(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const big_float &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::div(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator* (const sparse_matrix<T,Param> &sm, const rational<big_int> &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::mul(res.rep,sm.rep,val);
    return res;
}

template<typename T, typename Param>
inline sparse_matrix<T,Param> operator/ (const sparse_matrix<T,Param> &sm, const rational<big_int> &val)
{
    sparse_matrix<T,Param> res;
    spmt_comp::div(res.rep,sm.rep,val);
    return res;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_vector<T2,Param2> operator* (const sparse_matrix<T1,Param1> &sm, const sparse_vector<T2,Param2> &sv)
{
    sparse_vector<T2,Param2> res;
    spmt_comp::mul_right(res,sm.rep,sv);
    return res;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_vector<T2,Param2> operator* (const sparse_vector<T2,Param2> &sv, const sparse_matrix<T1,Param1> &sm)
{
    sparse_vector<T2,Param2> res;
    spmt_comp::mul_left(res,sv,sm.rep);
    return res;
}

template<typename T1, typename Param1, typename T2, typename Param2>
inline sparse_matrix<T1,Param1> operator* (const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    sparse_matrix<T1,Param1> res;
    spmt_comp::mul_matrices(res.rep,sm1.rep,sm2.rep);
    return res;
}

template<typename T, typename T1, typename Param1>
inline Arageli::vector<T> operator* (const sparse_matrix<T1,Param1> &sm, const Arageli::vector<T> &vec)
{
    Arageli::vector<T> res;
    spmt_comp::mul_right(res,sm.rep,vec);
    return res;
}

template<typename T, typename T1, typename Param1>
inline Arageli::vector<T> operator* (const Arageli::vector<T> &vec, const sparse_matrix<T1,Param1> &sm)
{
    Arageli::vector<T> res;
    spmt_comp::mul_left(res,vec,sm.rep);
    return res;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename T,typename Param>
std::ostream & operator <<
(
    std::ostream & s,
    const sparse_matrix<T,Param> & x
)
{
    return output_list(s,x);
}

template<typename T,typename Param>
std::istream & operator >>
(
    std::istream & s,
    sparse_matrix<T,Param> & x
)
{
    return input_list(s,x);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template <typename T,typename Param>
std::ostream& output_list
(
    std::ostream& out,
    const sparse_matrix<T,Param>& x,
    const char* first_bracket,
    const char* second_bracket,
    const char* row_separator,
    const char* first_row_bracket,
    const char* second_row_bracket,
    const char* col_separator
)
{
    out << first_bracket;

    if(!x.is_empty())
    {
        ARAGELI_ASSERT_0(x.ncols()&& x.nrows());

        sparse_matrix<T,Param> temp;

        sparse_matrix<T,Param>::const_row_indraw_iterator it, it_end;

        out << first_row_bracket;

        it = x.begin_row(0);
        it_end = x.end_row(0);

        if(it!=it_end && it.ind()==0)
        {
            out << it.el();
            ++it;
        }
        else
        {
            out << 0;
        }

        for(size_t i = 1;i < x.ncols();++i)
        {
            out << col_separator;
            if(it!=it_end && it.ind()==i)
            {
                out << it.el();
                ++it;
            }
            else
            {
                out << 0;
            }
        }

        out << second_row_bracket;

        for(size_t i = 1;i < x.nrows();++i)
        {
            it = x.begin_row(i);
            it_end = x.end_row(i);
            out << row_separator << first_row_bracket;
            if(it!=it_end && it.ind()==0)
            {
                out << it.el();
                ++it;
            }
            else
            {
                out << 0;
            }
            for(size_t j = 1;j < x.ncols();++j)
            {
                out << col_separator;
                if(it!=it_end && it.ind()==j)
                {
                    out << it.el();
                    ++it;
                }
                else
                {
                    out << 0;
                }
            }
            out << second_row_bracket;
        }
    }

    out << second_bracket;
    return out;
}


template <typename T,typename Param>
std::istream& input_list
(
    std::istream& in,
    sparse_matrix<T,Param>& x,
    const char* first_bracket,
    const char* second_bracket,
    const char* row_separator,
    const char* first_row_bracket,
    const char* second_row_bracket,
    const char* col_separator
)
{
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(first_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(second_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(row_separator));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(first_row_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(second_row_bracket));
    ARAGELI_ASSERT_0(_Internal::is_not_contains_spaces(col_separator));

    if(!_Internal::read_literal(in, first_bracket))
    {
        in.clear(std::ios_base::failbit);
        return in;
    }

    if(*second_bracket && _Internal::read_literal(in, second_bracket))
    { // empty matrix
        x.resize(0,0);
        return in;
    }

    typedef typename sparse_matrix<T,Param>::LineParam LineParam;
    typedef typename std::list<sparse_vector<T,LineParam> > Buf;// temporary buffer for rows
    Buf buf;
    std::size_t cols;

    do
    {
        sparse_vector<T,LineParam> tmp;

        input_list
        (
            in,
            tmp,
            first_row_bracket,
            second_row_bracket,
            col_separator
        );

        if(!in)
        {
            in.clear(std::ios_base::badbit);
            return in;
        }

        if(buf.empty())
        {
            cols = tmp.size();
        }
        else if(cols != tmp.size())
        {
            in.clear(std::ios_base::badbit);
            return in;
        }

        buf.push_back(tmp);

    }
    while(_Internal::read_literal(in, row_separator));


    if(!_Internal::read_literal(in, second_bracket))
    {
        in.clear(std::ios_base::badbit);
        return in;
    }

    sparse_matrix<T,Param> temp;

    if(temp.colform())
    {
        temp.resize(cols, buf.size());
    }
    else
    {
        temp.resize(buf.size(), cols);
    }

    size_t k=0;
    for(typename Buf::iterator i = buf.begin();i != buf.end();++i)
    {
        ARAGELI_ASSERT_1(i->size()== cols);
        temp.set_line(k++,*i);
    }

    if(x.colform())
    {
        transpose(x,temp);
    }
    else
    {
        x.swap(temp);
    }

    return in;
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void copy_to(sparse_matrix<T1,Param1> &dst, const sparse_matrix<T2,Param2> &src)
{
    spmt_comp::copy_to(dst.rep,src.rep);
}

template <typename T,typename T1,typename Param1>
inline void copy_to(sparse_matrix<T1,Param1> &dst, const Arageli::matrix<T> &src)
{
    spmt_comp::copy_to(dst.rep,src);
}

template <typename T,typename T1,typename Param1>
inline void copy_to(Arageli::matrix<T> &dst, const sparse_matrix<T1,Param1> &src)
{
    spmt_comp::copy_to(dst,src.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline bool equal(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    return spmt_comp::equal(sm1.rep,sm2.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline bool non_equal(const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    return spmt_comp::non_equal(sm1.rep,sm2.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void transpose(sparse_matrix<T1,Param1> &res, const sparse_matrix<T2,Param2> &sm)
{
    spmt_comp::transpose(res.rep,sm.rep);
}

template <typename T,typename Param>
inline sparse_matrix<T,Param> transpose(const sparse_matrix<T,Param> &sm)
{
    sparse_matrix<T,Param> res;
    spmt_comp::transpose(res.rep,sm.rep);
    return res;
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void add(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    spmt_comp::add(res.rep,sm1.rep,sm2.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void sub(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    spmt_comp::sub(res.rep,sm1.rep,sm2.rep);
}

template <typename T,typename T1,typename Param1>
inline void mul(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm, const T &val)
{
    spmt_comp::mul(res.rep,sm1.rep,val);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void mul(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    spmt_comp::mul_matrices(res.rep,sm1.rep,sm2.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void mul(sparse_vector<T2,Param2> &res, const sparse_matrix<T1,Param1> &sm, const sparse_vector<T2,Param2> &sv)
{
    spmt_comp::mul_right(res.rep,sm.rep,sv);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void mul(sparse_vector<T2,Param2> &res, const sparse_vector<T2,Param2> &sv, const sparse_matrix<T1,Param1> &sm)
{
    spmt_comp::mul_left(res.rep,sv,sm.rep);
}

template <typename T,typename T1,typename Param1>
inline void div(sparse_matrix<T1,Param1> &res, const sparse_matrix<T1,Param1> &sm, const T &val)
{
    spmt_comp::div(res.rep,sm1.rep,val);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void add_in(sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    spmt_comp::add_in(sm1.rep,sm2.rep);
}

template <typename T1,typename Param1,typename T2,typename Param2>
inline void sub_in(sparse_matrix<T1,Param1> &sm1, const sparse_matrix<T2,Param2> &sm2)
{
    spmt_comp::sub_in(sm1.rep,sm2.rep);
}

template <typename T,typename T1,typename Param1>
inline void mul_in(sparse_matrix<T1,Param1> &sm, const T &val)
{
    spmt_comp::sub_in(sm.rep,val);
}

template <typename T,typename T1,typename Param1>
inline void div_in(sparse_matrix<T1,Param1> &sm, const T &val)
{
    spmt_comp::sub_in(sm.rep,val);
}

template <typename T1,typename Param1,typename T2,typename T3>
inline void mul(Arageli::vector<T2> &res, const sparse_matrix<T1,Param1> &sm, const Arageli::vector<T3> &vec)
{
    spmt_comp::mul_right(res,sm.rep,vec);
}

template <typename T1,typename Param1,typename T2,typename T3>
inline void mul(Arageli::vector<T3> &res, const Arageli::vector<T2> &vec, const sparse_matrix<T1,Param1> &sm)
{
    spmt_comp::mul_left(res,vec,sm.rep);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


}// namespace Arageli


#else    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)|| ...


namespace Arageli
{

const char* sparse_matrix_output_list_first_bracket_default = "(";
const char* sparse_matrix_output_list_second_bracket_default = ")";
const char* sparse_matrix_output_list_row_separator_default = ", ";
const char* sparse_matrix_output_list_first_row_bracket_default = "(";
const char* sparse_matrix_output_list_second_row_bracket_default = ")";
const char* sparse_matrix_output_list_col_separator_default = ", ";
const char* sparse_matrix_input_list_first_bracket_default = "(";
const char* sparse_matrix_input_list_second_bracket_default = ")";
const char* sparse_matrix_input_list_row_separator_default = ",";
const char* sparse_matrix_input_list_first_row_bracket_default = "(";
const char* sparse_matrix_input_list_second_row_bracket_default = ")";
const char* sparse_matrix_input_list_col_separator_default = ",";
const char* sparse_matrix_output_aligned_left_col_default = "||";
const char* sparse_matrix_output_aligned_right_col_default = "||";
const char* sparse_matrix_output_aligned_inter_col_default = " ";

}


#endif    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE)|| ...
