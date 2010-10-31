/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

// Copyright (C) 2005--2010 N.Yu.Zolotykh
//University of Nizhni Novgorod, Russia

/**
    \file
    
    Double-description method for polyhedal cone.
    (Motzkin-Raiffa-Thompson-Thrall algorithm with N.Yu.Zolotykh's variation)
 */


#ifndef DDM_HPP_
#define DDM_HPP_

#include <arageli/arageli.hpp>

#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include <limits.h>

using namespace Arageli;


#define ARITHBIGINT       0
#define ARITHINT          1
#define ARITHFLOAT        2
#define ARITHRATIONAL     3
 
#define DDM_MIN_INDEX     0
#define DDM_MAX_INDEX     1
#define DDM_LEX_MIN       2
#define DDM_LEX_MAX       3
#define DDM_MIN_CUT_OFF   4
#define DDM_MAX_CUT_OFF   5
#define DDM_RANDOM        6
#define DDM_MIN_EDGES     7
#define DDM_MAX_EDGES     8

#define INE_FACET             0
#define INE_IMPLICIT_EQUATION 1
#define INE_REDUNDANT         2

#define THNRAYS 100
#define THNRAYSPAIRS 10000000
#define THNEDGES 5000
//#define INT64 __int64
#define INT64 long long

#define WRITELOG(msg) if(logonstdout) {std::cout << msg;} if(logfileflag) {logfile << msg;}

const char* ddm_mod_names[] = {"--minindex", "--maxindex", "--lexmin", "--lexmax", "--mincutoff",
  "--maxcutoff", "--random", "--minpairs", "--maxpairs"};



template <class T>
struct ray
{
  vector<T> coordinates;
  vector<T> discrepancy;
  size_t no;
};

template <class T>
void dispose_ray_list(typename std::list<ray<T>*>& rays)
{
    for (typename std::list<ray<T>*>::iterator r = rays.begin();
      r != rays.end(); ++r)
    {
      delete *r;
    }
    rays.clear();
}

template <class T>
struct edge
{
  ray<T>* ray1;
  ray<T>* ray2;
};


// abstract
double gcd(double a, double b)
{
  //error
  return 1.0;
}

// abstract
double gcd(vector<double> a)
{
  //error
  return 1.0;
}


template <class T>
  void gauss(const matrix<T>& a, matrix<T>& f, matrix<T>& e, matrix<T>& q, size_t& r,
    vector<size_t>& perm, const int intarith, const T& eps)
{
/*Элементарными преобразованиями над строками матрицы transpose(a)
  и перестановкой ее столбцов приводит ее к диагональному виду q:
    f * transpose(a) * P = q.*/

  size_t m = a.nrows();
  size_t n = a.ncols();


  f.assign_eye(n);
  e.resize(0, n);

  q = transpose(a);

  /* main loop */
  for (size_t i = 0; i < std::min(q.ncols(), q.nrows()); )
  {
    /* find non-zero entry in the i-th row beginning from i-th entry */
    T q_pivot = std::abs(q(i, i));
    size_t j_pivot = i;
    for (size_t j = i + 1; j < m; j++)
    {
      if (std::abs(q(i, j)) > q_pivot) 
      {
         j_pivot = j;
         q_pivot = std::abs(q(i, j));
      }
    }
    if (q_pivot <= eps)
    {
      /* it's a zero row */
      q.erase_row(i);
      e.insert_row(e.nrows(), f.take_row(i));
      continue; // main loop
    }

    if (i != j_pivot)
    {
      q.swap_cols(i, j_pivot);
      perm.swap_els(i, j_pivot);
    }

    if (q(i, i) < 0)
    {
      q.mult_row(i, -1);
      f.mult_row(i, -1);
    }

    /* В i-м столбце образуем все нули */
    if (intarith)
    {
      T b = q(i, i);
      size_t mm = q.nrows();
      for (size_t ii = 0; ii < mm; ii++)
        if (ii != i)
        {
          T b_ii = q(ii, i);
          T alpha = gcd(b, b_ii);
          T b_i = b / alpha;
          b_ii = -b_ii / alpha;
          q.mult_row(ii, b_i);
          q.addmult_rows(ii, i, b_ii);
          f.mult_row(ii, b_i);
          f.addmult_rows(ii, i, b_ii);

          alpha = gcd(gcd(q.copy_row(ii)), gcd(f.copy_row(ii))); 
          q.div_row(ii, alpha);
          f.div_row(ii, alpha);
        }
    }
    else
    {
      T b = q(i, i);
      q.div_row(i, b);
      f.div_row(i, b);

      size_t mm = q.nrows();
      for (size_t ii = 0; ii < mm; ii++)
        if (ii != i)
        {
      T b_ii = -q(ii, i);
          q.addmult_rows(ii, i, b_ii);
          f.addmult_rows(ii, i, b_ii);
        }
    }
    i++; //next row
  }

  r = std::min(q.ncols(), q.nrows());
  for (size_t i = r; i < q.nrows(); i++)
    e.insert_row(e.nrows(), f.take_row(r));
}



/*
int adjacent(ray* ray1, ray* ray2, size_t r, 
  const std::list<ray*>& rays, const std::list<size_t>& inequalities)
{
  std::list<size_t> common_zeros;
  size_t common_zeros_len = 0;

  for (std::list<size_t>::iterator ine = inequalities.begin(); ine.is_correct(); ine++)
    if (ray1->discrepancy[*ine] == T(0) && ray2->discrepancy[*ine] == T(0))
    {
      common_zeros.push_back(*ine);
      common_zeros_len++;
    } 

  if (common_zeros_len < r - 2)
    return 0;

  for (std::list<ray*>::iterator cur = rays.begin(); cur.is_correct(); cur++)
  {
    if (*cur != ray1 && *cur != ray2)
    {
      int other_in_the_face = 1; // flag
      for (std::list<size_t>::iterator ine = common_zeros.begin(); ine.is_correct(); ine++)
        if ((*cur)->discrepancy[*ine] != T(0))
        {
          other_in_the_face = 0;
          break;
        }
      if (other_in_the_face)
        return 0;
    }
  }

  return 1;
}
*/

//without bit ops

int adjacent(size_t ray1, size_t ray2, size_t nrays, size_t m, size_t r, int* inc, 
  size_t iine, size_t* common_zeros)
{
  size_t common_zeros_len = 0;

  int* pray1 = inc + m * ray1;
  int* pray2 = inc + m * ray2;

  for (size_t j = 0; j <= iine; j++)
    if (pray1[j] == 0 && pray2[j] == 0)
      common_zeros[common_zeros_len++] = j;

  if (common_zeros_len < r - 2)
    return 0;

  for (size_t i = 0; i < nrays; i++)
  {
    int* pi = inc + m * i;
    if (i != ray1 && i != ray2)
    {
      for (size_t j = 0; j < common_zeros_len; j++)
        if (pi[common_zeros[j]] != 0)
          goto next_ray;

      return 0;
    }
next_ray: ;
  }

  return 1;
}


typedef unsigned long ulong;
#define SETBITS (sizeof(ulong) * CHAR_BIT)
#define MASK_ALL_ONES (~ulong(0))

template <ulong VALUE, size_t SIZE>
struct constant_filler;

template <ulong VALUE>
struct constant_filler<VALUE, size_t(0)>
{
    static const ulong value = 0;
};

template <ulong VALUE, size_t SIZE = sizeof(ulong)>
struct constant_filler
{
    static const ulong value = (constant_filler<VALUE, SIZE - 1>::value << 8) | VALUE;  // assume(number of bits in char = 8)
};

template <size_t ulong_size>
struct bit_weight_helper;
          
template <>
struct bit_weight_helper<4>
{
    inline int bit_weight(ulong b)
    {
      b = b - ((b >> 1) & constant_filler<0x55>::value);
      b = (b & constant_filler<0x33>::value) + ((b >> 2) & constant_filler<0x33>::value);
      b = (b + (b >> 4)) & constant_filler<0x0F>::value;
      b = b + (b >> 8);
      b = b + (b >> 16);
      return b & 0x03FUL;
    }
};

template <>
struct bit_weight_helper<8>
{
    inline int bit_weight(ulong b)
    {
      b = b - ((b >> 1) & constant_filler<0x55>::value);
      b = (b & constant_filler<0x33>::value) + ((b >> 2) & constant_filler<0x33>::value);
      b = (b + (b >> 4)) & constant_filler<0x0F>::value;
      b = b + (b >> 8);
      b = b + (b >> 16);
      b = b + (b >> (sizeof(ulong)*4)); // sizeof(ulong)*4 == 32, but 32 is a bad constant for 32-bit system shift (avoiding warning)
      return b & 0x07FUL;
    }
};
          
size_t set_blocks(size_t n)
{
  return (n - 1)/SETBITS + 1;
}

void set_init(ulong *set, size_t n)
{
  size_t blocks = set_blocks(n);
  set = (ulong*) calloc(blocks, sizeof(ulong));
  for (size_t i = 0; i < blocks; i++)
    set[i] = ulong(0);
}

void set_setbit(ulong *set, size_t i)
{
  size_t ii = i / SETBITS;
  size_t jj = i % SETBITS;
  ulong tmp = ulong(1) << jj;
  set[ii] |= tmp;
}

bool set_member(ulong* set, size_t i)
{
  size_t ii = i / SETBITS;
  size_t jj = i % SETBITS;

  if (set[ii] & (ulong(1) << jj))
    return true;
  
  return false;
}

size_t set_card_old(ulong* set, size_t n)
{
  size_t card = 0;

  for (size_t j = 0; j < n; j++) 
    if (set_member(set, j)) 
      card++;

  return card;
}

inline int bit_weight(ulong b)
{
  bit_weight_helper<sizeof(ulong)> helper;
  return helper.bit_weight(b);
}

inline size_t set_card(ulong* set, size_t n)
{
  size_t blocks = set_blocks(n);
  size_t card = 0;

  for (size_t j = 0; j < blocks; j++)
    card += bit_weight(set[j]);

  return card;
}

void print_bits(ulong* set, size_t n)
{
      for (size_t j = 0; j < n; j++)
      {
        if (set_member(set, j))
          std::cout << "1 ";
        else
          std::cout << "0 ";
      }
      std::cout << "\n";
}

// With bit ops

int adjacent_bit(size_t ray1, size_t ray2, size_t nrays, size_t m, size_t r, ulong* inc_bit, 
  size_t iine, ulong* uni)
{
  size_t blocks = set_blocks(iine + 1);

  ulong* pray1 = inc_bit + blocks * ray1;
  ulong* pray2 = inc_bit + blocks * ray2;

  for (size_t j = 0; j < blocks; j++)
    uni[j] = pray1[j] | pray2[j];
    //uni[j] = pray1[j] & pray2[j];

  if (set_card(uni, iine + 1) > iine + 1 - r + 2)
  {//if (set_card(uni, iine + 1) < r - 2)
    return 0;
  }

  for (size_t i = 0; i < nrays; i++)
  {
    ulong* pi = inc_bit + blocks * i;

    if (i != ray1 && i != ray2)
    {
      for (size_t j = 0; j < blocks; j++)
      {
        if ((uni[j] | ~pi[j]) != MASK_ALL_ONES)
        //if ((uni[j] & ~pi[j]) != ulong(0))
          goto next_ray;
      }

      return 0;
    }

next_ray: ;
  }

  return 1;
}

// Algebraic test - работает очень медленно и только с usual_order

template <typename T>
int adjacent_alg(const matrix<T>& a, size_t ray1, size_t ray2, size_t m, size_t n, size_t r, 
  int* inc, size_t iine, size_t* common_zeros, int usual_order)
{
  size_t common_zeros_len = 0;

  int* pray1 = inc + m * ray1;
  int* pray2 = inc + m * ray2;

  for (size_t j = 0; j <= iine; j++)
    if (pray1[j] == 0 && pray2[j] == 0)
      common_zeros[common_zeros_len++] = j;

  if (common_zeros_len < r - 2)
    return 0;

  matrix<T> suba(common_zeros_len, n);

  if (usual_order)
  {
    for (size_t i = 0; i < common_zeros_len; i++)
      for (size_t j = 0; j < n; j++)
      {
        suba(i, j) = a(common_zeros[i], j);
      }
  }

  if (rank(suba) < r - 2)
    return 0;
  else
    return 1;

}

/*

int plus_plus(size_t ray1, size_t ray2, size_t m, T* inc, size_t iine, int prefixed_order)
{
  T* pray1 = inc + m * ray1;
  T* pray2 = inc + m * ray2;

  for (size_t ine = iine + 1; ine < m; ine++)
  {
    T a = pray1[ine];
    T b = pray2[ine];
    if (prefixed_order)
    {
      if (a < T(0) && b <= T(0) || a <= T(0) && b < T(0))
        return 1;
      if (a < T(0) || b < T(0)) // a < 0 && b > 0 or viceversa
        return 0;
    }
    else
    {
      if (a < T(0) && b > T(0) || a > T(0) && b < T(0))
        return 0;
    }
  }
  return 1;
}
*/

int plus_plus(size_t ray1, size_t ray2, size_t m, int* inc, size_t iine, int prefixed_order)
{
  int* pray1 = inc + m * ray1;
  int* pray2 = inc + m * ray2;

  if (prefixed_order)
  {
    for (size_t ine = iine + 1; ine < m; ine++)
    {
      int a = pray1[ine];
      int b = pray2[ine];
      if (a < 0 && b <= 0 || a <= 0 && b < 0)
        return 1;
      if (a < 0 || b < 0) // a < 0 && b > 0 or viceversa
        return 0;
    }
  }
  else
  {
    for (size_t ine = iine + 1; ine < m; ine++)
    {
      int a = pray1[ine];
      int b = pray2[ine];
      if (a < 0 && b > 0 || a > 0 && b < 0)
        return 0;
    }
  }
  return 1;
}


template <class T>
size_t min_edges(typename std::list<size_t>::iterator ine, typename std::list<ray<T>*>& rays)
{
  size_t plus = 0;
  size_t minus = 0;
  for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
  {
    T a = (*v)->discrepancy[*ine]; 
    if (a > 0) 
      plus++;
    else if (a < 0) 
      minus++;
  }
  return -plus * minus;
}

template <class T>
size_t max_edges(typename std::list<size_t>::iterator ine, typename std::list<ray<T>*>& rays)
{
  size_t plus = 0;
  size_t minus = 0;
  for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
  {
    T a = (*v)->discrepancy[*ine]; 
    if (a > 0) 
      plus++;
    else if (a < 0) 
      minus++;
  }
  return plus * minus;
}

template <class T>
size_t max_cut_off(typename std::list<size_t>::iterator ine, typename std::list<ray<T>*>& rays)
{
  size_t minus = 0;
  for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
  {
    T a = (*v)->discrepancy[*ine]; 
    if (a < 0) 
      minus++;
  }
  return minus;
}

template <class T>
size_t min_cut_off(typename std::list<size_t>::iterator ine, typename std::list<ray<T>*>& rays)
{
  size_t plus = 0;
  for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
  {
    T a = (*v)->discrepancy[*ine]; 
    if (a > 0) 
      plus++;
  }
  return plus;
}


template <class T>
typename std::list<size_t>::iterator select_cur_ine(typename std::list<ray<T>*>& rays,
  typename std::list<size_t>& inequalities, int order)
{
  typename std::list<size_t>::iterator ine_selected = inequalities.begin();
  size_t f_selected;
  switch (order)
  {  
     case DDM_MIN_EDGES:
       f_selected = min_edges(ine_selected, rays);
       break;
     case DDM_MAX_EDGES:
       f_selected = max_edges(ine_selected, rays);
       break;
     case DDM_MIN_CUT_OFF:
       f_selected = min_cut_off(ine_selected, rays);
       break;
     case DDM_MAX_CUT_OFF:
       f_selected = max_cut_off(ine_selected, rays);
       break;
  }

  typename std::list<size_t>::iterator ine = inequalities.begin();
  ine++;
  
  for (; ine != inequalities.end(); ++ine)
  {
    size_t f;
    switch (order)
    {  
       case DDM_MIN_EDGES:
         f = min_edges(ine, rays);
         break;
       case DDM_MAX_EDGES:
         f = max_edges(ine, rays);
         break;
       case DDM_MIN_CUT_OFF:
         f = min_cut_off(ine, rays);
         break;
       case DDM_MAX_CUT_OFF:
         f = max_cut_off(ine, rays);
         break;
    }
    if (f > f_selected)
    {
      f_selected = f;
      ine_selected = ine;
    }  
  }
  return ine_selected;
}

template <class T>
void sort_inverse_order(matrix<T>& Q, vector<size_t>& perm, size_t r)
{
  size_t m = Q.ncols();
  size_t mid = (m + r)/2;
  for (size_t j = r; j < mid; j++)
  {
      Q.swap_cols(j, m - 1 - j + r);
      perm.swap_els(j, m - 1 - j + r);
  }
}

template <class T>
int lexless(const matrix<T>& Q, size_t j1, size_t j2)
{
  size_t n = Q.nrows();
    
  for (size_t ii = 0; ii < n; ii++)
  {
    if (Q(ii, j1) < Q(ii, j2))
      return 1;
    else if (Q(ii, j1) > Q(ii, j2)) 
      return 0;
  }

  return 0; // columns equal
}

template <class T>
size_t lexmin(const matrix<T>& Q, size_t j, size_t m)
// return the size_t of lexmin column in matrix Q among columns j, j + 1, ... m - 1
{
  size_t jmin = j;
  for (size_t jj = j + 1; jj < m; jj++)
  {
    if (lexless(Q, jj, jmin))
      jmin = jj;
  }
  
  return jmin;
}

// reorder columns r, r + 1, ... , m - 1 in Q

template <class T>
void sort_lex_min(matrix<T>& Q, vector<size_t>& perm, size_t r)
{
  size_t m = Q.ncols();
  for (size_t j = r; j < m; j++)
  {
    size_t jmin = lexmin(Q, j, m);
    if (j != jmin)
    {
      Q.swap_cols(j, jmin);
      perm.swap_els(j, jmin);
    }
  }
}

template <class T>
size_t lexmax(const matrix<T>& Q, size_t j, size_t m)
// return the size_t of lexmax column in matrix Q among columns j, j + 1, ... m - 1
{
  size_t jmin = j;
  for (size_t jj = j + 1; jj < m; jj++)
  {
    if (lexless(Q, jmin, jj))
      jmin = jj;
  }
  
  return jmin;
}

template <class T>
void sort_lex_max(matrix<T>& Q, vector<size_t>& perm, size_t r)  
{
  size_t m = Q.ncols();
  for (size_t j = r; j < m; j++)
  {
    size_t jmax = lexmax(Q, j, m);
    if (j != jmax)
    {
      Q.swap_cols(j, jmax);
      perm.swap_els(j, jmax);
    }
  }
}
/*
unsigned int x_random_int = 0;

unsigned int random_int()
{
  const unsigned int a = 25173;
  const unsigned int c = 13849;

  return (x_random_int = a * x_random_int + c);
}

void init_random(unsigned int x)
{
  x_random_int = x;
}

void init_random()
{
  init_random((unsigned)time(NULL));
}
*/

namespace
{
  unsigned int x_random_int = 0;
  unsigned long x_random_long = 0;
}

unsigned int random_int()
{
  const unsigned int a = 25173;
  const unsigned int c = 13849;

  return (x_random_int = a * x_random_int + c);

}

unsigned long random_long()
{
  const unsigned long a = 3141592621ul;
  const unsigned long c = 907633385ul;

  return (x_random_long = a * x_random_long + c);
}

void init_random(unsigned long x)
{
  x_random_int = (unsigned int) x;
  x_random_long = x;
}

void init_random()
{
  init_random((unsigned)time(NULL));
}

size_t random_index(size_t j, size_t m)
{
  return j + random_int() % (m - j);
}

template <class T>
void sort_random(matrix<T>& Q, vector<size_t>& perm, size_t r)  
{
  size_t m = Q.ncols();
  for (size_t j = r; j < m; j++)
  {
    size_t jnew = random_index(j, m);
    if (j != jnew)
    {
      Q.swap_cols(j, jnew);
      perm.swap_els(j, jnew);
    }
  }
}

//time_t time_pre, time_main;

template <typename T>
  int sign(const T& a, const T& eps)
{
    if (a < -eps)
      return -1;
    if (a > eps)
      return 1;
    return 0;
}

template <class T>
void check_adjacency(const matrix<T>& a, size_t m, size_t n, size_t r, 
  size_t iine, std::list<size_t>& treated_inequalities,
  std::list<size_t>& untreated_inequalities, std::list<ray<T>*>& rays, 
  std::list<edge<T> >& edges, size_t usual_order, size_t& new_edges,
  int plusplus, int prefixed_order, int graphadj, int logonstdout, int logfileflag, 
  std::ostream& logfile, const T& eps)
{
  // first of all for accelerating we store all data in an array
  size_t nrays = rays.size();

  if (nrays == 0)
  {
    new_edges = 0;
    return;
  } 

  int* inc = (int*)calloc(nrays*m, sizeof(int));
  size_t blocks = set_blocks(iine + 1);
  ulong* inc_bit = (ulong*)calloc(blocks * nrays, sizeof(ulong));
  for (size_t i = 0; i < blocks * nrays; i++)
    inc_bit[i] = ulong(0);

  ray<T>** vp = (ray<T>**)calloc(nrays, sizeof(ray<T>*));
  size_t* common_zeros = (size_t*)calloc(m, sizeof(size_t));
  ulong* uni = (ulong*)calloc(blocks, sizeof(ulong));
  size_t* adj_rays = (size_t*) calloc(nrays - 1, sizeof(size_t));

  size_t i = 0;
  size_t j = 0;

  if (usual_order)
  {
    for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
    {
      for (size_t jj = 0; jj < m; jj++)
      {
        inc[i] = sign((*v)->discrepancy[jj], eps);
        i++;
      }
      vp[j++] = *v;
    }
  }
  else
  {
    for (typename std::list<ray<T>*>::iterator v = rays.begin(); v != rays.end(); ++v)
    {
      for (typename std::list<size_t>::iterator ine = treated_inequalities.begin(); 
          ine != treated_inequalities.end(); ++ine)
        inc[i++] = sign((*v)->discrepancy[*ine], eps);
      for (typename std::list<size_t>::iterator ine = untreated_inequalities.begin(); 
          ine != untreated_inequalities.end(); ++ine)
      {
        inc[i] = sign((*v)->discrepancy[*ine], eps);
        i++;
      }
      vp[j++] = *v;
    }
  }

  for (size_t i = 0; i < nrays; i++)
  {
    ulong* p_inc_bit = inc_bit + i*blocks;
    
    size_t jj = m*i;
    for (size_t j = 0; j <= iine; j++)
    {
      if (inc[jj++] != 0)
        set_setbit(p_inc_bit, j);
    }
  }

  // second we check all pairs of extreme rays

  new_edges = 0;

  size_t criterion = iine + 1 - r + 2;

  if (graphadj)
  {
    WRITELOG("        All " << nrays << " rays being scanned" << "\n")

    for (size_t ray1 = 0; ray1 < nrays; ray1++)
    {
      if (nrays > THNRAYS)
        if (ray1 == 0 || ray1 % THNRAYS == THNRAYS - 1 || ray1 == nrays - 1)
        { 
           WRITELOG("            Ray " << ray1 + 1 << " / " << nrays
                         << " (" << (ray1 + 1) * 100 / nrays << " %)"
                         << " on iteration " << iine + 1 << " / " << m << "\n")
        }

      ulong* pray1 = inc_bit + blocks * ray1;
      size_t n_adj_rays = 0;

      //time_pre -= clock();                                              
      for (size_t ray2 = 0; ray2 < nrays; ray2++)
        if (ray2 != ray1)
        {
           ulong* pray2 = inc_bit + blocks * ray2;
           for (size_t j = 0; j < blocks; j++)
             uni[j] = pray1[j] | pray2[j];
           if (set_card(uni, iine + 1) <= criterion)
           {
             adj_rays[n_adj_rays++] = ray2;
           }
        }
      //time_pre += clock();                                              

      //time_main -= clock();                                              
      for (size_t iray2 = 0; iray2 < n_adj_rays; iray2++)
      {
        size_t ray2 = adj_rays[iray2];

        if (ray1 >= ray2)
          continue;

        if (plusplus && plus_plus(ray1, ray2, m, inc, iine, prefixed_order))
          continue;

        int adj_flag = 1;
        ulong* pray2 = inc_bit + blocks * ray2;
        for (size_t j = 0; j < blocks; j++)
          uni[j] = pray1[j] | pray2[j];
        for (size_t iray3 = 0; iray3 < n_adj_rays; iray3++)
        {
          size_t ray3 = adj_rays[iray3];

          if (ray2 != ray3)
          {
            ulong* pray3 = inc_bit + blocks * (ray3);
            for (size_t j = 0; j < blocks; j++)
            {
              if ((uni[j] | ~pray3[j]) != MASK_ALL_ONES)
              //if ((uni[j] & ~pray3[j]) != ulong(0))
                goto next_ray;
            }
            adj_flag = 0;
          }
next_ray: ;
        }
        if (adj_flag)
        {
          new_edges++;

          edge<T> e;
          e.ray1 = vp[ray1];
          e.ray2 = vp[ray2];
          edges.push_back(e);
        }
      }

      //time_main += clock();                                              
    }
  }
  else  //nographadj
  {
    INT64 npairs = ((INT64)(nrays) - 1)*nrays/2;
    INT64 ipairs = 0;

    WRITELOG("        " << nrays << " rays\n")
    WRITELOG("        All " << npairs << " pairs of rays being scanned\n")

    for (size_t ray1 = 0; ray1 < nrays; ray1++)
    {
      for (size_t ray2 = ray1 + 1; ray2 < nrays; ray2++)
      {
        ipairs = ipairs + 1;

        if (npairs > THNRAYSPAIRS)
          if (ipairs == 1 || ipairs % THNRAYSPAIRS == 0 || ipairs == npairs)
          {  
             WRITELOG("            Step " << ipairs << " / " << npairs 
                      << " (" << ipairs * 100 / npairs << " %)"
                      << " on iteration " << iine + 1 << " / " << m << "\n")
          }

        if ((!plusplus || !plus_plus(ray1, ray2, m, inc, iine, prefixed_order))
          && 
             //adjacent(ray1, ray2, nrays, m, r - 1, inc, iine, common_zeros)
             adjacent_bit(ray1, ray2, nrays, m, r - 1, inc_bit, iine, uni)
             //adjacent_alg(a, ray1, ray2, m, n, r, inc, iine, common_zeros, usual_order)
            )
        {
          new_edges++;

          edge<T> e;
          e.ray1 = vp[ray1];
          e.ray2 = vp[ray2];
          edges.push_back(e);
        }
      }
    } 
  }
  free(inc);                    
  free(inc_bit);
  free(uni);
  free(vp);
  free(common_zeros);
  free(adj_rays);
}

template <typename T>
  T firstnonzero(const vector<T>& a, const T& eps)
{
  for (size_t i = 0; i < a.length(); i++)
    if (std::abs(a[i]) > eps)
      return std::abs(a[i]);

  return eps;
}

template <typename T>
  void clearzeros(vector<T>& a, const T& eps)
{
  for (size_t i = 0; i < a.length(); i++)
    if (std::abs(a[i]) < eps)
      a[i] = T(0);
}

/***************************/
/*                         */
/*           ddm           */
/*                         */
/***************************/


template <class T>
void do_ddm(const matrix<T>& a,
  matrix<T>& v, matrix<T>& b, matrix<T>& inc, matrix<size_t>& edges_ind, 
  INT64& totalnumrays,  INT64& totalnumedges,
  int edgesflag,
  int order, int prefixed_order, int graphadj, int plusplus, int intarith, const T& eps, 
  int logonstdout, int logfileflag, std::ostream& logfile); // forward

template <class T>
void ddm(const matrix<T>& ineq, const matrix<T>& eq, 
  matrix<T>& ext, matrix<T>& bas, matrix<T>& dis, matrix<size_t>& edges_ind, 
  INT64& totalnumrays,  INT64& totalnumedges,
  int edgesflag,
  int order, int prefixed_order, int graphadj, int plusplus, int intarith, const T& eps, 
  int logonstdout, int logfileflag, std::ostream& logfile)
{
    if (eq.nrows() == 0)
    {
        do_ddm(ineq, ext, bas, dis, edges_ind, totalnumrays, totalnumedges, edgesflag,
            order, prefixed_order, graphadj, plusplus, intarith, eps, 
            logonstdout, logfileflag, logfile);
    }
    else
    {
        matrix<T> bas_sub, f, q;
        size_t m = ineq.nrows();
        vector<size_t> perm(m, fromsize);
        //perm.range(0, m - 1); 
        for (size_t j = 0; j < m; j++)
            perm[j] = j;

        size_t r;

        gauss(eq, f, bas_sub, q, r, perm, intarith, eps);

        do_ddm(ineq*transpose(bas_sub), 
        ext, bas, dis, edges_ind, totalnumrays, totalnumedges, edgesflag,
            order, prefixed_order, graphadj, plusplus, intarith, eps, 
            logonstdout, logfileflag, logfile);

        ext = ext*bas_sub;
        bas = bas*bas_sub;

        size_t ext_nrows = ext.nrows();
        size_t bas_nrows = bas.nrows();

        for (size_t ii = 0; ii < ext_nrows; ii++)
        {
            T delta;
            if (intarith)
            {
               delta = gcd(ext.copy_row(ii));
            }
            else
            {
               delta = firstnonzero(ext.copy_row(ii), eps);
            }
            ext.div_row(ii, delta);
            dis.div_row(ii, delta);
        }

        for (size_t ii = 0; ii < bas_nrows; ii++)
        {
            T delta;
            if (intarith)
            {
               delta = gcd(bas.copy_row(ii));
            }
            else
            {
               delta = firstnonzero(bas.copy_row(ii), eps);
            }
            bas.div_row(ii, delta);
        }

    }
}


template <class T>
void do_ddm(const matrix<T>& a,
  matrix<T>& v, matrix<T>& b, matrix<T>& inc, matrix<size_t>& edges_ind, 
  INT64& totalnumrays,  INT64& totalnumedges,
  int edgesflag,
  int order, int prefixed_order, int graphadj, int plusplus, int intarith, const T& eps, 
  int logonstdout, int logfileflag, std::ostream& logfile)
{

/*
std::vector<std::list<size_t> > adj; 
std::vector<std::list<size_t> > extinc;
std::vector<std::list<size_t> > ineinc; 
int extincflag = 0;
int ineincflag = 0; 
int adjacencyflag = 0;
*/

  //time_pre = time_main = 0;

  typename std::list<edge<T> > edges; // linked list of all edges
  typename std::list<ray<T>*> rays, rays_plus, rays_minus, rays_zero;
  //std::list<ray<T>*> basis;
  typename std::list<size_t> untreated_inequalities, treated_inequalities;

  size_t m = a.nrows();
  size_t n = a.ncols();
  //size_t k = 0; // number of rays

  matrix<T> f, q;
  vector<size_t> perm(m, fromsize);
  //perm.range(0, m - 1); 
  for (size_t j = 0; j < m; j++)
    perm[j] = j;
  
  size_t r; // = f.nrows(); // == n - e.nrows();

  gauss(a, f, b, q, r, perm, intarith, eps);

/*
std::cout << "\n a = \n"; matrix_output(std::cout, a);
std::cout << "\n f = \n"; matrix_output(std::cout, f);
std::cout << "\n b = \n"; matrix_output(std::cout, b);
std::cout << "\n q = \n"; matrix_output(std::cout, q);
std::cout << "\n r = \n" << r;      
std::cout << "\n perm = " << perm;
*/

  // basis = eye
/*
  for (size_t j = 0; j < n; j++)
  {
    ray<T>* v = new ray<T>;
    v->coordinates.zeros(n);
    v->coordinates[j] = T(1);
    v->discrepancy = a.copy_col(j);
    basis.push_back(v);
  }
*/

  int usual_order = 1; // 0 or 1  usual_order == 1 means that inequalities will be
                       // considered in usual order and lists of (un)treated
                       // inequalities are not used

  switch (order)
  {  
     case DDM_MIN_INDEX:
       usual_order = 1;
       break;
     case DDM_MAX_INDEX:
       usual_order = 1;
       sort_inverse_order(q, perm, r); // reorder columns r, r + 1, ... , m - 1 in Q
       break;
     case DDM_LEX_MIN:
       usual_order = 1;
       sort_lex_min(q, perm, r); 
       break;
     case DDM_LEX_MAX:     
       usual_order = 1;
       sort_lex_max(q, perm, r);  
       break;
     case DDM_RANDOM:
       usual_order = 1;
       sort_random(q, perm, r);
       break;
     case DDM_MIN_CUT_OFF: 
     case DDM_MAX_CUT_OFF: 
     case DDM_MIN_EDGES: 
     case DDM_MAX_EDGES: 
       prefixed_order = 0;
       usual_order = 0;
       break;
  }

  WRITELOG("rank = " << r << "\n")
  WRITELOG("Initial set of inequalities:\n")

  for (size_t i = 0; i < r; i++)
  {
    ray<T>* v = new ray<T>;
    v->coordinates = f.copy_row(i);
    v->discrepancy = q.copy_row(i);
    if (!intarith)
    {
       T alpha = firstnonzero(v->coordinates, eps);
       v->coordinates /= alpha;
       v->discrepancy /= alpha;
       clearzeros(v->coordinates, eps);
       clearzeros(v->discrepancy, eps);
    }
    rays.push_back(v);
    WRITELOG(perm[i] + 1 << " ");
  }
  WRITELOG("\n")

  // Its edges
  for (typename std::list<ray<T>*>::iterator ray1 = rays.begin(); ray1 != rays.end(); ++ray1)
  {
    typename std::list<ray<T>*>::iterator ray2 = ray1;
    ray2++;
    for (; ray2 != rays.end(); ++ray2)
    {
        edge<T> e;
        e.ray1 = *ray1;
        e.ray2 = *ray2;
        edges.push_back(e);
    }
  }

  totalnumrays = rays.size();
  totalnumedges = edges.size();

  // (Un)treated inequalities preliminaries
  if (!usual_order)
  {
    // all inequalities are untreated
    for (size_t i = 0; i < r; i++)
      treated_inequalities.push_back(i);
    for (size_t i = r; i < m; i++)
      untreated_inequalities.push_back(i);
  }   

  // Main loop: consider all inequalities one by one
  for (size_t iine = r; iine < m; iine++)
  {

    logfile.flush();

    WRITELOG("\n---------------------------------------------------------------- \n")
    WRITELOG("Iteration "<< iine + 1 << " / " << m)
    WRITELOG(" (Inequality No ")

    size_t cur_ine;

    if (usual_order)
    {
      cur_ine = iine;
    }
    else
    {
      //cur_ine = untreated_inequalities.pop_front();
      //cur_ine = select_cur_ine(rays, untreated_inequalities, order).pop();
      typename std::list<size_t>::iterator cur_ine_it = 
        select_cur_ine(rays, untreated_inequalities, order);
      cur_ine = *cur_ine_it;
      untreated_inequalities.erase(cur_ine_it);
      
      treated_inequalities.push_back(cur_ine);
    }

    WRITELOG(perm[cur_ine] + 1 << " in the original system)\n")

    size_t nraysplus = 0;
    size_t nraysminus = 0;
    size_t nrayszero = 0;

    WRITELOG("    Rays classifying...\n")
    while (!rays.empty())
    {
      typename std::list<ray<T>*>::iterator v = rays.begin();
      if ((*v)->discrepancy[cur_ine] > eps)
      {
        nraysplus++;
        rays_plus.push_back(rays.front());
        rays.pop_front();
      }
      else if ((*v)->discrepancy[cur_ine] < -eps)
      {
        nraysminus++;
        rays_minus.push_back(rays.front());
        rays.pop_front();
      }
      else
      {
        nrayszero++;
        rays_zero.push_back(rays.front());
        rays.pop_front();
      }
    }

    WRITELOG("        The number of rays inside the current cone =   " << nraysplus << "\n")
    WRITELOG("        The number of rays outside the current cone =  " << nraysminus << "\n")
    WRITELOG("        The number of rays in the current hyperplane = " << nrayszero << "\n")
    WRITELOG("        Total number of rays =                         " 
              << nraysplus + nraysminus + nrayszero << "\n")

    if (rays_plus.empty() && rays_zero.empty())    
    {
      // system has only zero solution
      WRITELOG("The cone is subspace\n")
      dispose_ray_list(rays_minus);
      break;
    }

    if (rays_minus.empty())
    {    
      // redundant inequality
      WRITELOG("The current inequality is redundant and it follows from previous ones\n")
      rays.splice(rays.begin(), rays_plus);
      rays.splice(rays.begin(), rays_zero);
      continue;  
    }

    size_t new_rays = 0;
    size_t iedges = 0;
    size_t nedges = edges.size(); 
    WRITELOG("    All edges (" << nedges << ") being scanned for generating new rays\n")
    for(typename std::list<edge<T> >::iterator e = edges.begin(); e != edges.end(); )
    {
      iedges++;
      if (nedges > THNEDGES)
        if (iedges == 1 || iedges % THNEDGES == 0 || iedges == nedges)
        {  
           WRITELOG("          Edge " << iedges << " / " << nedges 
                    << " (" << iedges * 100 / nedges << " %)"
                    << " on iteration " << iine + 1 << " / " << m << "\n")
        }
      ray<T>* ray1 = (*e).ray1;
      ray<T>* ray2 = (*e).ray2;
      T alpha = ray1->discrepancy[cur_ine];
      T beta = ray2->discrepancy[cur_ine];
      if (alpha < -eps && beta > eps || alpha > eps && beta < -eps)
      {
        // +- edge
        new_rays++;
        ray<T>* v = new ray<T>;
        T abs_alpha = std::abs(alpha);
        T abs_beta = std::abs(beta);
        v->coordinates = ray1->coordinates * abs_beta + ray2->coordinates * abs_alpha; 
        v->discrepancy = ray1->discrepancy * abs_beta + ray2->discrepancy * abs_alpha; 
        if (intarith)
        {
           T delta = gcd(v->coordinates);
           v->coordinates /= delta;
           v->discrepancy /= delta;
        }
        else
        {
           T delta = firstnonzero(v->coordinates, eps);
           v->coordinates /= delta;
           v->discrepancy /= delta;
           clearzeros(v->coordinates, eps);
           clearzeros(v->discrepancy, eps);
        }
        rays_zero.push_back(v);
        if (alpha < -eps)
          (*e).ray1 = v;
        else
          (*e).ray2 = v;
        e++;
      }
      else if (alpha <= eps && beta <= eps)
      {
        // --, -0, 00 edge
        e = edges.erase(e);
      }
      else // ++ edge
        ++e;
    }

    dispose_ray_list(rays_minus);

    WRITELOG("            New " << new_rays << " rays constructed\n")

    totalnumrays += new_rays;

    if (iine < m - 1)
    {
      WRITELOG("    Constructing new edges in the current hyperplane...\n")
      size_t new_edges;
      check_adjacency(a, m, n, r, iine, treated_inequalities, untreated_inequalities, rays_zero, 
        edges, usual_order, new_edges, plusplus, prefixed_order, graphadj, 
        logonstdout, logfileflag, logfile, eps);
      WRITELOG("                " << new_edges << " new edges constructed\n")
      totalnumedges += new_edges;
    }

    rays.splice(rays.begin(), rays_plus);
    rays.splice(rays.begin(), rays_zero);

  } // main loop

  size_t nrays = rays.size();
  size_t nedges = 0;

  if (edgesflag || r <= 4) //|| adjacencyflag 
  {
    edges.clear();

    WRITELOG("\n****************************************************************\n")
    WRITELOG("    Final constructing of all edges in the cone...\n")
    check_adjacency(a, m, n, r, m - 1, treated_inequalities, untreated_inequalities, rays, edges, 
        usual_order, nedges, 0, prefixed_order, graphadj, logonstdout, logfileflag, logfile, eps);
    WRITELOG("                Final number of edges = " << nedges << "\n")
    totalnumedges += nedges;
  }
  
  v.resize(nrays, n);
  inc.resize(nrays, m);

  size_t i = 0;
  for (typename std::list<ray<T>*>::iterator cur = rays.begin();
    cur != rays.end(); ++cur)
  {
    if (!intarith)
    {
       T delta = firstnonzero((*cur)->coordinates, eps); // Это съедает много времени, тем более, если поставить в цикл
       (*cur)->coordinates /= delta;                     // Но если это приводит к ошибкам, то нужно оставить там      
       (*cur)->discrepancy /= delta;                     // (см. 2 закомментированных куска
       clearzeros((*cur)->coordinates, eps);             // Да, это действительно приводит к ошибкам
       clearzeros((*cur)->discrepancy, eps);
    }
    v.assign_row(i, (*cur)->coordinates);
    //inc.assign_row(i, (*cur)->discrepancy);
    for (size_t j = 0; j < m; j++)
      inc(i, perm[j]) = (*cur)->discrepancy[j];
    (*cur)->no = i; 
    i++;
  }

  if (edgesflag)
  {
    edges_ind.resize(nedges, 2);
    i = 0;
    for (typename std::list<edge<T> >::iterator cur = edges.begin();
      cur != edges.end(); ++cur)
    {
      edges_ind(i, 0) = (*cur).ray1->no;
      edges_ind(i, 1) = (*cur).ray2->no;
      i++;
    }
  }
  
/*
  if (r <= 4) // отдельно нужно рассмотреть случай 2d
    facets3d(adj, ineinc);
*/  
  dispose_ray_list(rays);

  //WRITELOG(" time_pre = " << time_pre <<
  //         " time_main = " << time_main)
}


void post_construct_ineinc(const matrix<int>& dis_int, std::vector<std::list<size_t> >& ineinc)
{
  size_t k = dis_int.nrows();
  size_t m = dis_int.ncols();

  ineinc.resize(m);

  for (size_t i = 0; i < m; i++)
  {
    for (size_t j = 0; j < k; j++)
      if (!dis_int(j, i))
        ineinc[i].push_back(j);
  }
}

void post_construct_extinc(const matrix<int>& dis_int, std::vector<std::list<size_t> >& extinc)
{
  size_t k = dis_int.nrows();
  size_t m = dis_int.ncols();

  extinc.resize(k);

  for (size_t i = 0; i < k; i++)
  {
    for (size_t j = 0; j < m; j++)
      if (!dis_int(i, j))
        extinc[i].push_back(j);
  }
}

void post_construct_adjacency(size_t nrays, const matrix<size_t>& edges_ind, 
  std::vector<std::list<size_t> >& adj)
{
  size_t nedges = edges_ind.nrows();

  adj.resize(nrays);

  for (size_t i = 0; i < nedges; i++)
  {
    adj[edges_ind(i, 0)].push_back(edges_ind(i, 1));
    adj[edges_ind(i, 1)].push_back(edges_ind(i, 0));
  }
}

/////////////////////////


template <class T>
  void disarray2intarray(const matrix<T>& dis, matrix<int>& dis_int, const T& eps)
{
  size_t m = dis.nrows();
  size_t n = dis.ncols();

  //dis_int.clear();
  dis_int.resize(m, n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      if (dis(i, j) <= eps)
        dis_int(i, j) = 0;
      else
        dis_int(i, j) = 1;
}

template <class T>
  void disarray2intarraytranspose(const matrix<T>& dis, matrix<int>& dis_int, const T& eps)
{
  size_t m = dis.nrows();
  size_t n = dis.ncols();

  //dis_int.clear();
  dis_int.resize(n, m);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      if (dis(i, j) <= eps)
        dis_int(j, i) = 0;
      else
        dis_int(j, i) = 1;
}

void post_construct_inetype_and_ridges(
  const matrix<int>& dis_int, std::vector<int>& ine_type, matrix<size_t>& ridges)
{
  ridges.resize(0, 2);

  // first of all for accelerating we store all data in a bit array

  size_t m = dis_int.ncols();
  size_t k = dis_int.nrows();

  ine_type.clear();
  ine_type.resize(m, INE_FACET);


  size_t blocks = set_blocks(k + 1);
  ulong* dis_bit = (ulong*)calloc(blocks * m, sizeof(ulong));
  for (size_t i = 0; i < blocks*m; i++)
    dis_bit[i] = ulong(0);

  size_t jj = 0;
  for (size_t i = 0; i < m; i++)
  {
    ulong* p_dis_bit = dis_bit + i*blocks;
    
    for (size_t j = 0; j < k; j++)
    {
      if (dis_int(j, i) != 0)
        set_setbit(p_dis_bit, j);
    }
  }

  // second we find all implicit equations

  for (size_t i = 0; i < m; i++)
  {
    ulong* p_dis_bit = dis_bit + i*blocks;
    for (size_t j = 0; j < blocks; j++)
    {
      if (p_dis_bit[j] != ulong(0))
        goto next_i;
    }
    ine_type[i] = INE_IMPLICIT_EQUATION;
    next_i: ;
  }

  // 3d we find all redundant inequalities

  for (size_t i = m; i > 0; )
  {
    i--;

    ulong* p_i = dis_bit + i*blocks;

    if (ine_type[i] == INE_IMPLICIT_EQUATION) 
      continue;
    // Добавить: if num_zeros < d_U - d_L - 1 // then redundant
    for (size_t ii = 0; ii < m; ii++)
    {
      if (ii != i && ine_type[ii] == INE_FACET)
      {
        ulong* p_ii = dis_bit + ii*blocks; //int* p_ii = inc_int + ii*k;
        bool overlap_flag = true;
        for (size_t j = 0; j < blocks; j++) //for (size_t j = 0; j < k; j++)
        {
          // if (inc(j, i) <= eps && inc(j, ii) > eps)
          if (~p_i[j] & p_ii[j]) // if (!p_i[j] && p_ii[j])
          {
            overlap_flag = false;
            break;
          }
        } 
        if (overlap_flag)
        {
          ine_type[i] = INE_REDUNDANT;
          break;
        }
      }
    }
  }

  // 4th we find ridges (adjacency of facets)

  ulong* uni = (ulong*)calloc(blocks, sizeof(ulong));

  for (size_t i = 0; i < m; i++)
    if (ine_type[i] == INE_FACET)
    {
      ulong* p_i = dis_bit + i*blocks;
      for (size_t ii = i + 1; ii < m; ii++)
        if (ine_type[ii] == INE_FACET)
        {
          ulong* p_ii = dis_bit + ii*blocks; 
          for (size_t j = 0; j < blocks; j++)
            uni[j] = p_i[j] | p_ii[j];
            //uni[j] = p_i[j] & p_ii[j];
          // if (set_card(uni, iine + 1) > iine + 1 - r + 2) // ДОбавить!
          // {//if (set_card(uni, iine + 1) < r - 2)
          //   return 0;
          // }
          for (size_t iii = 0; iii < m; iii++)
          {
            if (iii != i && iii != ii && ine_type[iii] == INE_FACET)
            {
              ulong* p_iii = dis_bit + blocks*iii;
              for (size_t j = 0; j < blocks; j++)
              {
                if ((uni[j] | ~p_iii[j]) != MASK_ALL_ONES)
                //if ((uni[j] & ~pi[j]) != ulong(0))
                  goto next_iii;
              }

              goto next_pair;
            }

            next_iii: ;
          }

          //std::cout << i + 1 << " " << ii + 1 << "\n";
          ridges.insert_row(ridges.nrows());
          ridges(ridges.nrows() - 1, 0) = i;
          ridges(ridges.nrows() - 1, 1) = ii;

          next_pair: ;
        }
    }

  free(dis_bit);
  free(uni);
}

/////////////////////////

void facets3d(std::vector<std::list<size_t> >& adj, std::vector<std::list<size_t> >& ineinc)
{
  // Of course, there is much quicker algorithm. Some day it will appear here

  size_t m = ineinc.size();

  for (size_t i = 0; i < m; i++)
  {
    if (ineinc[i].size() <= 3)
      continue;

    std::list<size_t> new_inc;

    size_t first_ray = *(ineinc[i].begin());
    size_t prev_ray = first_ray;
    size_t cur_ray = first_ray;
    size_t next_ray = first_ray;

    do
    {
      new_inc.push_back(cur_ray);

      for (std::list<size_t>::iterator iadjray = adj[cur_ray].begin();
        iadjray != adj[cur_ray].end(); ++iadjray) // for all adjacent rays
      {
         std::list<size_t>::iterator iiadjray = 
           std::find(ineinc[i].begin(), ineinc[i].end(), *iadjray);
         if (iiadjray != ineinc[i].end() && *iadjray != prev_ray) // if ray is in the current facet
         {
           next_ray = *iadjray;
           break;
         }
      }

      //assert(cur_ray != next_ray);

      prev_ray = cur_ray;
      cur_ray = next_ray;

    } while(next_ray != first_ray);

    ineinc[i].swap(new_inc);
  }
}

                                                                                           
#undef WRITELOG

#endif
