/*****************************************************************************

    compare_polyhedra.hpp

    This file is a part of phdiff tool, comparator for polyhedrons.
    Phdiff tool is a part of the Arageli library.

    Copyright (C) 2012 Anastasya A. Ryzhova

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

#ifndef _ARAGELI_APP_PHDIFF_compare_polyhedra_hpp_
#define _ARAGELI_APP_PHDIFF_compare_polyhedra_hpp_

#include "driver.hpp"
#include "kdtree.hpp"

namespace Arageli
{
namespace app
{
namespace phdiff
{

extern double threshold;

template <typename T>
bool vector_less ( Arageli::vector<T> const & vector1, Arageli::vector<T> const & vector2)
{
	Arageli::vector<T> vector;
	vector = vector1 - vector2;

	for ( int i = 0; i < vector.size(); ++i )
	{
		T eps1 = vector1[i]*threshold;
		T eps2 = vector2[i]*threshold;
		T eps = std::max(eps1, eps2);
		if ( vector[i] > eps )
		{
			return false;
		}

		if ( vector[i] < -eps )
		{
			return true;
		}
	}

	return false;
}

template <typename T>
void convert_to_vector(std::ifstream& pol, Arageli::vector<Arageli::vector<T>>& matrix, int n, int m)
{
	Arageli::vector<T> row(m);
	matrix.reserve(n);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			pol >> row[j]; 
		}

		matrix.push_back(row);
	}
}

template <typename T>
void print_matrix (std::ostream& out, const Arageli::vector<Arageli::vector<T>>& matrix)
{
	out << std::setprecision(10);
	size_t nrows = matrix.size();
	out << nrows << ' ';
	if(is_null(nrows))
	{
		out << "0\n";
		return;
	}

	size_t ncols = matrix[0].size();
	out << ncols << '\n';
	for(size_t i = 0; i < nrows; ++i)
	{
		for(size_t j = 0; j < ncols; ++j)
			out << matrix[i][j] << ' ';
		out << '\n';
	}
}

template <typename T>
void reduce_common_factor_field(Arageli::vector<Arageli::vector<T>>& matrix)
{
	int nrows = matrix.size();
	int dim = matrix[0].size();
	for ( int i = 0; i < nrows; ++i )
	{
#if 0
		for ( int j = 0; j < dim; ++j )
		{
			if ( matrix[i][j] != 0 )
			{
				T factor = abs(matrix[i][j]);
				for ( int k = j; k < dim; ++k )
				{
					matrix[i][k] /= factor;
				}
				break;
			}
		}
#else

        Arageli::vector<T> row = abs(matrix[i]);
        T factor = *std::max_element(row.begin(), row.end());
        matrix[i] /= factor;

#endif
	}
}

template <typename T>
void reduce_common_factor_integer(Arageli::vector<Arageli::vector<T>>& matrix)
{
	int nrows = matrix.size();
	int dim = matrix[0].size();
	for ( int i = 0; i < nrows; ++i )
	{
		T factor = gcd(Arageli::vector<T>(matrix[i]));

		if ( factor > 1 || factor < -1 )
		{
			for ( int j = 0; j < dim; ++j )
			{
				matrix[i][j] /= factor;
			}
		}
	}
}


template <typename T>
void transform_matrix_to_aabb(const Arageli::vector<Arageli::vector<T>>& matrix, Arageli::vector<Arageli::vector<T>>& matrix_min, Arageli::vector<Arageli::vector<T>>& matrix_max, double epsilon, bool absolute = false)
{
	size_t nrows = matrix.size();
	size_t ncolumns = matrix[0].size();

	matrix_min = matrix;
	matrix_max = matrix;

	for ( size_t i = 0; i < nrows; ++i )
	{
		for ( size_t j = 0; j < ncolumns; ++j )
		{
			T delta;
            delta = absolute ? epsilon : abs(matrix_min[i][j]*epsilon);
			matrix_min[i][j] = matrix_min[i][j] - delta;
			delta = absolute ? epsilon : abs(matrix_max[i][j]*epsilon);
			matrix_max[i][j] = matrix_max[i][j] + delta;
		}
	}

}


template <typename T>
bool traverse_matrix(node<T>* tree,  const Arageli::vector<Arageli::vector<T>>& matrix_min, const Arageli::vector<Arageli::vector<T>>& matrix_max)
{
	assert(tree);

	size_t nrows = matrix_min.size();

	Arageli::vector<node<T>*> intersect_nodes;
	intersect_nodes.reserve(nrows);

	for ( size_t i = 0; i < nrows; ++i )
	{
		node<T>* intersect_node = traverse(tree, matrix_min[i], matrix_max[i]);
		if ( intersect_node )
		{
			 intersect_nodes.push_back(intersect_node);
		}
		else
		{
            std::cout << "[ ERROR ] Cannot find a companion for " << matrix_min[i] << " - " << matrix_max[i] << "\n";
			return false;
		}
	}

	std::sort(intersect_nodes);

	node<T>* tmp_node = intersect_nodes[0];
	for ( size_t i = 1; i < nrows; ++i )
	{
		if ( intersect_nodes[i] == tmp_node )
		{
			std::cout << "Incorrect intersection!\n";
			return false;
		}
		else
		{
			tmp_node = intersect_nodes[i];
		}

	}

	return true;
}


bool compare_polyhedra(std::ifstream& pol1, std::ifstream& pol2, CmdArgs& cmdargs);


}
}
}

#endif
