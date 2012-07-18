#ifndef _ARAGELI_APP_PHDIFF_kdtree_hpp_
#define _ARAGELI_APP_PHDIFF_kdtree_hpp_


#include "stdafx.hpp"

namespace Arageli
{
namespace app
{
namespace phdiff
{

template<typename T>
struct node
{
	Arageli::vector<T> min;
	Arageli::vector<T> max;

	node* left;
	node* right;

	int split_dim;
	T split_val;
};

template<typename T>
bool find_split(node<T>* resnode, const Arageli::vector<Arageli::vector<T>>& matrix_min, const Arageli::vector<Arageli::vector<T>>& matrix_max, const Arageli::vector<size_t>& prohibited_dim)
{
	size_t nrows = matrix_min.size();
	size_t ncolumns = matrix_min[0].size();

	Arageli::vector<T> size = resnode->max - resnode->min;
	size_t prohibited_dim_num = prohibited_dim.size();

	size_t splitdim;
	T maxsize;

	// Choose start splitdim

	if ( prohibited_dim_num ) 
	{
		bool is_start_splitdim = true;

		for	( size_t i = 0; i < ncolumns; ++i )
		{
			is_start_splitdim = true;
			for ( size_t j = 0; j < prohibited_dim_num; ++j )
			{
				if ( i == prohibited_dim[j] )
				{
					is_start_splitdim = false;
					break;
				}
			}

			if ( is_start_splitdim )
			{
				splitdim = i;
				break;
			}
		}

		if ( !is_start_splitdim )
		{
			for ( size_t i = 0; i < nrows; ++i )
			{
				for ( size_t j = 0; j < ncolumns; ++j )
				{
					std::cout << matrix_min[i][j] << " ";
				}

				std::cout << '\n';

				for ( size_t j = 0; j < ncolumns; ++j )
				{
					std::cout << matrix_max[i][j] << " ";
				}

				std::cout << "\n\n";

			}
			throw "Can't find split dim.\n";
		}
	}
	else
	{
		splitdim = 0;
	}

	maxsize = size[splitdim];

	for ( size_t i = 0; i < ncolumns; ++i )
	{
		if ( size[i] > maxsize )
		{
			bool true_max = true;
			for ( size_t j = 0; j < prohibited_dim_num; ++j )
			{
				if ( i == prohibited_dim[j] )
				{
					true_max = false;
					break;
				}
			}

			if ( true_max )
			{
				maxsize = size[i];
				splitdim = i;
			}
		}
	}

	Arageli::vector<T> vector_min(nrows);
	Arageli::vector<T> vector_max(nrows);

	for ( size_t i = 0; i < nrows; ++i )
	{
		vector_min[i] = matrix_min[i][splitdim];
		vector_max[i] = matrix_max[i][splitdim];
	}

	std::sort(vector_min);
	std::sort(vector_max);
	
	T splitvalue;

	bool end = false;
	bool split_exist = false;

	size_t middle = nrows / 2;

	size_t j = 0;
    for ( size_t i = 0; i < nrows - 1; ++i )
	{
		while ( vector_min[j] <= vector_max[i] )
		{
			++j;
			if ( j == nrows )
			{
				end = true;
				break;
			}
		}

		if ( end )
		{
			break;
		}

		size_t num_max = i + 1;
		size_t num_min = j;

		if ( num_max - num_min == 0 )
		{
			split_exist = true;
			splitvalue = vector_max[i];
			if ( num_max >= middle )
			{
				break;
			}
		}
	}

	if ( split_exist )
	{
		resnode->split_val = splitvalue;
	}

	resnode->split_dim = splitdim;

	return split_exist;
}

template<typename T>
Arageli::vector<T> find_min(const Arageli::vector<Arageli::vector<T>>& matrix)
{
	size_t rows = matrix.size();
	size_t columns = matrix[0].size();

	Arageli::vector<T> vector = matrix[0];

	for ( size_t i = 1; i < rows; ++i )
	{
		for ( size_t j = 0; j < columns; ++j )
		{
			if ( matrix[i][j] < vector[j] )
			{
				vector[j] = matrix[i][j];
			}
		}
	}

	return vector;
}

template<typename T>
Arageli::vector<T> find_max(const Arageli::vector<Arageli::vector<T>>& matrix)
{
	size_t rows = matrix.size();
	size_t columns = matrix[0].size();

	Arageli::vector<T> vector = matrix[0];

	for ( size_t i = 1; i < rows; ++i )
	{
		for ( size_t j = 0; j < columns; ++j )
		{
			if ( matrix[i][j] > vector[j] )
			{
				vector[j] = matrix[i][j];
			}
		}
	}
	return vector;
}

template<typename T>
bool do_split(node<T>* tree, 
			  const Arageli::vector<Arageli::vector<T>>& matrix_min,
			  const Arageli::vector<Arageli::vector<T>>& matrix_max,
			  Arageli::vector<Arageli::vector<T>>& matrix_left_min,
			  Arageli::vector<Arageli::vector<T>>& matrix_left_max,
			  Arageli::vector<Arageli::vector<T>>& matrix_right_min,
			  Arageli::vector<Arageli::vector<T>>& matrix_right_max)

{
		size_t nrows = matrix_min.size();
		size_t ncolumns = matrix_min[0].size();

        Arageli::vector<size_t> prohibited_split_dim;
		prohibited_split_dim.reserve(ncolumns);

		while ( !find_split(tree, matrix_min, matrix_max, prohibited_split_dim) )
		{
			if ( prohibited_split_dim.size() == ncolumns )
			{
				return false;
			}

 			prohibited_split_dim.push_back(tree->split_dim);
		}
		
		for ( size_t i = 0; i < nrows; ++i )
		{
			if ( matrix_min[i][tree->split_dim] < tree->split_val )
			{
				matrix_left_min.push_back(matrix_min[i]);
				matrix_left_max.push_back(matrix_max[i]);
			}

			if ( matrix_max[i][tree->split_dim] > tree->split_val )
			{
				matrix_right_min.push_back(matrix_min[i]);
				matrix_right_max.push_back(matrix_max[i]);
			}
		}

        assert((matrix_left_min.size() < nrows)&&(matrix_right_min.size() < nrows));

		return true;
}

template<typename T>
node<T>* build(const Arageli::vector<Arageli::vector<T>>& matrix_min, const Arageli::vector<Arageli::vector<T>>& matrix_max)
{
	assert(matrix_min.size() == matrix_max.size());
	size_t rows = matrix_min.size();
	size_t ncolumns = matrix_min[0].size();

	if ( rows == 0 )
	{
		return 0;
	}

	node<T>* res = new node<T>;

	res->min = find_min(matrix_min);
	res->max = find_max(matrix_max);

	if ( rows > 1 )
	{
		Arageli::vector<Arageli::vector<T>> matrix_min_left_tmp;
		matrix_min_left_tmp.reserve(rows);
		Arageli::vector<Arageli::vector<T>> matrix_max_left_tmp;
		matrix_max_left_tmp.reserve(rows);
		Arageli::vector<Arageli::vector<T>> matrix_min_right_tmp;
		matrix_min_right_tmp.reserve(rows);
		Arageli::vector<Arageli::vector<T>> matrix_max_right_tmp;
		matrix_max_right_tmp.reserve(rows);

		if ( !do_split(	res,
						matrix_min,
						matrix_max,
						matrix_min_left_tmp,
						matrix_max_left_tmp,
						matrix_min_right_tmp,
						matrix_max_right_tmp) )
		{
            /*output_aligned(std::cerr, matrix_min);
            std::cerr << "\n";
            output_aligned(std::cerr, matrix_max);*/
			throw "Bad matrix from input1! Can't build tree!\n";
		}

		res->left = build(matrix_min_left_tmp, matrix_max_left_tmp);
		res->right = build(matrix_min_right_tmp, matrix_max_right_tmp);
	}
	else
	{
		res->left = 0;
		res->right = 0;
	}

	return res;
}

template <typename T>
bool interval_intersect(T x1_min, T x1_max, T x2_min, T x2_max)
{
	return (  x2_max >= x1_min && x1_max >= x2_min );
}

template <typename T>
bool aabb_intersect(Arageli::vector<T> vector1_min, Arageli::vector<T> vector1_max, Arageli::vector<T> vector2_min, Arageli::vector<T> vector2_max)
{
	size_t dim = vector1_min.size();
	for ( size_t i = 0; i < dim; ++i )
	{
		if ( !interval_intersect(vector1_min[i], vector1_max[i], vector2_min[i], vector2_max[i]) )
		{
			return false;
		}
	}

	return true;
}


template <typename T>
node<T>* traverse(node<T>* tree, Arageli::vector<T> vector_min, Arageli::vector<T> vector_max)
{
	if ( !tree )
	{
		return 0;
	}

	if ( !aabb_intersect(tree->min, tree->max, vector_min, vector_max) )
	{
		return 0;
	}

	if ( !tree->right )
	{
		assert(!tree->left);
		return tree;
	}

	if ( tree->split_val < vector_min[tree->split_dim] )
	{
		return traverse(tree->right, vector_min, vector_max);
	}
	else
	{ 
		if ( tree->split_val < vector_max[tree->split_dim] )
		{
			node<T>* node_left = traverse(tree->right, vector_min, vector_max);
			node<T>* node_right = traverse(tree->left, vector_min, vector_max);

			if ( node_left == node_right )
			{
				return node_left;
			}
			else
			{
				throw("Bad matrix!\n");
			}
		}
		else
		{
			return traverse(tree->left, vector_min, vector_max);
		}
	}
	
}

}
}
}

#endif