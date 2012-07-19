/*****************************************************************************

    compare_polyhedra.cpp

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
#include "stdafx.hpp"

#include "driver.hpp"
#include "compare_polyhedra.hpp"

namespace Arageli
{
namespace app
{
namespace phdiff
{

double threshold = 0.0001;

bool compare_polyhedra(std::ifstream& pol1, std::ifstream& pol2, CmdArgs& cmdargs)
{
	int rows1, rows2, columns1, columns2;

	pol1 >> rows1;
	pol1 >> columns1;
	pol2 >> rows2;
	pol2 >> columns2;

	if ( rows1 != rows2 || columns1 != columns2 )
	{
		std::cerr << "Number of vertices or dimensions are different. ";
		return false;
	}

	if ( cmdargs.type.getValue() == "int" )
	{
		Arageli::vector<Arageli::vector<int>> matrix1;
		convert_to_vector(pol1, matrix1, rows1, columns1);
		reduce_common_factor_integer(matrix1);
		std::sort(matrix1);
		print_matrix(std::cout, matrix1);

		Arageli::vector<Arageli::vector<int>> matrix2;
		convert_to_vector(pol2, matrix2, rows2, columns2);
		reduce_common_factor_integer(matrix2);
		std::sort(matrix2);
		print_matrix(std::cout, matrix2);

		if ( matrix1 == matrix2 ) return true;
	}

	if ( cmdargs.type.getValue() == "bigint" )
	{
		Arageli::vector<Arageli::vector<Arageli::big_int>> matrix1;
		convert_to_vector(pol1, matrix1, rows1, columns1);
		reduce_common_factor_integer(matrix1);
		std::sort(matrix1);
		print_matrix(std::cout, matrix1);

		Arageli::vector<Arageli::vector<Arageli::big_int>> matrix2;
		convert_to_vector(pol2, matrix2, rows2, columns2);
		reduce_common_factor_integer(matrix2);
		std::sort(matrix2);
		print_matrix(std::cout, matrix2);

		if ( matrix1 == matrix2 ) return true;
	}

	if ( cmdargs.type.getValue() == "double" )
	{
		if ( cmdargs.epsilon.getValue() )
		{
			threshold = cmdargs.epsilon.getValue();
		}
		else
		{
			threshold = 0.0001;
		}


		Arageli::vector<Arageli::vector<double>> matrix1;
		convert_to_vector(pol1, matrix1, rows1, columns1);
		reduce_common_factor_field(matrix1);
		
		Arageli::vector<Arageli::vector<double>> matrix_min1, matrix_max1;
		node<double>* matrix_tree1;

        for(;;)
        {
            try
            {
		        transform_matrix_to_aabb(matrix1, matrix_min1, matrix_max1, threshold);
		        matrix_tree1 = build(matrix_min1, matrix_max1);
                break;
            }
            catch(const cannot_find_split_dim&)
            {
                double new_threshold = threshold/2;
                if(!cmdargs.epsilon_adjust.getValue() || new_threshold == threshold)
                    throw;
                std::cout
                    << "Cannot find split dim with current epsilon = " << threshold
                    << "; continue with new epsilon = " << new_threshold << "\n";
                threshold = new_threshold;
            }
        }
        
		std::cout << "Tree building is completed!\n";

		Arageli::vector<Arageli::vector<double>> matrix2;
		convert_to_vector(pol2, matrix2, rows2, columns2);
		reduce_common_factor_field(matrix2);

		Arageli::vector<Arageli::vector<double>> matrix_min2, matrix_max2;
		transform_matrix_to_aabb(matrix2, matrix_min2, matrix_max2, threshold);

		return traverse_matrix(matrix_tree1, matrix_min2, matrix_max2);

		/*Arageli::vector<Arageli::vector<double>> matrix1;
		convert_to_vector(pol1, matrix1, rows1, columns1);
		reduce_common_factor_field(matrix1);
		std::sort(matrix1, vector_less<double>);
		print_matrix(std::cout, matrix1);

		Arageli::vector<Arageli::vector<double>> matrix2;
		convert_to_vector(pol2, matrix2, rows2, columns2);
		reduce_common_factor_field(matrix2);
		std::sort(matrix2, vector_less<double>);
		print_matrix(std::cout, matrix2);

		matrix1 -= matrix2;

 		int nrows = matrix1.size();
		int dim = matrix1[0].size();
		for ( int i = 0; i < nrows; ++i )
		{
			for ( int j = 0; j < dim; ++j )
			{
				if ( matrix1[i][j] != 0 )
				{
					if ( std::abs(matrix1[i][j]) < matrix2[i][j]*threshold )
					{
						matrix1[i][j] = matrix2[i][j];
					}
					else
					{
						matrix1[i][j] += matrix2[i][j]; 
					}
				}
				else
				{
					matrix1[i][j] = matrix2[i][j];
				}
			}
		}

		if ( matrix1 == matrix2 ) return true;*/
	}

	return false;
}


}
}
}
