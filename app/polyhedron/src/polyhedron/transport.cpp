/*****************************************************************************

    transport.cpp

    This file is a part of Polyhedron Software, a generator of various classes
    of polyhedra.

    The Polyhedron Software is a part of the Arageli library.

    Copyright (C) 2010 Sergey S. Lyalin

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
#include "transport.hpp"
#include "utility.hpp"


namespace
{

typedef Arageli::vector<int> task_def;
typedef Arageli::matrix<int> A;
using std::size_t;
using namespace Arageli;

size_t matrix_cols (task_def task)
{
    //cout << product(task) << '\n';
    return product(task);
}

size_t matrix_rows (task_def task)
{
    //cout << matrix_cols(task)/task << '\n';
    return sum(matrix_cols(task)/task);
}

size_t variable (task_def index, task_def task)
{
    size_t res = 0;
    size_t stride = 1;

    for(int i = 0; i < task.size(); ++i)
    {
        res += stride*index[i];
        stride *= task[i];
    }

    return res;
}

bool next (task_def& index, task_def local)
{
    //cout << index << '\n';
    for(int i = 0; i < local.size(); ++i)
        if(++index[i] >= local[i])
            index[i] = 0;
        else
            return true;

    return false;
}

A form_matrix (task_def task)
{
    A result(matrix_rows(task), matrix_cols(task));

    size_t cur_row = 0;

    for(size_t b = 0; b < task.size(); ++b)
    {
        task_def local = task;
        local.erase(b);
        task_def index(local.size(), 0);

        for(;;)
        {
            index.insert(size_t(b), 0);
            for(int i = 0; i < task[b]; ++i)
            {
                index[b] = i;
                // OK, index is index of 1 in the matrix
                // Make linear address and put 1 in result at that position

                result(cur_row, variable(index, task)) = 1;
            }
            index.erase(size_t(b));

            ++cur_row;
            if(!next(index, local))
                break;
            //cout << index << "\n";
        }
    }

    return result;
}

}


namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Transport::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    using namespace Arageli;
    using Arageli::vector;

    vector<int> vn = cmdargs.vecnumber.getValue();
    if(vn.size() < 2)
    {
        std::cerr << "[ ERROR ] You should provide at least two parameters for transport problem with key --vecnumber.";
        return;
    }

    for(size_t i = 0; i < vn.size(); ++i)
    {
        if(vn[i] < 2)
        {
            std::cerr << "[ ERROR ] You should provide parameters for transport problem which are greater than 1 with key --vecnumber.";
            return;
        }
    }

    //std::cout << vn << "\n";

    output_matrix(out, transpose(form_matrix(vn)));
}


}
}
}
