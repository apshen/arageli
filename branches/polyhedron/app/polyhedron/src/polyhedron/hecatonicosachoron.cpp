/*****************************************************************************

    hecatonicosachoron.cpp

    This file is a part of Polyhedron Software, a generator of various classes
    of polyhedra.

    The Polyhedron Software is a part of the Arageli library.

    Copyright (C) 2010 Anastasya A. Ryzhova

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
#include "hecatonicosachoron.hpp"
#include "permute.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Hecatonicosachoron::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    if(cmdargs.dim.isSet() && cmdargs.dim.getValue() != 4)
    {
        throw "This type of polyhedron is supported for dimension 4 only.";
    }
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        throw "Random generation is not supported for this type";
    }
    else
    {
        unsigned int numVertices = 600;
        out << numVertices << ' ' << 5 << '\n';

        double phi = (1 + std::sqrt(5.0)) * 0.5;
        double phiSq = phi * phi;
        double phiInv = 1 / phi;
        double phiInvSq = phiInv * phiInv;

        out << std::setprecision(15);

        {/*****(0, 0, ±2, ±2) all permutations*****/
            double val1[1] = {0};
            double val2[2] = {2, -2};
            for (unsigned int i = 0; i < 1; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[i]);
                    val.push_back(val1[i]);
                    val.push_back(val2[j]);
                    val.push_back(val2[j]);
                    AllPermutations(out, val);
                    //out << '\n';
                }
            }

            {
                std::vector<double> val;
                val.reserve(4);
                val.push_back(0);
                val.push_back(0);
                val.push_back(2);
                val.push_back(-2);
                AllPermutations(out, val);
                //out << '\n';
            }
        }

        {/*****(±1, ±1, ±1, ±sqrt(5)) all permutations*****/
            double val1[2] = {1, -1};
            double val2[2] = {std::sqrt(5.0), -std::sqrt(5.0)};
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[0]);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    val.push_back(val1[1]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }
            }
        }

        {/*****(±phiInvSq, ±phi, ±phi, ±phi) all permutations*****/

            double val1[2] = {phi, -phi};
            double val2[2] = {phiInvSq, -phiInvSq};

            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val2[i]);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val2[i]);
                    val.push_back(val1[0]);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val2[i]);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    val.push_back(val1[1]);
                    AllPermutations(out, val);
                    //out << '\n';
                }
            }
        }

        {/*****(±phiInv, ±phiInv, ±phiInv, ±phiSq) all permutations*****/

            double val1[2] = {phiInv, -phiInv};
            double val2[2] = {phiSq, -phiSq};

            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    val.push_back(val1[j]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[0]);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }

                {
                    std::vector<double> val;
                    val.reserve(4);
                    val.push_back(val1[0]);
                    val.push_back(val1[1]);
                    val.push_back(val1[1]);
                    val.push_back(val2[i]);
                    AllPermutations(out, val);
                    //out << '\n';
                }
            }
        }

        {/*****(0, ±phiInvSq, ±1, ±phiSq) even permutations*****/

            double val1[1] = {0};
            double val2[2] = {phiInvSq, -phiInvSq};
            double val3[2] = {1, -1};
            double val4[2] = {phiSq, -phiSq};

            for (unsigned int i = 0; i < 1; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        for (unsigned int t = 0; t < 2; ++t)
                        {
                            out << "1 " << val1[i] << ' ' << val2[j] << ' ' << val3[k] << ' ' << val4[t] << '\n';
                            out << "1 " << val2[j] << ' ' << val3[k] << ' ' << val1[i] << ' ' << val4[t] << '\n';
                            out << "1 " << val3[k] << ' ' << val1[i] << ' ' << val2[j] << ' ' << val4[t] << '\n';
                            out << "1 " << val1[i] << ' ' << val3[k] << ' ' << val4[t] << ' ' << val2[j] << '\n';
                            out << "1 " << val2[j] << ' ' << val1[i] << ' ' << val4[t] << ' ' << val3[k] << '\n';
                            out << "1 " << val3[k] << ' ' << val2[j] << ' ' << val4[t] << ' ' << val1[i] << '\n';
                            out << "1 " << val1[i] << ' ' << val4[t] << ' ' << val2[j] << ' ' << val3[k] << '\n';
                            out << "1 " << val2[j] << ' ' << val4[t] << ' ' << val3[k] << ' ' << val1[i] << '\n';
                            out << "1 " << val3[k] << ' ' << val4[t] << ' ' << val1[i] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val3[k] << ' ' << val2[j] << ' ' << val1[i] << '\n';
                            out << "1 " << val4[t] << ' ' << val1[i] << ' ' << val3[k] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val2[j] << ' ' << val1[i] << ' ' << val3[k] << '\n';
                            //out << '\n';
                        }
                    }
                }
            }

        }

        {/*****(0, ±phiInv, ±phi, ±sqrt(5)) even permutations*****/

            double val1[1] = {0};
            double val2[2] = {phiInv, -phiInv};
            double val3[2] = {phi, -phi};
            double val4[2] = {std::sqrt(5.0), -std::sqrt(5.0)};

            for (unsigned int i = 0; i < 1; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        for (unsigned int t = 0; t < 2; ++t)
                        {
                            out << "1 " << val1[i] << ' ' << val2[j] << ' ' << val3[k] << ' ' << val4[t] << '\n';
                            out << "1 " << val2[j] << ' ' << val3[k] << ' ' << val1[i] << ' ' << val4[t] << '\n';
                            out << "1 " << val3[k] << ' ' << val1[i] << ' ' << val2[j] << ' ' << val4[t] << '\n';
                            out << "1 " << val1[i] << ' ' << val3[k] << ' ' << val4[t] << ' ' << val2[j] << '\n';
                            out << "1 " << val2[j] << ' ' << val1[i] << ' ' << val4[t] << ' ' << val3[k] << '\n';
                            out << "1 " << val3[k] << ' ' << val2[j] << ' ' << val4[t] << ' ' << val1[i] << '\n';
                            out << "1 " << val1[i] << ' ' << val4[t] << ' ' << val2[j] << ' ' << val3[k] << '\n';
                            out << "1 " << val2[j] << ' ' << val4[t] << ' ' << val3[k] << ' ' << val1[i] << '\n';
                            out << "1 " << val3[k] << ' ' << val4[t] << ' ' << val1[i] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val3[k] << ' ' << val2[j] << ' ' << val1[i] << '\n';
                            out << "1 " << val4[t] << ' ' << val1[i] << ' ' << val3[k] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val2[j] << ' ' << val1[i] << ' ' << val3[k] << '\n';
                            //out << '\n';
                        }
                    }
                }
            }

        }

        {/*****(±phiInv, ±1, ±phi, ±2) even permutations*****/

            double val1[2] = {phiInv, -phiInv};
            double val2[2] = {1, -1};
            double val3[2] = {phi, -phi};
            double val4[2] = {2, -2};

            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        for (unsigned int t = 0; t < 2; ++t)
                        {
                            out << "1 " << val1[i] << ' ' << val2[j] << ' ' << val3[k] << ' ' << val4[t] << '\n';
                            out << "1 " << val2[j] << ' ' << val3[k] << ' ' << val1[i] << ' ' << val4[t] << '\n';
                            out << "1 " << val3[k] << ' ' << val1[i] << ' ' << val2[j] << ' ' << val4[t] << '\n';
                            out << "1 " << val1[i] << ' ' << val3[k] << ' ' << val4[t] << ' ' << val2[j] << '\n';
                            out << "1 " << val2[j] << ' ' << val1[i] << ' ' << val4[t] << ' ' << val3[k] << '\n';
                            out << "1 " << val3[k] << ' ' << val2[j] << ' ' << val4[t] << ' ' << val1[i] << '\n';
                            out << "1 " << val1[i] << ' ' << val4[t] << ' ' << val2[j] << ' ' << val3[k] << '\n';
                            out << "1 " << val2[j] << ' ' << val4[t] << ' ' << val3[k] << ' ' << val1[i] << '\n';
                            out << "1 " << val3[k] << ' ' << val4[t] << ' ' << val1[i] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val3[k] << ' ' << val2[j] << ' ' << val1[i] << '\n';
                            out << "1 " << val4[t] << ' ' << val1[i] << ' ' << val3[k] << ' ' << val2[j] << '\n';
                            out << "1 " << val4[t] << ' ' << val2[j] << ' ' << val1[i] << ' ' << val3[k] << '\n';
                            //out << '\n';
                        }
                    }
                }
            }

        }

    }
}


}
}
}
