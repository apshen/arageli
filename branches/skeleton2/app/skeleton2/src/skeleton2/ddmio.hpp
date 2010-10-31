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

#ifndef DDMIO_HPP_
#define DDMIO_HPP_

#include <arageli/arageli.hpp>

#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <list>
#include <vector>
#include "ddm.hpp"

using namespace Arageli;

#define SIMPLEFORMAT      0
#define SKELETONFORMAT    1
#define AVISFUKUDAFORMAT  2

#define WRITEOUTPUT(msg)               \
  if (outputonstdout)                  \
  {                                    \
    std::cout << msg;                  \
  }                                    \
  if (outputfileflag)                  \
  {                                    \
    outputfile << msg;                 \
  }

#define WRITEOUTPUTA(msg, a)           \
  if (outputonstdout)                  \
  {                                    \
    std::cout << msg;                  \
    matrix_output(std::cout, a);       \
  }                                    \
  if (outputfileflag)                  \
  {                                    \
    outputfile << msg;                 \
    matrix_output(outputfile, a);      \
  }

#define WRITELOG(msg) if(logonstdout) {std::cout << msg;} if(logfileflag) {logfile << msg;}
#define WRITESUMMARY(msg) if(summaryonstdout) {std::cout << msg;} if(summaryfileflag) {summaryfile << msg;}

template <class T>
void print_rays(const typename std::list<ray<T>*>& rays)
{
    for (typename std::list<ray<T>*>::iterator cur = rays.begin(); cur != rays.end(); ++cur)
    {
      std::cout << (*cur)->coordinates << "\n";
      std::cout << (*cur)->discrepancy << "\n";
    }
}


void print_no(const char* title,
  const std::vector<int>& ine_type, int type, 
  int outputonstdout, int outputfileflag, std::ostream& outputfile)
{
  WRITEOUTPUT(title)
  size_t s = 0;
  size_t m = ine_type.size();
  for (size_t i = 0; i < m; i++)
    if (ine_type[i] == type)
      s++;

  WRITEOUTPUT(s << "\n")

  for (size_t i = 0; i < m; i++)
    if (ine_type[i] == type)
    {
      WRITEOUTPUT(i + 1 << " ")
    }

  WRITEOUTPUT("\n")
}

void print_adjacency(const char* title, const std::vector<std::list<size_t> >& adj,
  int outputonstdout, int outputfileflag, std::ostream& outputfile)
{
  WRITEOUTPUT(title)
  for (size_t i = 0; i < adj.size(); i++)
  {
    WRITEOUTPUT(i + 1 << ":")
    for (std::list<size_t>::const_iterator cur = adj[i].begin();
      cur != adj[i].end(); ++cur)
    {
       WRITEOUTPUT(" " << (*cur) + 1)
    }
    WRITEOUTPUT("\n")
  }
}

void print_edges(const char* title, const matrix<size_t>& edges_ind,
  int outputonstdout, int outputfileflag, std::ostream& outputfile)
{
  WRITEOUTPUT(title)
  WRITEOUTPUT(edges_ind.nrows() << "\n")
  
  for (size_t i = 0; i < edges_ind.nrows(); i++)
  {
      WRITEOUTPUT(edges_ind(i, 0) + 1 << " " << edges_ind(i, 1) + 1 << "\n")
  }
}


                                                                                           
template <typename A_type>                                                                 
void input_using_avis_fukuda_format(                                               
        std::istream& s,                                                                   
        A_type& A,                                                                         
        int logonstdout,                                                                   
        int logfileflag,                                                                 
        std::ostream& logfile)                                                          
{                                                                                        
  std::string str;                                                                       
                                                                                         
  WRITELOG("Avis-Fukuda format for input\n")

  while (true) // process options until "begin" found
  {
    if (!s.good()) 
    {
      WRITELOG("Error: Something wrong with input file!")
      return;
    }
    str = "";
    s >> str;
    if (str == "begin")
    {
      break;
    }
    if (str.size() > 0 && str[0] == '*') // skip the line beginning with *
    {
      getline(s, str);
      continue;
    }
    WRITELOG("Option = " << str)
    if (str != "H-representation" && str != "V-representation")
    {
      WRITELOG(" -> unsupported")
    }
    WRITELOG(std::endl)
  }

  matrix_input(s, A);

  //matrix_output(std::cout, A);                                                       

  s >> str; // str must be "end"

  while (s.good()) // process options until eof
  {
    str = "";
    s >> str;
    if (str.size() > 0 && str[0] == '*') // skip the line beginning with *
    {
      getline(s, str);
      continue;
    }
    if (str == "")
      continue;

    WRITELOG("Option = " << str)
    WRITELOG(" -> ignored")
    WRITELOG(std::endl)
  }

}
                                                                                           
template <typename A_type>                                                                 
void input_using_skeleton_format(                                               
        std::istream& s,                                                                   
        A_type& ine,
        A_type& equ,                                                                         
        int logonstdout,                                                                   
        int logfileflag,                                                                 
        std::ostream& logfile)                                                          
{                                                                                        
  bool equ_flag = false;                                                      
  bool ine_flag = false;                                                      
                                                                                         
  WRITELOG("SKELETON format for input\n")

  while (s.good())
  {
    char ch;
    s >> ch;
    if (ch == '*')
    {
        std::string str = "";
        s >> str;
        if (str == "Equations:")
        {
           WRITELOG("Equations found in input file\n")
           if (equ_flag)
           {
               WRITELOG("Warning: at least 2 equations field in the file\n")
           }
           matrix_input(s, equ);
           equ_flag = true;
        }
        else if (str == "Inequalities:")
        {
           WRITELOG("Inequalities found in input file\n")
           if (ine_flag)
           {
               WRITELOG("Warning: at least 2 inequalities field in the file\n")
           }
           matrix_input(s, ine);
           ine_flag = true;
        }
    }
    std::string str = "";
    getline(s, str); // skip the remaining part of the line
  }

  if (!equ_flag && !ine_flag)
  {
      WRITELOG("Warning: no equation nor inequality found\n")
  }
  else if (equ_flag && !ine_flag)
  {
      WRITELOG("No inequality in input file\n")
      ine = A_type(0, equ.ncols(), fromsize);                                             
  }
  else if (!equ_flag && ine_flag)
  {
      WRITELOG("No equation in input file\n")
      equ = A_type(0, ine.ncols(), fromsize);                                             
  }
  else if (equ.ncols() != ine.ncols())
  {
      WRITELOG("Warning: the numbers of columns in equation matrix and\n")
      WRITELOG("inequality matrix are different\n")
      WRITELOG("Equations ignored\n")
      equ = A_type(0, ine.ncols(), fromsize);                                             
  }
std::cout << "ine = " << ine << "\n";
std::cout << "equ = " << equ << "\n";
}

template <class T>
void input_matrices(matrix<T>& ine, matrix<T>& eq, 
    int inputoutputformat, int inputfromstdin, std::istream& inputfile, 
    int logonstdout, int logfileflag, std::ostream& logfile)
{ 
    if (inputoutputformat == SKELETONFORMAT)
    {
        if (inputfromstdin)                                                                          
            input_using_skeleton_format(std::cin, ine, eq, logonstdout, logfileflag, logfile);          
        else                                                                                         
            input_using_skeleton_format(inputfile, ine, eq, logonstdout, logfileflag, logfile);         
    }
    else
    {                                                                                                
        if (inputoutputformat == AVISFUKUDAFORMAT)                                                     
        {                                                                                              
            if (inputfromstdin)                                                                          
                input_using_avis_fukuda_format(std::cin, ine, logonstdout, logfileflag, logfile);          
            else                                                                                         
                input_using_avis_fukuda_format(inputfile, ine, logonstdout, logfileflag, logfile);         
        }                                                                                              
        else /* SIMPLEFORMAT */                                                                                          
        {                                                                                              
            if (inputfromstdin)                                                                          
                matrix_input(std::cin, ine);                                                               
            else                                                                                         
                matrix_input(inputfile, ine);                                                              
        }
        eq = matrix<T>(0, ine.ncols(), fromsize);                                             
    }                                                                                              
}  

template <class T>
void output_results(const matrix<T>& ine, const matrix<T>& equ, 
    const matrix<T>& ext, const matrix<T>& bas, const matrix<T>& inc,
    const matrix<size_t>& edges_ind, 
    int output_ine, int output_equ, int output_ext, int output_bas, int output_inc, int edgesflag,
    int outputonstdout, int outputfileflag, std::ostream& outputfile)
{
                                                                                             
  if (output_ine)                                                                                  
  {                                                                                                
    WRITEOUTPUTA("* Inequalities:\n", ine)                                                         
  }                                                                                                
  if (output_equ)                                                                                  
  {                                                                                                
    WRITEOUTPUTA("* Equations:\n", equ)                                                         
  }                                                                                                
  if (output_bas)                                                                                  
  {                                                                                                
    WRITEOUTPUTA("* Basis:\n", bas)                                                                
  }                                                                                                
  if (output_ext)                                                                                  
  {                                                                                                
    WRITEOUTPUTA("* Extreme rays:\n", ext)                                                         
  }                                                                                                
  if (output_inc)                                                                                  
  {                                                                                                
    WRITEOUTPUTA("* Discrepancies:\n", inc)                                                     
  }
  if (edgesflag)                                                       
  {                                                                    
    print_edges("* Edges:\n", edges_ind, outputonstdout, outputfileflag, outputfile);
  } 
}

template <class T>
void output_summary(
    const matrix<T>& ine, const matrix<T>& equ, 
    const matrix<T>& ext, const matrix<T>& bas, const matrix<T>& inc,
    INT64 totalnumrays, INT64 totalnumedges,
    int modification, int prefixed_order, int graphadj, int plusplus, int arith, const T& eps,
    int argc, char *argv[], 
    char* inputfilename, char* outputfilename, char* logfilename, char* summaryfilename, 
    int summaryonstdout, int summaryfileflag, std::ostream& summaryfile)
{                                                                                                
    WRITESUMMARY("\n\n")                                                                           
    WRITESUMMARY("----------------------\n")                                                       
    WRITESUMMARY("   Computation done   \n")                                                       
    WRITESUMMARY("----------------------\n")                                                       
    WRITESUMMARY("\n")                                                                             
    WRITESUMMARY("  Command line: skeleton")                                                       
    for (int i = 1; i < argc; i++)                                                                 
    {                                                                                              
        WRITESUMMARY(" ")                                                                          
        WRITESUMMARY(argv[i])                                                                      
    }                                                                                              
    WRITESUMMARY("\n")                                                                             
    WRITESUMMARY("\n")                                                                             
    WRITESUMMARY("  Input file:   " << inputfilename << "\n")                                      
    WRITESUMMARY("  Output file:  " << outputfilename << "\n")                                     
    WRITESUMMARY("  Log file:     " << logfilename << "\n")                                        
    WRITESUMMARY("  Summary file: " << summaryfilename << "\n")                                    
    WRITESUMMARY("\n")                                                                             
    WRITESUMMARY("  ine sizes = " << ine.nrows() << " x " << ine.ncols() << "\n")                  
    WRITESUMMARY("  equ sizes = " << equ.nrows() << " x " << equ.ncols() << "\n")                  
    WRITESUMMARY("  bas sizes = " << bas.nrows() << " x " << bas.ncols() << "\n")                  
    WRITESUMMARY("  ext sizes = " << ext.nrows() << " x " << ext.ncols() << "\n")                  
    WRITESUMMARY("  inc sizes = " << inc.nrows() << " x " << inc.ncols() << "\n")                  
    WRITESUMMARY("\n")                                                                             

    WRITESUMMARY("  Order    = " << ddm_mod_names[modification] << "\n")                        
    if (prefixed_order)                                                                            
    {                                                                                              
       WRITESUMMARY("  Prefixed = --prefixedorder" << "\n")                                        
    }                                                                                              
    else                                                                                           
    {                                                                                              
       WRITESUMMARY("  Prefixed = --noprefixedorder" << "\n")                                      
    }                                                                                              
    if (graphadj)                                                                                  
    {                                                                                              
       WRITESUMMARY("  Graphadj = --graphadj" << "\n")                                             
    }                                                                                              
    else                                                                                           
    {                                                                                              
       WRITESUMMARY("  Graphadj = --nographadj" << "\n")                                           
    }                                                                                              
    if (arith == ARITHBIGINT)                                                                      
    {                                                                                              
       WRITESUMMARY("  Arith    = --bigint" << "\n")                                               
    }                                                                                              
    else if (arith == ARITHINT)                                                                    
    {                                                                                              
       WRITESUMMARY("  Arith    = --int" << "\n")                                                  
    }                                                                                              
    else if (arith == ARITHFLOAT)                                                                  
    {                                                                                              
       WRITESUMMARY("  Arith    = --float" << "       zerotol = " << eps << "\n")                  
    }                                                                                              
    else                                                                                           
    {                                                                                              
       WRITESUMMARY("  Arith    = --rational" << "\n")                                             
    }                                                                                              
    if (plusplus)                                                                                  
    {                                                                                              
       WRITESUMMARY("  Plusplus = --plusplus" << "\n")                                             
    }                                                                                              
    else                                                                                           
    {                                                                                              
       WRITESUMMARY("  Plusplus = --noplusplus" << "\n")                                           
    }                                                                                              
    WRITESUMMARY("\n")                                                                   
    WRITESUMMARY("  Total number of rays  generated (all iterations) = " << totalnumrays << "\n")  
    WRITESUMMARY("  Total number of edges generated (all iterations) = " << totalnumedges << "\n") 
    WRITESUMMARY("\n")                                                                   
}

void output_time(time_t begin_time, time_t end_time, timer time_sec,
    int summaryonstdout, int summaryfileflag, std::ostream& summaryfile)
{          
                                                                                                   
    WRITESUMMARY("  Computation starts:     " << asctime(localtime(&begin_time)))                  
    WRITESUMMARY("  Computation terminates: " << asctime(localtime(&end_time)))
    double t = time_sec.time();                    
    WRITESUMMARY("\n  Time elapsed = " << t << " s (")                             
    WRITESUMMARY(floor(t/3600) << " h  ")                                                            
    t = fmod(t, 3600);                                                                                  
    WRITESUMMARY(floor(t/60) << " m  ")                                                              
    t = fmod(t, 60);                                                                                  
    WRITESUMMARY(t << " s)\n")                                                             
}                                                                                                   
                                                                                              

template <typename A_type>
    void matrix_input(std::istream& s, A_type& A)
{
    typedef typename A_type::difference_type index;
    typedef typename A_type::value_type T;

    index m, n;
    std::string str;

    s >> m;
    s >> n;

    getline(s, str);

    A.resize(m, n);

    for (index i = 0; i < m; i++)
        for (index j = 0; j < n; j++)
            s >> A(i, j);
}

template <typename A_type>
    void matrix_output(std::ostream& s, const A_type& A)
{
    typedef typename A_type::difference_type index;
    typedef typename A_type::value_type T;

    index m = A.nrows();
    index n = A.ncols();

    s << m << " " << n << "\n";

    for (index i = 0; i < m; i++)
    {
        for (index j = 0; j < n; j++)
            s << A(i, j) << " ";
        s << "\n";
    }
}

template <typename T>
void visual_matlab(matrix<T> ext, std::vector<std::list<size_t> >& ineinc, std::ostream& mfile)
{
  size_t n = ext.ncols();
  size_t k = ext.nrows();
  size_t m = ineinc.size();

  if (n != 4)
  {
    mfile << "% dim != 4. Visualization is not possible\n";
    return;
  }

  mfile << "V = [\n";

  for (size_t i = 0; i < k; i++)
  {
    double d = ext(i, 0);
    if (d == 0) // This is an extreme ray. What should I do?
      for (size_t j = 1; j < n; j++)
        mfile << "  NaN";
    else
      for (size_t j = 1; j < n; j++)
        mfile << "  " << double(ext(i, j)) / d;

    mfile << "\n";
  }
  mfile << "];\n\n";

  size_t maxs = 0;

  for (size_t i = 0; i < m; i++)
    if (ineinc[i].size() > maxs)
       maxs = ineinc[i].size();

  mfile << "F = [\n";

  for (size_t i = 0; i < m; i++)
  {
    size_t s = 0;
    for (std::list<size_t>::iterator cur = ineinc[i].begin();
      cur != ineinc[i].end(); ++cur)
    {
       mfile << " " << (*cur) + 1;
       ++s;
    }
    for (size_t j = s; j < maxs; j++)
       mfile << " NaN";
    mfile << "\n";
  }
  mfile << "];\n\n";

  mfile << "h = polydraw(V, F);\n";

}  

#undef WRITELOG
#undef WRITEOUTPUT
#undef WRITEOUTPUTA
#undef WRITESUMMARY

#endif