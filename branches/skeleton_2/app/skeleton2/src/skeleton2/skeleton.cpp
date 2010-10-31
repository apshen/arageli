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

/********************************************************************

 Skeleton 02.01.02

 Copyright (C) 2005--2010 N.Yu.Zolotykh
 University of Nizhni Novgorod, Russia
********************************************************************/

/**
    \file
    All extreme rays and basis of a cone
 */


#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "ddm.hpp"
#include "ddmio.hpp"

#include <arageli/arageli.hpp>

using namespace Arageli;

void print_banner()
{
  std::cout << "                                                                         \n";          
  std::cout << "  &&&&&&&    &&                  &&             &&                       \n";          
  std::cout << " &&&   &&&   &&                  &&             &&                       \n";          
  std::cout << " &&&         &&   &&   &&&&&&&   &&   &&&&&&&  &&&&   &&&&&&   &&&&&&&   \n";          
  std::cout << "  &&&&&&&&   &&  &&   &&     &&  &&  &&     &&  &&   &&    &&  &&    &&  \n";          
  std::cout << "        &&&  &&&&&    &&&&&&&&&  &&  &&&&&&&&&  &&   &&    &&  &&    &&  \n";          
  std::cout << " &&     &&&  &&  &&   &&&        &&  &&&        &&   &&    &&  &&    &&  \n";          
  std::cout << "  &&&&&&&&   &&   &&   &&&&&&&   &&   &&&&&&&    &&&  &&&&&&   &&    &&  \n";          
  std::cout << "                                                                         \n";          
  std::cout << "                                                                         \n";          
  std::cout << "                         by Nikolai Yu. Zolotykh                         \n";          
  std::cout << "                                                                         \n";          
  std::cout << "    A lot of help and suggestions from Sergey Lobanov & Sergey Lyalin    \n";          
  std::cout << "                                                                         \n";          
  std::cout << "                  http://www.uic.nnov.ru/~zny/skeleton                   \n";          
  std::cout << "                Skeleton on-line: http://www.arageli.org                 \n";          
  std::cout << "                                                                         \n";          
}

void print_version()
{
  std::cout << "                                                                         \n";
  std::cout << "               ******************************************                \n";
  std::cout << "              ***                                      ***               \n";
  std::cout << "             ****          Skeleton 02.01.02           ****              \n";
  std::cout << "            *****   Compiled:  " << __DATE__ << "  " << __TIME__ << "   *****\n";
  std::cout << "            *****                                      *****             \n";
  std::cout << "            *****         Nikolai Yu. Zolotykh         *****             \n";
  std::cout << "             **** http://www.uic.nnov.ru/~zny/skeleton ****              \n";
  std::cout << "              ***                                      ***               \n";
  std::cout << "               ******************************************                \n";
  std::cout << "                                                                         \n";
}

void print_short_help()
{
  std::cout << "                                                                         \n";          
  std::cout << " Skeleton is an implementation of the Double Description Method          \n";          
  std::cout << " for vertex/facet enumeration of polyhedra                               \n";
  std::cout << " (or extreme rays/facet enumeration for polyhedral cones)                \n";
  std::cout << "                                                                         \n";          
  std::cout << " The program is under the GNU General Public License 2                   \n";          
  std::cout << "                                                                         \n";          
  std::cout << " Type skeleton --version to check the version                            \n";          
  std::cout << "      skeleton --help    for help (try to read the manual too)           \n";          
  std::cout << "      skeleton --copying for more details about license                  \n";          
  std::cout << "                                                                         \n";          
}

void print_copying()
{
  std::cout << "                                                                         \n";          
  std::cout << "  This program is free software; you can redistribute it and/or modify   \n";
  std::cout << "  it under the terms of the GNU General Public License as published by   \n";
  std::cout << "  the Free Software Foundation; either version 2 of the License, or      \n";
  std::cout << "  (at your option) any later version.                                    \n";
  std::cout << "                                                                         \n";
  std::cout << "  This program is distributed in the hope that it will be useful,        \n";
  std::cout << "  but WITHOUT ANY WARRANTY; without even the implied warranty of         \n";
  std::cout << "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          \n";
  std::cout << "  GNU General Public License for more details.                           \n";
  std::cout << "                                                                         \n";
  std::cout << "  You should have received a copy of the GNU General Public License      \n";
  std::cout << "  along with this program; if not, write to the Free Software            \n";
  std::cout << "  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              \n";
  std::cout << "                                                                         \n";
  std::cout << "  To read the license see file copying.                                  \n";          
  std::cout << "                                                                         \n";          
}

void print_help()
{
  std::cout << "Double description method for polyhedral cone\n";
  std::cout << "USAGE: skeleton [inputfilename] [options]\n";
  std::cout << "OPTIONS:\n";
  std::cout << "Order:      {--minindex} / --maxindex / --lexmin / --lexmax / --random /\n";
  std::cout << "            --mincutoff / --maxcutoff / --minpairs / --maxpairs\n";
  std::cout << "Prefixed:   {--prefixedorder} / --noprefixedorder  \n";
  std::cout << "Graphadj:   {--graphadj} / --nographadj  \n";
  std::cout << "Plusplus:   {--plusplus} / --noplusplus  \n";
  std::cout << "Arithmetic: {--bigint} / --int / --float / --rational\n";
  std::cout << "            --zerotol value (by default 1e-8)\n";
  std::cout << "Find edges/adjacency: \n";
  std::cout << "            --edges / {--noedges}\n";
  std::cout << "            --adjacency / {--noadjacency}\n";
  std::cout << "            --ridges / {--noridges}\n";
  std::cout << "            --facetadjacency / {--nofacetadjacency}\n";
  std::cout << "            --verifyine / {--noverifyine}\n";
  std::cout << "File names: --inputfile inputfilename \n";
  std::cout << "            -o outputfilename | --outputfile outputfilename  \n";
  std::cout << "            --logfile logfilename  \n";
  std::cout << "            --summaryfile summaryfilename  \n";
  std::cout << "Inp.option: --inputfromstdin / {--noinputfromstdin}  \n";
  std::cout << "            {--simpleformat} / --skeletonformat / --avisfukudaformat \n";
  std::cout << "Output options:\n";
  std::cout << "            {--outputinfile} / --nooutputinfile \n";
  std::cout << "            --outputonstdout / {--nooutputonstdout}  \n";
  std::cout << "            {--loginfile} / --nologinfile \n";
  std::cout << "            {--logonstdout} / --nologonstdout  \n";
  std::cout << "            {--summaryinfile} / --nosummaryinfile \n";
  std::cout << "            {--summaryonstdout} / --nosummaryonstdout  \n";
  std::cout << "            --silence  \n";
  std::cout << "Out.options:--ine / {--noine}  \n";
  std::cout << "            --equ / {--noequ}  \n";
  std::cout << "            {--ext} / --noext  \n";
  std::cout << "            {--bas} / --nobas  \n";
  std::cout << "            --dis / {--nodis}  \n";
  std::cout << "            --extinc / {--noextinc}  \n";
  std::cout << "            --ineinc / {--noineinc}  \n";
  std::cout << "            --matrices / --nomatrices  \n";
  std::cout << "Out.options:{--log} / --nolog  \n";
  std::cout << "            {--summary} / --nosummary  \n";
//  std::cout << "Matlab format:  \n";
// std::cout << "            --visualmatlab / {--novisualmatlab}  \n";
  std::cout << "Other:      -h | --help     print this help message and exit\n";
  std::cout << "            -v | --version  print program version and exit\n";
  std::cout << "            --copying       license info\n";
}

#define FILENAMELEN 200
#define FILENAMEEXTLEN 4
#define DEFAULTOUTPUTFILENAME "skeleton.out"
#define DEFAULTLOGFILENAME "skeleton.log"
#define DEFAULTSUMMARYFILENAME "skeleton.sum"
#define DEFAULTOUTPUTFILENAMEEXT ".out"
#define DEFAULTLOGFILENAMEEXT ".log"
#define DEFAULTSUMMARYFILENAMEEXT ".sum"
#define DEFAULTMFILENAMEEXT ".m"

template <class T>
void read_matrices_call_ddm_and_print_results(
  matrix<T>& ine, matrix<T>& equ, 
  matrix<T>& ext, matrix<T>& bas, matrix<T>& dis, matrix<size_t>& edges_ind, 
  INT64& totalnumrays,  INT64& totalnumedges,
  int edgesflag, int adjacencyflag,
  int ridgesflag, int facetadjacencyflag,
  int extincflag, int ineincflag, int inetypeflag,
  int modification, int prefixed_order, int graphadj, int plusplus, int arith, 
  int intarith, const T& eps, 
  int argc, char *argv[], 
  char* inputfilename, char* outputfilename, char* logfilename, char* summaryfilename, 
  int output_ine, int output_equ, int output_ext, int output_bas, int output_dis,
  int inputoutputformat, int inputfromstdin, std::istream& inputfile, 
  int outputonstdout, int outputfileflag, std::ostream& outputfile,
  int logonstdout, int logfileflag, std::ostream& logfile,
  int summaryonstdout, int summaryfileflag, std::ostream& summaryfile)
{
  try
  {
      input_matrices(ine, equ,
          inputoutputformat, inputfromstdin, inputfile, 
          logonstdout, logfileflag, logfile);
  }
  catch(...)                                     
  {                                              
    std::cout << "Input error. Skeleton terminated.";
    return;                                    
  } 

  time_t begin_time;                                                                               
  time(&begin_time);                                                                               
  timer time_sec;                                                                                  
                                                                                                   
  ddm(ine, equ, ext, bas, dis, edges_ind, totalnumrays, totalnumedges, edgesflag | adjacencyflag,                                            
    modification, prefixed_order, graphadj, plusplus, intarith, eps,                               
    logonstdout, logfileflag, logfile);                                                            
                                                                                                  
  time_sec.stop();                                                                                 
  time_t end_time;                                                                                 
  time(&end_time);

  output_results(ine, equ, ext, bas, dis, edges_ind,
    output_ine, output_equ, output_ext, output_bas, output_dis, edgesflag,
    outputonstdout, outputfileflag, outputfile);

  output_summary(ine, equ, ext, bas, dis,
    totalnumrays, totalnumedges,
    modification, prefixed_order, graphadj, plusplus, arith, eps,
    argc, argv, 
    inputfilename, outputfilename, logfilename, summaryfilename, 
    summaryonstdout, summaryfileflag, summaryfile);                                                                                                

  output_time(begin_time, end_time, time_sec,
    summaryonstdout, summaryfileflag, summaryfile);

  if (adjacencyflag)
  {
    std::vector<std::list<size_t> > adj;
    post_construct_adjacency(ext.nrows(), edges_ind, adj);
    print_adjacency("* Adjacency:\n", adj, outputonstdout, outputfileflag, outputfile);
  }
  
  if (extincflag || ineincflag || inetypeflag || ridgesflag || facetadjacencyflag)
  {
    matrix<int> dis_int;
    disarray2intarray(dis, dis_int, eps);
    if (ineincflag)
    {
      std::vector<std::list<size_t> > ineinc;
      post_construct_ineinc(dis_int, ineinc);
      print_adjacency("* Inequalities-to-rays incidence:\n", 
        ineinc, outputonstdout, outputfileflag, outputfile);
    }
    if (extincflag)
    {
      std::vector<std::list<size_t> > extinc;
      post_construct_extinc(dis_int, extinc);
      print_adjacency("* Rays-to-inequalities incidence:\n", 
         extinc, outputonstdout, outputfileflag, outputfile);
    }
    if (inetypeflag || ridgesflag || facetadjacencyflag)
    {
      std::vector<int> ine_type;
      matrix<size_t> ridges;
      post_construct_inetype_and_ridges(dis_int, ine_type, ridges);
      if (inetypeflag)
      {
        print_no("* Implicit equations:\n", ine_type, INE_IMPLICIT_EQUATION, 
          outputonstdout, outputfileflag, outputfile); 
        print_no("* Redundant inequalities:\n", ine_type, INE_REDUNDANT, 
          outputonstdout, outputfileflag, outputfile);
      }
      if (ridgesflag)
      {
        print_edges("* Ridges:\n", ridges, outputonstdout, outputfileflag, outputfile);
      }
      if (facetadjacencyflag)
      {
        std::vector<std::list<size_t> > facetadj;
        post_construct_adjacency(ine.nrows(), ridges, facetadj);
        print_adjacency("* Facet adjacency:\n", facetadj, 
          outputonstdout, outputfileflag, outputfile);
      }
    }
  }
}                                                                                          
                                                                                           
int main(int argc, char *argv[])
{
  int modification = DDM_MIN_INDEX;
  int output_ine = 0;
  int output_equ = 0;
  int output_bas = 1;
  int output_ext = 1;
  int output_dis = 0;
  int extincflag = 0;
  int ineincflag = 0;
  int edgesflag = 0;
  int adjacencyflag = 0;
  int ridgesflag = 0;
  int facetadjacencyflag = 0;
  int inetypeflag = 0;
  int prefixed_order = -1;
  int graphadj = 1;
  int inputfromstdin = 0;
  int inputoutputformat = SIMPLEFORMAT;
  int logfileflag = 1;
  int logonstdout = 1;
  int summaryfileflag = 1;
  int summaryonstdout = 1;
  int outputfileflag = 1;
  int outputonstdout = 0;
  int arith = ARITHBIGINT;
  int plusplus = 1;
  int visualmatlab = 0;
  int there_were_errors = 0;
  double zerotol = 1e-8;

  char inputfilename[FILENAMELEN];
  char outputfilename[FILENAMELEN+FILENAMEEXTLEN];
  char logfilename[FILENAMELEN+FILENAMEEXTLEN];
  char summaryfilename[FILENAMELEN+FILENAMEEXTLEN];
  char mfilename[FILENAMELEN+FILENAMEEXTLEN];

  inputfilename[0] = '\0';
  outputfilename[0] = '\0';
  logfilename[0] = '\0';
  summaryfilename[0] = '\0';
  mfilename[0] = '\0';

  if (argc == 1)
  {
    print_banner();
    print_short_help();
    return 0;
  }

  // arguments parsing:
  for(int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
    {
      print_help();
      return 0;
    }
    else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0)
    {
      print_version();
      return 0;
    }
    if (strcmp(argv[i], "--copying") == 0)
    {
      print_copying();
      return 0;
    }
    else if (strcmp(argv[i], "--ine") == 0)
      output_ine = 1; 
    else if (strcmp(argv[i], "--equ") == 0)
      output_equ = 1; 
    else if (strcmp(argv[i], "--ext") == 0)
      output_ext = 1; 
    else if (strcmp(argv[i], "--dis") == 0)
      output_dis = 1; 
    else if (strcmp(argv[i], "--bas") == 0)
      output_bas = 1; 
    else if (strcmp(argv[i], "--log") == 0)
    {
      logfileflag = 1; 
      logonstdout = 1; 
    }
    else if (strcmp(argv[i], "--summary") == 0)
    {
      summaryfileflag = 1; 
      summaryonstdout = 1; 
    }
    else if (strcmp(argv[i], "--extinc") == 0)
      extincflag = 1; 
    else if (strcmp(argv[i], "--ineinc") == 0)
      ineincflag = 1; 
    else if (strcmp(argv[i], "--noine") == 0)
      output_ine = 0; 
    else if (strcmp(argv[i], "--noequ") == 0)
      output_equ = 0; 
    else if (strcmp(argv[i], "--noext") == 0)
      output_ext = 0; 
    else if (strcmp(argv[i], "--nodis") == 0)
      output_dis = 0; 
    else if (strcmp(argv[i], "--nobas") == 0)
      output_bas = 0; 
    else if (strcmp(argv[i], "--nolog") == 0)
    {
      logfileflag = 0; 
      logonstdout = 0; 
    }
    else if (strcmp(argv[i], "--nosummary") == 0)
    {
      summaryfileflag = 0; 
      summaryonstdout = 0; 
    }
    else if (strcmp(argv[i], "--noextinc") == 0)
      extincflag = 0; 
    else if (strcmp(argv[i], "--noineinc") == 0)
      ineincflag = 0; 
    else if (strcmp(argv[i], "--matrices") == 0)
    {
      output_ine = 1;
      output_equ = 1;
      output_ext = 1;
      output_dis = 1;
      output_bas = 1;
    }
    else if (strcmp(argv[i], "--nomatrices") == 0)
    {
      output_ine = 0;
      output_equ = 0;
      output_ext = 0;
      output_dis = 0;
      output_bas = 0;
    }
    else if(strcmp(argv[i], "--adjacency") == 0)
      adjacencyflag = 1;
    else if(strcmp(argv[i], "--noadjacency") == 0)
      adjacencyflag = 0;
    else if(strcmp(argv[i], "--edges") == 0)
      edgesflag = 1;
    else if(strcmp(argv[i], "--noedges") == 0)
      edgesflag = 0;
    else if(strcmp(argv[i], "--ridges") == 0)
      ridgesflag = 1;
    else if(strcmp(argv[i], "--noridges") == 0)
      ridgesflag = 0;
    else if(strcmp(argv[i], "--facetadjacency") == 0)
      facetadjacencyflag = 1;
    else if(strcmp(argv[i], "--nofacetadjacency") == 0)
      facetadjacencyflag = 0;
    else if(strcmp(argv[i], "--verifyine") == 0)
      inetypeflag = 1;
    else if(strcmp(argv[i], "--noverifyine") == 0)
      inetypeflag = 0;
    else if(strcmp(argv[i], "--lexmin") == 0)
      modification = DDM_LEX_MIN;
    else if(strcmp(argv[i], "--lexmax") == 0)
      modification = DDM_LEX_MAX;
    else if(strcmp(argv[i], "--minindex") == 0)
      modification = DDM_MIN_INDEX;
    else if(strcmp(argv[i], "--maxindex") == 0)
      modification = DDM_MAX_INDEX;
    else if(strcmp(argv[i], "--mincutoff") == 0)
      modification = DDM_MIN_CUT_OFF;
    else if(strcmp(argv[i], "--maxcutoff") == 0)
      modification = DDM_MAX_CUT_OFF;
    else if(strcmp(argv[i], "--minpairs") == 0)
      modification = DDM_MIN_EDGES;
    else if(strcmp(argv[i], "--maxpairs") == 0)
      modification = DDM_MAX_EDGES;
    else if(strcmp(argv[i], "--random") == 0)
      modification = DDM_RANDOM;
    else if(strcmp(argv[i], "--prefixedorder") == 0)
      prefixed_order = 1;
    else if(strcmp(argv[i], "--noprefixedorder") == 0)
      prefixed_order = 0;
    else if(strcmp(argv[i], "--graphadj") == 0)
      graphadj = 1;
    else if(strcmp(argv[i], "--nographadj") == 0)
      graphadj = 0;
    else if(strcmp(argv[i], "--plusplus") == 0)
      plusplus = 1;
    else if(strcmp(argv[i], "--noplusplus") == 0)
      plusplus = 0;
    else if(strcmp(argv[i], "--inputfromstdin") == 0)
      inputfromstdin = 1;
    else if(strcmp(argv[i], "--noinputfromstdin") == 0)
      inputfromstdin = 0;
    else if(strcmp(argv[i], "--simpleformat") == 0)
      inputoutputformat = SIMPLEFORMAT;
    else if(strcmp(argv[i], "--skeletonformat") == 0)
      inputoutputformat = SKELETONFORMAT;
    else if(strcmp(argv[i], "--avisfukudaformat") == 0)
      inputoutputformat = AVISFUKUDAFORMAT;
    else if(strcmp(argv[i], "--inputfile") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "There is no name for input file\n";
        there_were_errors = 1;
      }
      else if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for input file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(inputfilename, argv[i]);
      }
    }
    else if(strcmp(argv[i], "--outputfile") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "Skeleton Error: There is no name for output file\n";
        there_were_errors = 1;
      }
      else if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for output file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(outputfilename, argv[i]);
      }
    }
    else if(strcmp(argv[i], "--logfile") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "Skeleton Error: There is no name for log file\n";
        there_were_errors = 1;
      }
      else if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for log file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(logfilename, argv[i]);
      }
    }
    else if(strcmp(argv[i], "--summaryfile") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "Skeleton Error: There is no name for summary file\n";
        there_were_errors = 1;
      }
      else if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for summary file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(summaryfilename, argv[i]);
      }
    }
    else if(strcmp(argv[i], "--mfile") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "Skeleton Error: There is no name for m-file\n";
        there_were_errors = 1;
      }
      else if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for m-file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(mfilename, argv[i]);
      }
    }
    else if(strcmp(argv[i], "--outputinfile") == 0)
      outputfileflag = 1;
    else if(strcmp(argv[i], "--nooutputinfile") == 0)
      outputfileflag = 0;
    else if(strcmp(argv[i], "--outputonstdout") == 0)
      outputonstdout = 1;
    else if(strcmp(argv[i], "--nooutputonstdout") == 0)
      outputonstdout = 0;
    else if(strcmp(argv[i], "--loginfile") == 0)
      logfileflag = 1;
    else if(strcmp(argv[i], "--nologinfile") == 0)
      logfileflag = 0;
    else if(strcmp(argv[i], "--logonstdout") == 0)
      logonstdout = 1;
    else if(strcmp(argv[i], "--nologonstdout") == 0)
      logonstdout = 0;
    else if(strcmp(argv[i], "--summaryinfile") == 0)
      summaryfileflag = 1;
    else if(strcmp(argv[i], "--nosummaryinfile") == 0)
      summaryfileflag = 0;
    else if(strcmp(argv[i], "--summaryonstdout") == 0)
      summaryonstdout = 1;
    else if(strcmp(argv[i], "--nosummaryonstdout") == 0)
      summaryonstdout = 0;
    else if(strcmp(argv[i], "--silence") == 0)
      outputonstdout = logonstdout = summaryonstdout = 0;
    else if(strcmp(argv[i], "--visualmatlab") == 0)
      visualmatlab = 1;
    else if(strcmp(argv[i], "--novisualmatlab") == 0)
      visualmatlab = 0;
    else if(strcmp(argv[i], "--bigint") == 0)
      arith = ARITHBIGINT;
    else if(strcmp(argv[i], "--int") == 0)
      arith = ARITHINT;
    else if(strcmp(argv[i], "--float") == 0)
      arith = ARITHFLOAT;
    else if(strcmp(argv[i], "--rational") == 0)
      arith = ARITHRATIONAL;
    else if(strcmp(argv[i], "--zerotol") == 0)
    {
      i++;
      if (i >= argc)
      {
        std::cout << "There is no zero tolerance\n";
        there_were_errors = 1;
      }
      else
      {
        zerotol = atof(argv[i]);
      }
    }
    else 
      if (strlen(argv[i]) >= FILENAMELEN)
      {
        std::cout << "Skeleton Error: name for input file is too long\n";
        there_were_errors = 1;
      }
      else
      {
        strcpy(inputfilename, argv[i]);
      }
    /*{
      std::cout << "Skeleton Error: Skeleton Error: unrecognized option: "  
                << argv[i] << "\n";
      there_were_errors = 1;
    }*/
  }

  if (there_were_errors)
  {
    std::cout << "There were errors. Skeleton terminated\n";
    return 1;
  }

  if (inputfilename[0] == '\0' && !inputfromstdin)
  {
    std::cout << "No input file name. Skeleton terminated\n";
    return 1;
    /* in versions before 02.00.04 we output from cin in this case:
    inputfromstdin = 1;
    if (outputfileflag && outputfilename[0] == '\0')
      strcpy(outputfilename, DEFAULTOUTPUTFILENAME);
    if (logfileflag && logfilename[0] == '\0')
      strcpy(logfilename, DEFAULTLOGFILENAME);
    if (summaryfileflag && summaryfilename[0] == '\0')
      strcpy(summaryfilename, DEFAULTSUMMARYFILENAME);
    */
  }
  if (outputfileflag && outputfilename[0] == '\0')
  {
    if (inputfilename[0] == '\0')
    {
      strcpy(outputfilename, DEFAULTOUTPUTFILENAME);
    }
    else
    {
      strcpy(outputfilename, inputfilename);
      strcat(outputfilename, DEFAULTOUTPUTFILENAMEEXT);
    }
  }
  if (logfileflag && logfilename[0] == '\0')
  {
    if (inputfilename[0] == '\0')
    {
      strcpy(logfilename, DEFAULTLOGFILENAME);
    }
    else
    {
      strcpy(logfilename, inputfilename);
      strcat(logfilename, DEFAULTLOGFILENAMEEXT);
    }
  }
  if (summaryfileflag && summaryfilename[0] == '\0')
  {
    if (inputfilename[0] == '\0')
    {
      strcpy(summaryfilename, DEFAULTSUMMARYFILENAME);
    }
    else
    {
      strcpy(summaryfilename, inputfilename);
      strcat(summaryfilename, DEFAULTSUMMARYFILENAMEEXT);
    }
  }
  if (visualmatlab && mfilename[0] == '\0')
  {
    strcpy(mfilename, inputfilename);
    char *ch;
    while (ch = strchr(mfilename, '.'))
      *ch = '_';
    strcat(mfilename, DEFAULTMFILENAMEEXT);
  }

  std::ifstream inputfile;
  std::ofstream outputfile;
  std::ofstream logfile;
  std::ofstream summaryfile;
  std::ofstream mfile;

  if (!inputfromstdin)
  {
    inputfile.open(inputfilename, std::ios::in);
    if (!inputfile)
    {
      std::cout << "Skeleton Error: can not open file " << inputfilename << " for input\n";
      std::cout << "                or incorrect option\n";
      return 1;
    }
  }

  if (outputfileflag)
  {
     outputfile.open(outputfilename, std::ios::out);
     if (!outputfile)
     {
       std::cout << "Skeleton Error: can not open file " << outputfilename << " for output\n";
       return 1;
     }
  }

  if (logfileflag)
  {
    logfile.open(logfilename, std::ios::out);
    if (!logfile)
    {
      std::cout << "Skeleton Error: can not open log file " << logfilename << " for output\n";
      return 1;
    }
  }

  if (summaryfileflag)
  {
    summaryfile.open(summaryfilename, std::ios::out);
    if (!summaryfile)
    {
      std::cout << "Skeleton Error: can not open summary file " << summaryfilename << " for output\n";
      return 1;
    }
  }

  if (visualmatlab)
  {
    mfile.open(mfilename, std::ios::out);
    if (!mfile)
    {
      std::cout << "Skeleton Error: can not open m-file " << summaryfilename << " for output\n";
      return 1;
    }
  }

  if (modification == DDM_MIN_CUT_OFF || modification == DDM_MAX_CUT_OFF 
    ||  modification == DDM_MIN_EDGES || modification == DDM_MAX_EDGES)
  {
    if (prefixed_order == 1)
      std::cout << "Skeleton Warning: for this modification prefixed order is not possible\n";
    prefixed_order = 0;
  }
  else if (prefixed_order == -1)
  {
    prefixed_order = 1; // by default for others modifications
  }

  matrix<size_t> edges_ind;
  //std::vector<std::list<size_t> > adj;
  //std::vector<std::list<size_t> > extinc;
  //std::vector<std::list<size_t> > ineinc;
  INT64 totalnumrays, totalnumedges;

  if (arith == ARITHBIGINT)
  {
    matrix<big_int> ine, equ, ext, dis, bas;
    int intarith = 1;
    big_int eps = 0;
    read_matrices_call_ddm_and_print_results(
        ine, equ, 
        ext, bas, dis, edges_ind, 
        totalnumrays, totalnumedges,
        edgesflag,  adjacencyflag,
        ridgesflag, facetadjacencyflag,
        extincflag, ineincflag, inetypeflag,
        modification, prefixed_order, graphadj, plusplus, arith,
        intarith, eps, 
        argc, argv,
        inputfilename, outputfilename, logfilename, summaryfilename, 
        output_ine, output_equ, output_ext, output_bas, output_dis,
        inputoutputformat, inputfromstdin, inputfile, 
        outputonstdout, outputfileflag, outputfile,
        logonstdout, logfileflag, logfile,
        summaryonstdout, summaryfileflag, summaryfile);
  }
  else if (arith == ARITHINT)
  {
    matrix<long> ine, equ, ext, dis, bas;
    int intarith = 1;
    long eps = 0;
    read_matrices_call_ddm_and_print_results(
        ine, equ, 
        ext, bas, dis, edges_ind, 
        totalnumrays, totalnumedges,
        edgesflag,  adjacencyflag,
        ridgesflag, facetadjacencyflag,
        extincflag, ineincflag, inetypeflag,
        modification, prefixed_order, graphadj, plusplus, arith,
        intarith, eps, 
        argc, argv,
        inputfilename, outputfilename, logfilename, summaryfilename, 
        output_ine, output_equ, output_ext, output_bas, output_dis,
        inputoutputformat, inputfromstdin, inputfile, 
        outputonstdout, outputfileflag, outputfile,
        logonstdout, logfileflag, logfile,
        summaryonstdout, summaryfileflag, summaryfile);
  }
  else if (arith == ARITHFLOAT)
  {
    matrix<double> ine, equ, ext, dis, bas;
    int intarith = 0;
    double eps = zerotol;
    read_matrices_call_ddm_and_print_results(
        ine, equ, 
        ext, bas, dis, edges_ind, 
        totalnumrays, totalnumedges,
        edgesflag,  adjacencyflag,
        ridgesflag, facetadjacencyflag,
        extincflag, ineincflag, inetypeflag,
        modification, prefixed_order, graphadj, plusplus, arith,
        intarith, eps, 
        argc, argv,
        inputfilename, outputfilename, logfilename, summaryfilename, 
        output_ine, output_equ, output_ext, output_bas, output_dis,
        inputoutputformat, inputfromstdin, inputfile, 
        outputonstdout, outputfileflag, outputfile,
        logonstdout, logfileflag, logfile,
        summaryonstdout, summaryfileflag, summaryfile);
  }
  else if (arith == ARITHRATIONAL)
  {
    matrix<rational<big_int> > ine, equ, ext, dis, bas;
    int intarith = 0;
    rational<big_int> eps = 0;
    read_matrices_call_ddm_and_print_results(
        ine, equ, 
        ext, bas, dis, edges_ind, 
        totalnumrays, totalnumedges,
        edgesflag,  adjacencyflag,
        ridgesflag, facetadjacencyflag,
        extincflag, ineincflag, inetypeflag,
        modification, prefixed_order, graphadj, plusplus, arith,
        intarith, eps, 
        argc, argv,
        inputfilename, outputfilename, logfilename, summaryfilename, 
        output_ine, output_equ, output_ext, output_bas, output_dis,
        inputoutputformat, inputfromstdin, inputfile, 
        outputonstdout, outputfileflag, outputfile,
        logonstdout, logfileflag, logfile,
        summaryonstdout, summaryfileflag, summaryfile);
  }

  inputfile.close();
  outputfile.close();
  logfile.close();
  summaryfile.close();
  mfile.close();

  return 0;
}

