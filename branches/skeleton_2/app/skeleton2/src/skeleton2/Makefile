##############################################################################
#
#    Makefile for GNU Make Utility to build the skeleton application.
#    Copyright (C) 2010 N. Yu. Zolotykh
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
##############################################################################
ifndef CXXFLAGS
	CXXFLAGS=-O2
endif
INCLUDEDIR=../tools/arageli/src/
LIBDIR=../tools/arageli/lib/

../bin/skeleton: skeleton.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) $^ -L$(LIBDIR) -larageli -o $@ 

clean:
	rm -f ../bin/skeleton
