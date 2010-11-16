##############################################################################
#
#    Makefile
#
#    Makefile for GNU Make Utility to build the library and
#    other programs from the Arageli Package and run some test and
#    examples.
#
#    This file is a part of the Arageli library.
#
#    Copyright (C) 2010 Anastasya A. Ryzhova
#
#    The Arageli Library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License version 2
#    as published by the Free Software Foundation.
#
#    The Arageli Library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
#
#    We are also open for dual licensing for the whole library or
#    for its particular part. If you are interested to get the library
#    in this way, i.e. not under the GNU General Public License,
#    please contact Arageli Support Service support.arageli@gmail.com.
#
##############################################################################

BIN_DIR = ../../../../bin
ARAGELI_SOURCES_DIR = ../../../../src
SOURCE_DIR = .
INCLUDE_DIR = $(SOURCE_DIR) $(ARAGELI_SOURCES_DIR)../../../../tool/tclap/include/tclap
ARAGELI_LIB_DIR =../../../../lib

PERFCXXFLAGS += -MD -pg -DARAGELI_SOURCE -g -O3 -DARAGELI_DEBUG_LEVEL=0
LINKFLAGS = - L"$(ARAGELI_LIB_DIR)" -lperfarageli

all: $(BIN_DIR)/polyhedron

$(BIN_DIR)/polyhedron: $(ARAGELI_LIB_DIR)libperftarageli.a \
	stdafx.hpp.gch \
	$(patsubst %.cpp,%.o, $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR))))
	make -C $(ARAGELI_SOURCES_DIR)/arageli $(ARAGELI_LIB_DIR)/libperfarageli.a
	$(CXX) -o $(BIN_DIR)/polyhedron $(patsubst %.cpp,%.o, $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR)))) $(LINKFLAGS)

$(INCLUDE_DIR)/stdafx.hpp.gch: $(INCLUDE_DIR)/stdafx.hpp
	$(CC) -x c++-header $(PERFCXXFLAGS) -o $@ $< -I"$(INCLUDE_DIR)"

%.o: %.cpp
	$(CXX) -c $(PERFCXXFLAGS) $< -o $@ -I"$(INCLUDE_DIR)"

$(LIBRARY_DIR)/libperfarageli.a:
	$(MAKE) -C $(ARAGELI_SOURCES_DIR)/arageli

TAGS:
	ctags $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR))) $(wildcard \
	$(addsuffix /*.hpp,$(SOURCE_DIR))) #-f TAGS

clean:
	rm -f $(BIN_DIR)/polyhedron
	rm -f $(wildcard $(addsuffix /*.o,$(SOURCE_DIR))) $(wildcard $(addsuffix /*.d,$(SOURCE_DIR)))

cleanall: clean
	rm -f stdafx.hpp.gch
	rm -f tags

include $(wildcard *.d)
