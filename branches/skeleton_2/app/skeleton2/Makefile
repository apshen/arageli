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

.PHONY: skeleton runexamples check clean

skeleton: arageli
	$(MAKE) -C src

arageli:
	$(MAKE) -C tools/arageli arageli

runexamples: skeleton
	cd examples && $(SHELL) runexamples.sh

check: runexamples
	for name in $(notdir $(wildcard $(addsuffix /*.out, examples/golden.out))); do \
		if diff -q --strip-trailing-cr $(addsuffix /$$name, examples) $(addsuffix /$$name, examples/golden.out); \
			then echo $$name - passed; \
			else echo $$name - not passed; \
		fi; \
	done

clean:
	$(MAKE) -C tools/arageli clean
	$(MAKE) -C src clean
