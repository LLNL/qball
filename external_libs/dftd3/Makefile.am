## Copyright (C) 2017 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: Makefile.am 7479 2011-02-18 11:33:08Z xavier $


noinst_LIBRARIES = libdftd3.a

dftd3dir = $(pkgdatadir)/dftd3

libdftd3_a_SOURCES =            \
	api.f90                 \
	common.f90              \
	core.f90                \
	pars.F90                \
	sizes.f90               \
	interface.f90

dist_dftd3_DATA =               \
	pars.dat

AM_FCFLAGS = -DSHARE_DIR='"$(pkgdatadir)"'

interface.o: interface.f90 api.o
api.o: api.f90 sizes.o core.o common.o
core.o: core.f90 pars.o sizes.o common.o
pars.o: pars.F90 common.o sizes.o

check:
	for f in $(dist_dftd3_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

clean-local:
	if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_dftd3_DATA); \
	fi

.PHONY: check

CLEANFILES = *~ *.bak *.mod *.MOD *.il *.d *.pc* ifc*
