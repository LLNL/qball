## Process this file with automake to produce Makefile.in

## Copyright (C) 2015 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

pdpbesolstrdir = $(pkgdatadir)/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbesol_stringent/

dist_pdpbesolstr_DATA =                            \
	Ag.upf                                  \
	Al.upf                                  \
	Ar.upf                                  \
	As.upf                                  \
	Au.upf                                  \
	Ba.upf                                  \
	Be.upf                                  \
	Bi.upf                                  \
	Br.upf                                  \
	B.upf                                   \
	Ca.upf                                  \
	Cd.upf                                  \
	Cl.upf                                  \
	Co.upf                                  \
	Cr.upf                                  \
	Cs.upf                                  \
	C.upf                                   \
	Cu.upf                                  \
	Fe.upf                                  \
	F.upf                                   \
	Ga.upf                                  \
	Ge.upf                                  \
	He.upf                                  \
	Hf.upf                                  \
	Hg.upf                                  \
	H.upf                                   \
	In.upf                                  \
	Ir.upf                                  \
	I.upf                                   \
	Kr.upf                                  \
	K.upf                                   \
	La.upf                                  \
	Li.upf                                  \
	Lu.upf                                  \
	Mg.upf                                  \
	Mn.upf                                  \
	Mo.upf                                  \
	Na.upf                                  \
	Nb.upf                                  \
	Ne.upf                                  \
	Ni.upf                                  \
	N.upf                                   \
	Os.upf                                  \
	O.upf                                   \
	Pb.upf                                  \
	Pd.upf                                  \
	Po.upf                                  \
	Pt.upf                                  \
	P.upf                                   \
	Rb.upf                                  \
	Re.upf                                  \
	Rh.upf                                  \
	Rn.upf                                  \
	Ru.upf                                  \
	Sb.upf                                  \
	Sc.upf                                  \
	Se.upf                                  \
	Si.upf                                  \
	Sn.upf                                  \
	Sr.upf                                  \
	S.upf                                   \
	Ta.upf                                  \
	Tc.upf                                  \
	Te.upf                                  \
	Ti.upf                                  \
	Tl.upf                                  \
	V.upf                                   \
	W.upf                                   \
	Xe.upf                                  \
	Y.upf                                   \
	Zn.upf                                  \
	Zr.upf

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_pdpbesolstr_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_pdpbesolstr_DATA); \
	fi

.PHONY: check

