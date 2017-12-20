#
# automake template for src/main directory
#
# $Id: Makefile.am 134 2007-06-07 20:05:52Z dkumar $
#

bin_PROGRAMS = titan

titan_SOURCES = forward_solver.C compare_key.C delete_tab.C datread.C \
		hilbert.C init_piles.C update_topo.C restart.C hpfem.C

titan_CXXFLAGS = -I$(srcdir)/../gisapi

titan_LDADD = ../adjoint/libadjoint.a ../gisapi/libgisapi.a ../adapt/libadapt.a @FLIBS@ \
  ../geoflow/libgeoflow.a ../datstr/libdatstr.a  ../useful/libuseful.a \
  ../repartition/librepartition.a ../tecplot/libtecplot.a 


