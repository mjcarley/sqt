## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libsqt.la

bin_PROGRAMS = \
	sqt-tests

sqt_tests_LDADD = \
	$(top_builddir)/src/.libs/libsqt.a 
