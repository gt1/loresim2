ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=Intv.hpp

MANPAGES = 

man_MANS = ${MANPAGES}

EXTRA_DIST = ${MANPAGES}

bin_PROGRAMS = longreadgen mapcheckstrict fastareformat

longreadgen_SOURCES = longreadgen.cpp
longreadgen_LDADD = ${LIBMAUS2LIBS}
longreadgen_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
longreadgen_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

mapcheckstrict_SOURCES = mapcheckstrict.cpp
mapcheckstrict_LDADD = ${LIBMAUS2LIBS}
mapcheckstrict_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
mapcheckstrict_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fastareformat_SOURCES = fastareformat.cpp
fastareformat_LDADD = ${LIBMAUS2LIBS}
fastareformat_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fastareformat_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}
