
WERR = #-Werror

AM_CFLAGS = -O2 -std=c99 -fno-common -Wextra -Wall \
    -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-align \
    -Wwrite-strings -Wnested-externs -Wundef -Wsign-compare \
    -Wmissing-declarations -Wmissing-format-attribute -Wold-style-definition \
    -Wdisabled-optimization -Wformat-y2k \
    -Winit-self -Wmissing-include-dirs -Wmissing-noreturn -Wredundant-decls \
    -Wvariadic-macros -Wvolatile-register-var -Wstrict-prototypes $(WERR)

AM_CPPFLAGS = $(CPL_INCLUDES) -DMOLECFIT_BASEDIR='"$(prefix)"' -DMOLECFIT_DATADIR='"$(profilesdir)"'
AM_LDFLAGS = $(CPL_LDFLAGS)

# link against static variant as most symbols are private
LDADD = $(top_builddir)/src/.libs/libmolecfit.a $(LIBCPLCORE) $(LIBCPLUI) $(LIBCPLDFS) $(LIBCPLDRS) $(MATHLIB)

SUBDIRS = . test

noinst_includedir =
noinst_include_HEADERS = mf_basic.h mf_par.h mf_molecfit.h mf_readspec.h mf_atm.h \
                         mf_mpfit.h mpfit.h mf_modsim.h mf_lib.h mf_lnfl.h \
                         mf_lblrtm.h mf_trans.h mf_conv.h mf_corr.h mf_plot.h

include_HEADERS = molecfit.h

EXTRA_DIST = DISCLAIMER

bin_PROGRAMS = molecfit calctrans corrfilelist preptable prepguitable \
               extract_gdas_profiles calctrans_lblrtm calctrans_convolution

lib_LTLIBRARIES = libmolecfit.la

libmolecfit_la_SOURCES = \
    mf_par.c \
    mf_basic.c \
    mf_molecfit.c \
    mf_readspec.c \
    mf_corr.c \
    mf_conv.c \
    mf_atm.c \
    mf_mpfit.c \
    mpfit.c \
    mf_modsim.c \
    mf_lib.c \
    mf_trans.c \
    mf_plot.c \
    mf_lnfl.c \
    mf_lblrtm.c

# explicit dependencies needed as we are linking the static library
molecfit_DEPENDENCIES = libmolecfit.la
molecfit_SOURCES = molecfit.c

calctrans_DEPENDENCIES = libmolecfit.la
calctrans_SOURCES = calctrans.c

calctrans_lblrtm_DEPENDENCIES = libmolecfit.la
calctrans_lblrtm_SOURCES = calctrans_lblrtm.c

calctrans_convolution_DEPENDENCIES = libmolecfit.la
calctrans_convolution_SOURCES = calctrans_convolution.c

corrfilelist_DEPENDENCIES = libmolecfit.la
corrfilelist_SOURCES = corrfilelist.c

preptable_DEPENDENCIES = libmolecfit.la
preptable_SOURCES = preptable.c

prepguitable_DEPENDENCIES = libmolecfit.la
prepguitable_SOURCES = prepguitable.c

extract_gdas_profiles_DEPENDENCIES = libmolecfit.la
extract_gdas_profiles_SOURCES = extract_gdas_profiles.c
