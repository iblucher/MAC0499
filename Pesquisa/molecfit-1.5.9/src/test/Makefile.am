
WERR = #-Werror

AM_CFLAGS = -O2 -std=c99 -fno-common -Wextra -Wall \
    -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-align \
    -Wwrite-strings -Wnested-externs -Wundef -Wsign-compare \
    -Wmissing-declarations -Wmissing-format-attribute -Wold-style-definition \
    -Wdisabled-optimization -Wformat-y2k \
    -Winit-self -Wmissing-include-dirs -Wmissing-noreturn -Wredundant-decls \
    -Wvariadic-macros -Wvolatile-register-var -Wstrict-prototypes $(WERR)

AM_CPPFLAGS = $(CPL_INCLUDES) -I$(top_srcdir)/src -DMOLECFIT_BASEDIR='"$(prefix)"' -DMOLECFIT_DATADIR='"$(profilesdir)"' -DMOLECFIT_SOURCEDIR='"$(srcdir)"'
AM_LDFLAGS = $(CPL_LDFLAGS)

LDADD = $(top_builddir)/src/.libs/libmolecfit.a $(LIBCPLCORE) $(LIBCPLUI) $(LIBCPLDFS) $(LIBCPLDRS) $(MATHLIB)

DISTCLEANFILES = \
                 output/* \
                 vgcore.* \
                 *.gcno \
                 *.gcda \
                 *.out \
                 *.log \
                 *.trs \
                 .logfile \
                 Makefile.in
                 
TESTS = test_mf_basic \
        test_mf_par \
        test_mf_readspec \
        test_mf_atm_readatm \
        test_mf_atm_readgdas \
        test_mf_atm_adaptenv \
        test_mf_atm_writeatm \
        test_mf_plot

check_PROGRAMS = test_mf_basic \
                 test_mf_par \
                 test_mf_readspec \
                 test_mf_atm_readatm \
                 test_mf_atm_readgdas \
                 test_mf_atm_adaptenv \
                 test_mf_atm_writeatm \
                 test_mf_plot


#installcheckdir=$(examplesdir)
#installcheck_PROGRAMS = \
#                        test_mf_atm_createatm \
#                        test_mf_molecfit \
#                        test_mf_trans \
#                        test_mf_corr \
#                        calc_pwv
                                  
EXTRA_DIST = input/2009-05-01T06.gdas \
             input/crires_spec_jitter_extracted_0000.fits \
             input/GDAS_t0_s1.atm \
             input/molecfit_test_crires_fit.atm \
             input/molecfit_test_crires_fit.res \
             input/molecfit_test_crires.fits \
             input/molecfit_test_crires_tac.fits \
             input/testdat_plot_xy.dat \
             input/testspec1.dat \
             input/testspec2.dat \
             config/exclude_crires.dat \
             config/molecfit_test_crires.par \
             config/test_crires.list

test_mf_basic_SOURCES = test_mf_basic.c

test_mf_par_SOURCES = test_mf_par.c

test_mf_readspec_SOURCES = test_mf_readspec.c

test_mf_atm_readatm_SOURCES = test_mf_atm_readatm.c

test_mf_atm_readgdas_SOURCES = test_mf_atm_readgdas.c

test_mf_atm_adaptenv_SOURCES = test_mf_atm_adaptenv.c

test_mf_atm_writeatm_SOURCES = test_mf_atm_writeatm.c

test_mf_plot_SOURCES = test_mf_plot.c

#test_mf_atm_createatm_SOURCES = test_mf_atm_createatm.c

#test_mf_molecfit_SOURCES = test_mf_molecfit.c

#test_mf_trans_SOURCES = test_mf_trans.c

#test_mf_corr_SOURCES = test_mf_corr.c

#calc_pwv_SOURCES = calc_pwv.c

#installcheck-local:
#	$(DESTDIR)${examplesdir}/test_mf_atm_createatm;\
#	$(DESTDIR)${examplesdir}/test_mf_molecfit;\
#	$(DESTDIR)${examplesdir}/test_mf_trans;\
#	$(DESTDIR)${examplesdir}/test_mf_corr
	