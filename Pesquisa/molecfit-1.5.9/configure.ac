#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([molecfit], [1.5.9], [usd-help@eso.org])
AC_CONFIG_AUX_DIR([admin])
AC_CONFIG_SRCDIR([src/mf_basic.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/test/Makefile
                 doc/Makefile
                 data/Makefile
                 examples/Makefile
                 examples/CRIRES_4889p5/Barnard.par
                 examples/CRIRES_4889p5/expert_fit.par
                 examples/CRIRES_4889p5/gammaGem.par
                 examples/CRIRES_4889p5/gammaGem_sky.par
                 examples/CRIRES_4889p5/HD134453_expert.par
                 examples/CRIRES_4889p5/HD20010.par
                 examples/CRIRES_HighSNR/CRIRES_HighSNR_Science_AM1p475.par
                 examples/CRIRES_HighSNR/CRIRES_HighSNR_Telluric_AM1p472.par
                 examples/FLAMES/ngc3603.par
                 examples/KMOS/N6822_36.par
                 examples/KMOS/N6822_40.par
                 examples/KMOS/N6822_8.par
                 examples/MUSE/molecfit_5326.par
                 examples/SINFONI/sinfoni_j_whole.par
                 examples/VISIR/molecfit_CERES_H2S2.par
                 examples/VISIR/molecfit_HD104237_H2S2.par
                 examples/XSHOOTER/molecfit_XSHOOTER_NIR_Pipeline_R71.par
                 examples/XSHOOTER/molecfit_XSHOOTER_NIR_TELL_2D.par
                 examples/XSHOOTER/molecfit_XSHOOTER_VIS_ADP_iPTF16abc.par
                 examples/XSHOOTER/molecfit_XSHOOTER_VIS_ADP_V1309Sco.par
                 examples/XSHOOTER/molecfit_XSHOOTER_VIS_HRCar_TELL.par
                 examples/XSHOOTER/molecfit_XSHOOTER_VIS_Pipeline_HIP079048_F9V.par
                 examples/XSHOOTER/molecfit_XSHOOTER_VIS_Pipeline_R71.par
                 examples/XSHOOTER/XSHOOTER_NIR_TELL.lis
                 examples/XSHOOTER/XSHOOTER_VIS_HRCar.lis
                 gui/Makefile
                 reflex/Makefile
                 reflex/molecfit_reflex.xml])
AM_INIT_AUTOMAKE([-Wall foreign])

# Checks for programmes.
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_PROG_RANLIB
AC_PROG_LN_S

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Check for CPL presence and usability
CPL_CHECK_LIBS

AC_CHECK_LIB([m],[pow],[MATHLIB=-lm],[
    AC_MSG_ERROR([Math library (libm) not found!])
])
AC_SUBST([MATHLIB])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h stdio.h ctype.h math.h regex.h],[],[
    AC_MSG_ERROR(["A standard header file was not found!"])
])

# Checks for typedefs, structures, and compiler characteristics.
ESO_PROG_CC_FLAG([std=c99], [CFLAGS="$CFLAGS -std=c99"])
ESO_PROG_CC_FLAG([Werror=implicit-function-declaration], [CFLAGS="$CFLAGS -Werror=implicit-function-declaration"])
ESO_PROG_CC_FLAG([Werror=incompatible-pointer-types], [CFLAGS="$CFLAGS -Werror=incompatible-pointer-types"])
ESO_PROG_CC_FLAG([fvisibility=hidden], [CFLAGS="$CFLAGS -fvisibility=hidden"])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([mkdir strchr strpbrk strrchr strspn strstr strtol])
AM_MAINTAINER_MODE([enable])

AC_CHECK_PROG(SVNVERSION, svnversion, found, not-found)
if test $SVNVERSION = found; then
  SVNREV=$(svnversion "$srcdir")
  echo $SVNREV > SVNREV
else
  SVNREV=$(cat SVNREV)
fi
AC_SUBST(SVNREV)

if test -z "$wkfextradir"; then
    wkfextradir='${prefix}/share/${PACKAGE_NAME}/reflex'
fi

if test -z "$configdir"; then
    configdir='${prefix}/share/${PACKAGE_NAME}/config'
fi

if test -z "$pyguidir"; then
    pyguidir='${prefix}/share/${PACKAGE_NAME}/gui'
fi

if test -z "$examplesdir"; then
    examplesdir='${prefix}/share/${PACKAGE_NAME}/examples'
fi

if test -z "$profilesdir"; then
    profilesdir='${prefix}/share/${PACKAGE_NAME}/data'
fi

OS=$(uname)
if test "$OS" = "Darwin"; then
    installdatahook_sediopt=' .bak'
else
    installdatahook_sediopt=''
fi

AC_ARG_ENABLE(molecfit-gui,
              AS_HELP_STRING([--disable-molecfit-gui],
                             [disable the molecfit gui]),
              molecfit_gui=$enableval,
              molecfit_gui=yes)

AC_SUBST(molecfit_gui)


AC_SUBST(configdir)
AC_SUBST(wkfextradir)
AC_SUBST(pyguidir)
AC_SUBST(examplesdir)
AC_SUBST(profilesdir)
AC_SUBST(installdatahook_sediopt)
AC_SUBST(docdir)

AC_OUTPUT
