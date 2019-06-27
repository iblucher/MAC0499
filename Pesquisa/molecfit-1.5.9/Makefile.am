AUTOMAKE_OPTIONS = 1.8 foreign

ACLOCAL_AMFLAGS = -I m4macros

DISTCLEANFILES = *~ SVNREV

EXTRA_DIST = m4macros/eso.m4 m4macros/cpl.m4 SVNREV

dist_bin_SCRIPTS = bin/update_gdas_database.sh \
		   bin/extract_grib.sh \
		   bin/get_gdas_profiles.sh

dist_noinst_SCRIPTS = scripts/install_package.sh \
		      scripts/molecfit_installer.run

dist_config_DATA = config/lblrtm_setup


SUBDIRS = src reflex gui examples data doc

# Files to be removed with command: make maintainer-clean
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in config.log \
        config.status configure stamp-h.in depcomp install-sh missing

install-data-hook:
	chmod -R uga+rw $(DESTDIR)$(configdir)
	