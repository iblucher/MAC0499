
AUTOMAKE_OPTIONS = foreign

nobase_profiles_DATA = \
		   profiles/lib/GDAS_t0_s1.atm \
		   profiles/lib/GDAS_t0_s2.atm \
		   profiles/lib/GDAS_t0_s3.atm \
		   profiles/lib/GDAS_t0_s4.atm \
		   profiles/lib/GDAS_t0_s5.atm \
		   profiles/lib/GDAS_t0_s6.atm \
		   profiles/mipas/equ.atm \
		   profiles/mipas/std.atm \
		   profiles/mipas/tro.atm \
		   profiles/gdas/gdas_profiles_C-70.4-24.6.tar.gz

EXTRA_DIST = $(nobase_profiles_DATA)

install-data-hook:
	chmod -R uga+rw $(DESTDIR)$(profilesdir)/
	