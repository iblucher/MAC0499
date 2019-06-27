
AUTOMAKE_OPTIONS = foreign
DISTCLEANFILES = *.xml

WORKFLOWS = molecfit_reflex.xml

#This installs the workflow and the OCA rules in reflexaux directory
wkfextra_DATA = $(WORKFLOWS)

wkfextra_SCRIPTS = create_parfile.sh

#Pack in the distribution the original workflow
EXTRA_DIST = $(WORKFLOWS).in $(wkfextra_SCRIPTS)

install-data-hook:
	chmod -R uga+rw $(DESTDIR)$(wkfextradir)
	chmod uga+x $(DESTDIR)$(wkfextradir)/create_parfile.sh
	