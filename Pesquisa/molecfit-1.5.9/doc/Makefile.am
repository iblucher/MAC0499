
AUTOMAKE_OPTIONS = foreign

# temp, precompiled manual
doc_DATA = VLT-MAN-ESO-19550-5772_Molecfit_User_Manual.pdf \
           VLT-MAN-ESO-19550-5928_Molecfit_GUI_and_Tutorial.pdf

# find figures/ -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> subfiles.mk
include $(srcdir)/subfiles.mk
dist_noinst_DATA = dmd-doc.sty \
                   eso-logo.eps \
                   MF_UM_acronyms.tex \
                   MF_UM_app_cheatsheets.tex \
                   MF_UM_appendix.tex \
                   MF_UM_evaluation.tex \
                   MF_UM_installation.tex \
                   MF_UM_introduction.tex \
                   MF_UM_layout.tex \
                   MF_UM_model.tex \
                   MF_UM_running.tex \
                   MF_UM_shortcut.tex \
                   MF_UM_skipped.tex \
                   molecfit_user_manual.tex \
                   VLT-MAN-ESO-19550-5928_Molecfit_GUI_and_Tutorial.docx \
                   $(subdir_files)

install-data-hook:
	ln -sf VLT-MAN-ESO-19550-5928_Molecfit_GUI_and_Tutorial.pdf $(DESTDIR)$(docdir)/Tutorial.pdf
	ln -sf VLT-MAN-ESO-19550-5772_Molecfit_User_Manual.pdf $(DESTDIR)$(docdir)/Manual.pdf

manual: $(dist_noinst_DATA) $(subdir_files)
	pdflatex molecfit_user_manual.tex
	pdflatex molecfit_user_manual.tex

CLEANFILES = rm -f *.dvi *.aux *.toc *.out *.log *.idx *.ps *.bbl *.blg *.brf \
	     molecfit_user_manual.pdf \
	     figures/*converted-to*.pdf \
	     *converted-to*.pdf

EXTRA_DIST = $(doc_DATA) $(dist_noinst_DATA)

.NOTPARALLEL:
