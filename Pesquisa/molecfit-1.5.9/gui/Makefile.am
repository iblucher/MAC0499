
AUTOMAKE_OPTIONS = foreign
DISTCLEANFILES = *.py[oc]

pygui_DATA = SM02GUI_Buttons.py \
	     SM02GUI_Data.py \
	     SM02GUI_Global.py \
	     SM02GUI_Help_Dialog_ToolTip_Texts.py \
	     SM02GUI_Help.py \
	     SM02GUI_Init.py \
	     SM02GUI_LogPopup.py \
	     SM02GUI_Main.py \
	     SM02GUI_Mask.py \
	     SM02GUI_Menu.py \
	     SM02GUI_Paramfile.py \
	     SM02GUI_Plot.py \
	     SM02GUI_Popups.py \
	     SM02GUI_FitKernelDialog.py


EXTRA_DIST = $(pygui_DATA) SM02GUI_FitKernelDialog.fbp

install-data-hook:
	if [ $(molecfit_gui) = yes ]; then\
		echo "Install molecfit_gui";\
		ln -sf $(pyguidir)/SM02GUI_Main.py $(DESTDIR)$(bindir)/molecfit_gui;\
		chmod -R uga+rw $(DESTDIR)$(pyguidir);\
		chmod uga+x $(DESTDIR)$(pyguidir)/SM02GUI_Main.py;\
	else\
		echo "Remove molecfit_gui";\
		rm -rf $(DESTDIR)$(pyguidir);\
	fi
