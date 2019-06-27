#!/usr/bin/python -B
#
# MENU structure of the MOLECFIT graphical user interface (GUI)
# 
#  Contact: ESO User Support Department usd-help@eso.org
#
#  This file is part of the MOLECFIT software package.
#  Copyright (C) 2009-2013 European Southern Observatory
#
#  This programme is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This programme is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this programme. If not, see <http://www.gnu.org/licenses/>.
#
import os
import subprocess
import sys
import wx
import SM02GUI_Global      as GUI_Global
import SM02GUI_Popups      as GUI_Popups
import SM02GUI_Mask        as GUI_Mask
import SM02GUI_Data        as GUI_Data
import SM02GUI_Paramfile   as GUI_Parameters
import SM02GUI_Help_Dialog_ToolTip_Texts as TextPool
import SM02GUI_LogPopup    as GUI_Console

HELPFILE   = str('Tutorial.pdf')
MANUALFILE = str('Manual.pdf')

#=============== START MENU ====================================================
def OnAbout(e):
    dlg = wx.MessageDialog( GUI_Global.frame, 
		 TextPool.dialog_001,
	     GUI_Global.VERSION, wx.OK)
    dlg.ShowModal() # Show it
    dlg.Destroy()   # finally destroy it when finished with OK

def OnShort(e):
    dlg = wx.MessageDialog( GUI_Global.frame, 
	TextPool.dialog_002,
	     GUI_Global.VERSION, wx.OK)
    dlg.ShowModal() # Show it
    dlg.Destroy()   # finally destroy it when finished with OK

def OnExit(e):
    if GUI_Global.param_changes == True:
        dlg = wx.MessageDialog(GUI_Global.frame, "Parameter changed\n"\
        "Do you want to SAVE before EXIT ?",
         GUI_Global.VERSION, wx.YES_NO | wx.ICON_QUESTION | wx.NO_DEFAULT)
        result = dlg.ShowModal()             # Show it
        dlg.Destroy()                        # finally destroy it when finished 
        if result == wx.ID_YES:
           OnSaveAs(None)  
        
    dlg = wx.MessageDialog(GUI_Global.frame, TextPool.dialog_003,
         GUI_Global.VERSION, wx.YES_NO | wx.ICON_QUESTION | wx.NO_DEFAULT)
    result = dlg.ShowModal()             # Show it
    dlg.Destroy()                        # finally destroy it when finished 
    if result == wx.ID_YES:
        GUI_Global.frame.Close(True)     # Close the frame.
        sys.exit(0)

def OnOpen(e):
    # Create a list of filters
    filters = 'Config Files (*.par)|*.par|All files (*.*)|*.*'
    dialog = wx.FileDialog ( None, message = 
             TextPool.dialog_004, 
             wildcard = filters,
             style = wx.FD_OPEN)
    dialog.SetDirectory(GUI_Global.cwd)
    if dialog.ShowModal() == wx.ID_OK:
        # We'll have to make room for multiple files here
        selected = dialog.GetPath()
        GUI_Parameters.read_param_file(selected)
    else:
        GUI_Parameters.my_message('Nothing was selected!')
    dialog.Destroy()

def OnSave(e):
    if GUI_Parameters.generate_and_write_mask_files() == False:
        dlg  = wx.MessageDialog( GUI_Global.frame, "ERROR: 'Save' failed "\
           "due to missing filenames for masks.\n\nGoto"\
           "\n   Settings->Working Directories and Files\n"\
                "and add filenames.\n",
           "ERROR", wx.OK)
        dlg.ShowModal() 
        dlg.Destroy()           
        return
    if GUI_Parameters.write_paramfile(GUI_Global.frame, \
                             GUI_Global.Current_Setup_File, False) == False:
        return

def OnSaveAs(e):
    if GUI_Parameters.generate_and_write_mask_files() == False:
        dlg  = wx.MessageDialog( GUI_Global.frame, "ERROR: 'Save As' failed "\
            "due to missing filenames for masks.\n\nGoto"\
            "\n   Settings->Working Directories and Files\n"\
                "and add filenames.\n",
            "ERROR", wx.OK)
        dlg.ShowModal() 
        dlg.Destroy()           
        return
    # Create a list of filters
    filters = 'Config Files (*.par)|*.par|All files (*.*)|*.*'
    files   = os.path.basename(GUI_Global.Current_Setup_File)
    dirs    = os.path.dirname(GUI_Global.Current_Setup_File)
    dialog = wx.FileDialog ( None, message = TextPool.dialog_005, 
             wildcard = filters,  defaultFile = files, defaultDir = dirs, 
             style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
    dialog.SetDirectory(GUI_Global.cwd)
    if dialog.ShowModal() == wx.ID_OK:
        # We'll have to make room for multiple files here
        selected = dialog.GetPath()
        GUI_Global.Current_Setup_File = selected
        text = ("Currently used parameter file:           %s" % 
            str(GUI_Global.Current_Setup_File))
        GUI_Global.current_par_file_t.SetLabel(text)
        # Export of the configuration setup to a new selection
        if GUI_Parameters.write_paramfile(GUI_Global.frame, selected, False) \
            == False:
            return
    else:
        dlg  = wx.MessageDialog( GUI_Global.frame, TextPool.dialog_006,
           TextPool.dialog_006, wx.OK)
        dlg.ShowModal() 
        dlg.Destroy()   
    dialog.Destroy()

def OnExport(e):
    # according to the mini workshop at ESO this feature is not used at 
    # the moment 
    # Create a list of filters
    filters = 'FITS Files (*.fits)|*.fits|All files (*.*)|*.*'
    dialog = wx.FileDialog ( None, message = TextPool.dialog_007, 
             wildcard = filters, 
             style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    dialog.SetDirectory(GUI_Global.cwd)

    dialog.Destroy()

def OnDefault(e):
    GUI_Parameters.set_default_parameters()
    GUI_Global.param_available = True
    if e != None:
        dlg  = wx.MessageDialog( GUI_Global.frame, 
           TextPool.dialog_010, "INFO", wx.OK)
        dlg.ShowModal() 
        dlg.Destroy()   

def OnImport(e):
    # Create a list of filters
    filters = 'FITS Files (*.fits)|*.fits|All files (*.*)|*.*'
    dialog = wx.FileDialog ( None, message = TextPool.dialog_009, 
             wildcard = filters, 
             style = wx.FD_OPEN)
    dialog.SetDirectory(GUI_Global.workdir_v)
    if dialog.ShowModal() == wx.ID_OK:
        GUI_Global.filename_v = dialog.GetPath()
    dialog.Destroy()
    

def mime_open(filename):
    """ open file with default application """
    if sys.platform.startswith('darwin'):
        subprocess.call(('open', filename))
    elif os.name == 'nt':
        os.startfile(filename)
    elif os.name == 'posix':
        subprocess.call(('xdg-open', filename))

def OnHelp(e):
    mime_open(os.path.join(GUI_Global.basedir_v, 'share/doc/molecfit', HELPFILE));

def OnManual(e):
    mime_open(os.path.join(GUI_Global.basedir_v, 'share/doc/molecfit', MANUALFILE));

def create_menu(parent):
    # Creating the menubar itself
    menuBar = wx.MenuBar()

	#---------------------------------------------------------------------------
    # Setting up the menu the FIRST pulldown.
    filemenu = wx.Menu()
    menuOpen  = filemenu.Append(wx.NewId(),TextPool.menu_002,TextPool.tip_002)
    menuSave  = filemenu.Append(wx.NewId(),TextPool.menu_001,TextPool.tip_001)
    menuSaveAs= filemenu.Append(wx.NewId(),TextPool.menu_001b,TextPool.tip_001b)
    filemenu.AppendSeparator()
    menuDefault= filemenu.Append(wx.NewId(), TextPool.menu_02c, TextPool.tip_02c)
    # according to the mini workshop at ESO June 2013  this feature is not 
    # used at the moment but kept in code
    #menuImport = filemenu.Append(wx.NewId(), TextPool.menu_02b, TextPool.tip_02b)
    #menuExport= filemenu.Append(wx.NewId(), TextPool.menu_003, TextPool.tip_003)
    filemenu.AppendSeparator()
    menuExit = filemenu.Append(wx.NewId(), TextPool.menu_004, TextPool.tip_004)
    # Adding the "filemenu" to the MenuBar
    menuBar.Append(filemenu,"File\t") 
	# Adding actions
    parent.Bind(wx.EVT_MENU, OnOpen,   menuOpen)
    parent.Bind(wx.EVT_MENU, OnSave,   menuSave)
    parent.Bind(wx.EVT_MENU, OnSaveAs, menuSaveAs)
    parent.Bind(wx.EVT_MENU, OnDefault,menuDefault)
    # parent.Bind(wx.EVT_MENU, OnImport, menuImport)
    # parent.Bind(wx.EVT_MENU, OnExport, menuExport)
    parent.Bind(wx.EVT_MENU, OnExit,   menuExit)

	#---------------------------------------------------------------------------
    # Setting up the menu the SECOND pulldown.
    settingsmenu = wx.Menu()
    menuExcl = settingsmenu.Append(wx.NewId(), TextPool.menu_005, 
            TextPool.tip_005)

    # Adding the "settings menu" to the MenuBar
    menuBar.Append(settingsmenu,"Masking\t") 
	# Adding actions
    parent.Bind(wx.EVT_MENU, GUI_Mask.show_popup_mask, menuExcl)

	#---------------------------------------------------------------------------
    # Setting up the menu the THIRD pulldown.
    settingsmenu = wx.Menu()
    # my own private ID and thus no 'wx.' in front
    menuMol =settingsmenu.Append(wx.NewId(),TextPool.menu_006, TextPool.tip_006)
    menuWave=settingsmenu.Append(wx.NewId(),TextPool.menu_007, TextPool.tip_007)
    menuRes =settingsmenu.Append(wx.NewId(),TextPool.menu_008, TextPool.tip_008)
    menuCont=settingsmenu.Append(wx.NewId(),TextPool.menu_009, TextPool.tip_009)
    menuBack=settingsmenu.Append(wx.NewId(),TextPool.menu_010, TextPool.tip_010)
    settingsmenu.AppendSeparator()
    menuPath=settingsmenu.Append(wx.NewId(),TextPool.menu_011, TextPool.tip_011)
    menuOther= settingsmenu.Append(wx.NewId(), TextPool.menu_012, 
        TextPool.tip_012)

    # Adding the "settings menu" to the MenuBar
    menuBar.Append(settingsmenu,"Settings\t") 
	# Adding actions
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_molec,       menuMol)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_wavelength,  menuWave)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_kernel,      menuRes)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_continuum,   menuCont)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_back,        menuBack)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_path,        menuPath)
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_other,       menuOther)

	#---------------------------------------------------------------------------
    # Setting up the menu the FOURTH pulldown.
    runmenu = wx.Menu()
    # my own private ID and thus no 'wx.' in front
    menuRun =runmenu.Append(wx.NewId(),TextPool.menu_013, TextPool.tip_013)
    runmenu.AppendSeparator()
    menuConsoleOn  =runmenu.Append(wx.NewId(),TextPool.menu_013a, 
                        TextPool.tip_013a)
    menuConsoleOff =runmenu.Append(wx.NewId(),TextPool.menu_013b, 
                        TextPool.tip_013b)

    # Adding the "settings menu" to the MenuBar
    menuBar.Append(runmenu,"RUN\t") 

	# Adding actions
    parent.Bind(wx.EVT_MENU, GUI_Popups.show_popup_fit,   menuRun)
	# Adding actions
    parent.Bind(wx.EVT_MENU, GUI_Console.unhide_console,  menuConsoleOn)
	# Adding actions
    parent.Bind(wx.EVT_MENU, GUI_Console.hide_console,    menuConsoleOff)

	#---------------------------------------------------------------------------
    # Setting up the menu the final help and about pulldown.
    helpmenu  = wx.Menu()
    menuShort  = helpmenu.Append(wx.NewId(),TextPool.menu_014, TextPool.tip_014)
    menuHelp   = helpmenu.Append(wx.NewId(),TextPool.menu_015, TextPool.tip_015)
    menuManual = helpmenu.Append(wx.NewId(),TextPool.menu_016, TextPool.tip_016)
    helpmenu.AppendSeparator()
    menuAbout = helpmenu.Append(wx.NewId(),TextPool.menu_017, TextPool.tip_017)
	# Adding actions
    parent.Bind(wx.EVT_MENU, OnShort, menuShort)
    parent.Bind(wx.EVT_MENU, OnHelp,  menuHelp)
    parent.Bind(wx.EVT_MENU, OnManual,menuManual)
    parent.Bind(wx.EVT_MENU, OnAbout, menuAbout)
    #  Adding the "helpmenu" to the MenuBar
    menuBar.Append(helpmenu,"?") 

	#---------------------------------------------------------------------------
    # Adding the MenuBar to the Frame content.
    parent.SetMenuBar(menuBar)  


