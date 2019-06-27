#!/usr/bin/python -B
#
# popup CONSOLE for display of data from underlying codes like LBLRTM 
# or molecfit
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
import SM02GUI_Paramfile   as GUI_Parameters
import SM02GUI_Data        as GUI_Data
#===============================================================================
# the main window class of the GUI

def refresh_my_console():
    # timer controlled refresh feedback function
    array = ''
    if os.path.exists(GUI_Global.console_text_file) == True:
        if  GUI_Global.console_text_file[\
            len(GUI_Global.console_text_file)-5:] == '.diff' or \
            GUI_Global.console_text_file[\
            len(GUI_Global.console_text_file)-10:] == '.temporary':
            if  GUI_Global.console_text_file[\
                len(GUI_Global.console_text_file)-10:] == '.temporary':
                array += 'GUI CONFIG:  current setup in local memory of GUI\n'
                array += '====================================='\
                         '=================================\n'
        else:
            array += 'FILE: ' + GUI_Global.console_text_file + '\n'
            array += '====================================='\
                         '=================================\n'
        fd = open( GUI_Global.console_text_file, 'r' )
        for line in fd:
            array += line.decode('UTF-8')
        fd.close()
        if len(array) == 0:
            array += str('** empty **')
            if GUI_Global.console_text_file[
                len(GUI_Global.console_text_file)-5:] == '.diff':
                array += str(' no difference was found ')
        
    # refresh the text
    GUI_Global.my_console_text.SetValue(unicode(array))
    # get position and set it to the end
    GUI_Global.my_console_text.ShowPosition(
        GUI_Global.my_console_text.GetLastPosition ())
    # set to END my_console_text.SetPosition
    if GUI_Global.reshoot_console == True:
        wx.FutureCall(500, refresh_my_console)

def create_console(parent):
    # create the CONSOLE

    # ----- BUG ----- #
    # The following code makes the dialog on focus.
    # It stays on focus event when the user click on "close"
    # since this button just *hide* the window, it does *NOT* close it.
    #GUI_Global.my_console = wx.Dialog(parent, wx.NewId(), 'Molecfit GUI CONSOLE',
    #                                  size=(650, 730))
    # ----- BUG ----- #
    GUI_Global.my_console = wx.Dialog(None, wx.NewId(), 'Molecfit GUI CONSOLE',
                                      size=(650, 760), style=wx.RESIZE_BORDER | wx.CLOSE_BOX | wx.DEFAULT_DIALOG_STYLE | wx.DIALOG_NO_PARENT)
    GUI_Global.my_console_text = wx.TextCtrl(GUI_Global.my_console,id=-1, style=
        wx.TE_MULTILINE | wx.HSCROLL | wx.TE_READONLY, size=(610,520),
        value=' ', pos=(20,20))
    
    font1 = wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    
    """
    font1 = wx.Font( 11, 
                     family    = wx.FONTFAMILY_TELETYPE,
                     style     = wx.FONTSTYLE_NORMAL,
                     weight    = wx.FONTWEIGHT_NORMAL,
                     encoding  = wx.FONTENCODING_DEFAULT, 
                     underline = False,
                     face='Monospace')
    """
    
    GUI_Global.my_console_text.SetFont(font1)
    wx.StaticText(GUI_Global.my_console, -1, 'logfiles:', pos=(20,554))


    #---------------------------------------------------------------------------
    # DONE BUTTON   
    GUI_Global.my_console.button_done = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="close", pos=(30,680))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, hide_console, 
        GUI_Global.my_console.button_done ) 
    GUI_Global.my_console.button_done.SetToolTip(wx.ToolTip(
        "click to if you are ready"))
    GUI_Global.my_console.button_done.Show(True)

    #---------------------------------------------------------------------------
    # FIT RESULT BUTTON   
    GUI_Global.my_console.button_fit_result = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="FIT RESULT", pos=(150,550))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_fit_console, 
        GUI_Global.my_console.button_fit_result ) 
    GUI_Global.my_console.button_fit_result.SetToolTip(wx.ToolTip(
        "click to show the last fit results (statistics and resulting "\
        "abundances of molecules)"))
    GUI_Global.my_console.button_fit_result.Show(True)

    #---------------------------------------------------------------------------
    # LOG of molecfit/LBLRTM   in fit mode
    GUI_Global.my_console.button_log = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="last FIT log", pos=(260,550))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_log_console, 
        GUI_Global.my_console.button_log ) 
    GUI_Global.my_console.button_log.SetToolTip(wx.ToolTip(
        "click to show the last logfile of a FIT action"))
    GUI_Global.my_console.button_log.Show(True)

    #---------------------------------------------------------------------------
    # LOG of molecfit/LBLRTM   in apply mode 
    GUI_Global.my_console.button_apply = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="last APPLY log", pos=(370,550))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_apply_console, 
        GUI_Global.my_console.button_apply ) 
    GUI_Global.my_console.button_apply.SetToolTip(wx.ToolTip(
        "click to show the last logfile of an APPLY action"))
    GUI_Global.my_console.button_apply.Show(True)

    #---------------------------------------------------------------------------
    # LOG of molecfit   in prepguitable
    GUI_Global.my_console.button_prep = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="prepguitable log", pos=(500,550))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_preguitable_console, 
        GUI_Global.my_console.button_prep ) 
    GUI_Global.my_console.button_prep.SetToolTip(wx.ToolTip(
        "click to show the last logfile of a prepguitable call "\
        "(used to load FITS files)"))
    GUI_Global.my_console.button_prep.Show(True)

    #---------------------------------------------------------------------------
    # display CONFIG file   
    wx.StaticLine(GUI_Global.my_console,-1, (20,595), (520,0) )
    wx.StaticText(GUI_Global.my_console, -1, 'parameters:   '
        'stored on current param. file             vs.                       '
        'local setting in GUI', pos=(20,604))
    GUI_Global.my_console.button_config = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="stored CONFIG", pos=(150,640))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_config_console, 
        GUI_Global.my_console.button_config ) 
    GUI_Global.my_console.button_config.SetToolTip(wx.ToolTip(
        "click to show the current stored configuration"))
    GUI_Global.my_console.button_config.Show(True)

    #---------------------------------------------------------------------------
    # display CONFIG setting currently in RAM of the GUI in same format
    # as if it would be saved/stored  
    GUI_Global.my_console.button_curconfig = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="current GUI config", pos=(490,640))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_curconfig_console, 
        GUI_Global.my_console.button_curconfig ) 
    GUI_Global.my_console.button_curconfig.SetToolTip(wx.ToolTip(
        "click to show the current parameter settings in the GUI"))
    GUI_Global.my_console.button_curconfig.Show(True)

    #---------------------------------------------------------------------------
    # display CONFIG difference between settings on DISK and in RAM
    GUI_Global.my_console.button_curconfig = wx.Button(
        GUI_Global.my_console, 
        id=-1, label="difference", pos=(320,640))
    GUI_Global.my_console.Bind(wx.EVT_BUTTON, show_diffconfig_console, 
        GUI_Global.my_console.button_curconfig ) 
    GUI_Global.my_console.button_curconfig.SetToolTip(wx.ToolTip(
        "click to show the diffrenece between current parameter settings"\
        " and the stored values in the config file"))
    GUI_Global.my_console.button_curconfig.Show(True)


def show_console(filename, respawn):
    # show the file with 'filename' in console - the boolean respawn controlls 
    # if a periodic refresh is done (e.g. while running FIT)
    GUI_Global.console_text_file = filename
    GUI_Global.reshoot_console   = respawn
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_fit_console(e):
    # display the FIT result in the CONSOLE
    GUI_Global.console_text_file = str(GUI_Data.output_dir_path() + '/' \
            + GUI_Global.output_name_v + 
                '_fit.res')
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_config_console(e):
    # display the disk stored configuration
    GUI_Global.console_text_file = GUI_Global.Current_Setup_File
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_curconfig_console(e):
    # display the current 'in RAM' stored configuration
    GUI_Parameters.write_paramfile(GUI_Global.frame,
        GUI_Global.Current_Setup_File+'.temporary', True)
    GUI_Global.console_text_file = GUI_Global.Current_Setup_File+'.temporary'
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_diffconfig_console(e):
    # display the difference of the current 'in RAM' stored configuration
    # with that on DISK
    GUI_Parameters.write_paramfile(GUI_Global.frame,
        GUI_Global.Current_Setup_File+'.temporary', True)
    cmd = 'diff ' + GUI_Global.Current_Setup_File + ' ' + \
          GUI_Global.Current_Setup_File+'.temporary' + \
          ' > ' + GUI_Global.Current_Setup_File+'.diff'
    os.system(cmd)
    GUI_Global.console_text_file = GUI_Global.Current_Setup_File+'.diff'
    #GUI_Global.console_text_file = filename
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_log_console(e):
    # display the last stored molcfit log
    GUI_Global.console_text_file = str('%s/%s_molecfit_out.txt' % (
        GUI_Data.output_dir_path(),GUI_Global.output_name_v))
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_preguitable_console(e):
    # display the last stored prepguitable log
    GUI_Global.console_text_file = str('%s/%s_prepguitable_out.txt' % (
        GUI_Data.output_dir_path(),GUI_Global.output_name_v))
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def show_apply_console(e):
    # display the last stored molcfit 'apply teluric fit' log
    GUI_Global.console_text_file = str('%s/%s_apply_out.txt' % (
        GUI_Data.output_dir_path(),GUI_Global.output_name_v))
    GUI_Global.reshoot_console   = False
    wx.FutureCall(500, refresh_my_console)
    GUI_Global.my_console.Show()

def hide_console(e):
    # hide console an delete possible refresh respan timers
    GUI_Global.reshoot_console   = False
    GUI_Global.my_console.Show(False)

def unhide_console(e):
    # unhide console
    GUI_Global.reshoot_console   = False
    GUI_Global.my_console.Show(True)

def stop_refresh_console():
    # stop (possibly working) timer refreshes
    GUI_Global.reshoot_console   = False


