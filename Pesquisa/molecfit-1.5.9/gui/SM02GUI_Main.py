#!/usr/bin/env python2.7
#
# MAIN code starting the MOLECFIT graphical user interface (GUI)
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

import matplotlib
# monkey patch matplotlib bug preventing mask tool selection
# fixed in >= 1.1.1, >= 1.2.1 and >= 1.3.0
# https://github.com/matplotlib/matplotlib/pull/1630
def _fixed_disconnect_110(self, cid):
    """
    disconnect the callback registered with callback id *cid*
    """
    for eventname, callbackd in self.callbacks.items():
        try:
            del callbackd[cid]
        except KeyError:
            continue
        else:
            for key, value in self._func_cid_map.items():
                if value == cid: 
                    del self._func_cid_map[key]
            return

def _fixed_disconnect_120(self, cid):
    """
    disconnect the callback registered with callback id *cid*
    """
    for eventname, callbackd in self.callbacks.items():
        try:
            del callbackd[cid]
        except KeyError:
            continue
        else:
            for category, functions in self._func_cid_map.items():
                for function, value in functions.items():
                    if value == cid: 
                        del functions[function]
            return

try:
    if matplotlib.__version__ == "1.1.0":
        matplotlib.cbook.CallbackRegistry.disconnect = _fixed_disconnect_110
    elif matplotlib.__version__ == "1.2.0":
        matplotlib.cbook.CallbackRegistry.disconnect = _fixed_disconnect_120
except:
    pass

matplotlib.use('WXAgg')

import wx
import os
import sys
import SM02GUI_Global      as GUI_Global
import SM02GUI_Buttons     as GUI_Button
import SM02GUI_Init        as GUI_Init
import SM02GUI_Menu        as GUI_Menu
import SM02GUI_Plot        as GUI_Plot
import SM02GUI_Paramfile   as GUI_Parameters
import SM02GUI_Popups      as GUI_Popups
import SM02GUI_Data        as GUI_Data
import SM02GUI_Mask        as GUI_Mask
import SM02GUI_Help        as GUI_Help
import SM02GUI_LogPopup    as GUI_Console

#===============================================================================
# the main window class of the GUI
class MainWindow(wx.Frame):
    
    def __init__(self, parent, title):
        
        GUI_Global.VERSION = 'Molecfit GUI v1.5.9'
        # open a frame in the main window
        # the size and the pos is given and used on screen if the 
        # screen and the window manager allows it. E.g. Gnome allows
        # no position above 27 due to its headline and if window is oversized 
        # at start time the maximum of the screen in that very direction is used
        
        my_frame = wx.Frame.__init__(self, parent, title=GUI_Global.VERSION, 
            pos=(40,40), size=(1300,800), style= wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | 
            wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX)
        
        GUI_Global.frame = my_frame

        if GUI_Init.init_path_and_environment() == False:
            sys.exit(-1)

        # for creation of the popup masks the molecule structure has to exist 
        # first as it defines some deafult class arrays for the molecules
        GUI_Parameters.create_and_fill_default_molecules()
        self.CreateStatusBar() 
        # create empty popup masks 
        GUI_Popups.create_popup_other(self)
        GUI_Popups.create_popup_wavelength(self)
        GUI_Popups.create_popup_background(self)
        GUI_Popups.create_popup_continuum(self)
        GUI_Popups.create_popup_kernel(self)
        GUI_Popups.create_popup_molec(self)
        GUI_Popups.create_popup_molec_sub(self)
        GUI_Help.create_popup_error(self)
        GUI_Help.create_popup_help(self)
        GUI_Popups.create_popup_path(self)
        # create a Menue and the main Buttons for the main window
        GUI_Menu.create_menu(self)
        GUI_Button.create_action_buttons(self)
        # create the windows for the submenues - but do NOT show them yet
        GUI_Mask.create_popup_mask(self)
        GUI_Popups.create_popup_fit(self)
        GUI_Popups.create_popup_apply(self)
        #GUI_Popups.create_popup_parameter(self)

        # or some dummy to be able to define mathplotlib plots
        GUI_Data.empty_data_define()         
        GUI_Plot.create_main_panel(self)
        GUI_Plot.draw_figure(False)

        GUI_Console.create_console(self)
        GUI_Menu.OnDefault(None)

        # interpretation of the first parameter (if given) as a config file 
        if len(sys.argv) > 1:
            path = GUI_Init.make_absolute_path(sys.argv[1])
            if GUI_Parameters.read_param_file(path) == False:
                dlg = wx.MessageDialog( GUI_Global.frame, 
		            "ERROR: parameter error - GUI not started",
	                GUI_Global.VERSION, wx.OK)
                dlg.ShowModal() 
                sys.exit(-1)    
        
        # show it
        self.Show(True) 
        GUI_Global.param_changes = False
        

def setup_library_path(base=None):
    # TODO should be called when basedir parameter changes
    # default to installation root
    if base is None:
        base = os.path.dirname(os.path.dirname(sys.argv[0]))
    libdir = os.path.abspath(os.path.join(base, 'lib/'))
    if sys.platform.startswith('darwin'):
        key = 'DYLD_LIBRARY_PATH'
    else:
        key = 'LD_LIBRARY_PATH'
    if key in os.environ:
        os.environ[key] = ':'.join([libdir, os.environ[key]])
    else:
        os.environ[key] = libdir

setup_library_path()
#  and then start everything defined above
app   = wx.App(False)
try:
    GUI_Global.frame = MainWindow(None, 'SM02 GUI - Molecfit')
    app.MainLoop()
except Exception as e:
    wx.MessageBox(str(e), 'Error', wx.OK | wx.ICON_ERROR)

