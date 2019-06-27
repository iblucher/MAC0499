#!/usr/bin/python -B
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
import wx
import SM02GUI_Global      as GUI_Global

#===============================================================================
#
#  the WARNING/ERROR sub-popup
#

def OnButton_error_done(e):
    GUI_Global.popup_error.Show(False)

def create_popup_error(parent):
    GUI_Global.popup_error = wx.Dialog(parent, wx.NewId(), 'WARNING / ERROR', 
        size=(300, 200))
    
    GUI_Global.error_text = wx.StaticText(GUI_Global.popup_error, 
        -1, 'no error', pos=(20,10))#, color='#ff0000')
    GUI_Global.error_text.SetForegroundColour((255,0,0)) # set text color
    GUI_Global.error_text.SetBackgroundColour((0,0,255)) # set text back color
    button_error_done = wx.Button(GUI_Global.popup_error, id=-1, 
        label="DISMISS", pos=(100,120))
    GUI_Global.popup_error.Bind(wx.EVT_BUTTON, OnButton_error_done, 
        button_error_done ) 

def refresh_error(text):
    GUI_Global.error_text.SetLabel(text)
    GUI_Global.popup_error.Show(True)

#===============================================================================
#
#  the HELP sub-popup
#

def OnButton_help_done(e):
    GUI_Global.popup_help.Show(False)

def create_popup_help(parent):
    GUI_Global.popup_help = wx.Dialog(parent, wx.NewId(), 'HELP', 
        size=(500, 600))
    
    GUI_Global.help_text = wx.StaticText(GUI_Global.popup_help, 
        -1, 'no help', pos=(20,10))#, color='#ff0000')
    GUI_Global.help_text.SetForegroundColour((0,0,255)) # set text color
    GUI_Global.help_text.SetBackgroundColour((0,0,255)) # set text back color
    button_help_done = wx.Button(GUI_Global.popup_help, id=-1, 
        label="DISMISS", pos=(200,550))
    GUI_Global.popup_help.Bind(wx.EVT_BUTTON, OnButton_help_done, 
        button_help_done ) 

def refresh_help(text):
    GUI_Global.help_text.SetLabel(text)
    GUI_Global.popup_help.Show(True)

