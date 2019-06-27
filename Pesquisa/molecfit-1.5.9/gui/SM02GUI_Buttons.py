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
import wx

import SM02GUI_Global      as GUI_Global
import SM02GUI_Popups      as GUI_Popups
import SM02GUI_Mask        as GUI_Mask
import SM02GUI_Menu        as GUI_Menu


def create_action_buttons(parent):
    
    # generations of the buttons and callbacks of the main window
    wx.StaticText(parent,-1,"Step 1", (30,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (90,30), (80,0) )
    button_mask_tool = wx.Button(parent, id=-1, label="Mask Tool", pos=(20,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Mask.show_popup_mask,  button_mask_tool) 
    button_mask_tool.SetToolTip(wx.ToolTip("the masking tool to select included and excluded regions"))
    button_mask_tool.Show(True)

    wx.StaticText(parent,-1,"Step 2", (180,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (240,30), (80,0) )
    button_molecules = wx.Button(parent, id=-1, label="Molecules", pos=(170,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_molec,  button_molecules) 
    button_molecules.SetToolTip(wx.ToolTip("the selector box for molecules"))
    button_molecules.Show(True)

    wx.StaticText(parent,-1,"Step 3", (330,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (390,30), (80,0) )
    button_wave = wx.Button(parent, id=-1, label="Wavelength\n  Correction", 
    pos=(320,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_wavelength,  button_wave) 
    button_wave.SetToolTip(wx.ToolTip("finetune the wavelength grid"))
    button_wave.Show(True)

    wx.StaticText(parent,-1,"Step 4", (480,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (540,30), (80,0) )
    button_resol = wx.Button(parent, id=-1, label="  Resolution  /\nFitting Kernel", pos=(470,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_kernel,  button_resol) 
    button_resol.SetToolTip(wx.ToolTip("the definition of the line spread function (LFS)"))
    button_resol.Show(True)

    wx.StaticText(parent,-1,"Step 5", (630,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (690,30), (80,0) )
    button_cont = wx.Button(parent, id=-1, label="Continuum", pos=(620,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_continuum,  button_cont) 
    button_cont.SetToolTip(wx.ToolTip("the way how to treat the continuum"))
    button_cont.Show(True)

    wx.StaticText(parent,-1,"Step 6", (780,15), (-1,-1), wx.ALIGN_LEFT)
    wx.StaticLine(parent,-1, (840,30), (150,0) )
    button_back = wx.Button(parent, id=-1, label="Background", pos=(770,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_back,  button_back) 
    button_back.SetToolTip(wx.ToolTip("the telescope background / emissivity"))
    button_back.Show(True)

    wx.StaticText(parent,-1,"Step 7", (1000,15), (-1,-1), wx.ALIGN_LEFT)
    button_fit = wx.Button(parent, id=-1, label="FIT", pos=(990,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_fit,  button_fit) 
    button_fit.SetToolTip(wx.ToolTip("finally execute a fit"))
    button_fit.Show(True)

    button_apply = wx.Button(parent, id=-1, label="Apply telluric\n  Correction", pos=(1150,40), size=(110,60))
    parent.Bind(wx.EVT_BUTTON, GUI_Popups.show_popup_apply,  button_apply) 
    button_apply.SetToolTip(wx.ToolTip("and apply the correction to the data"))
    button_apply.Show(True)

    button_exit = wx.Button(parent, id=-1, label="EXIT", pos=(1150,650), size=(110,40))
    parent.Bind(wx.EVT_BUTTON, GUI_Menu.OnExit,  button_exit) 
    button_exit.SetToolTip(wx.ToolTip("exit the GUI"))
    button_exit.Show(True)

    GUI_Global.current_plot_file_t = wx.StaticText(parent, -1, 
        "Currently displayed working file:     -- None --", pos=(120,640), size=(-1,-1))
    GUI_Global.current_plot_file_t.SetForegroundColour((0,0,255))
    
    GUI_Global.current_par_file_t = wx.StaticText(parent, -1, 
        "Currently used parameter file:           *** not defined yet ***", pos=(120,660), size=(-1,-1))
    GUI_Global.current_par_file_t.SetForegroundColour((0,0,255))
    
    GUI_Global.current_rescale_t = wx.StaticText(parent, -1, 
        ('Data rescale factor for display/plot: %2.2g' % (GUI_Global.rescale+GUI_Global.rescale*1e-10)), 
        pos=(120,680), size=(-1,-1))
    GUI_Global.current_rescale_t.Show(False)
    GUI_Global.current_rescale_t.SetForegroundColour((255,0,0))
  
