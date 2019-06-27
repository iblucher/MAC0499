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
import os
import subprocess
from functools import partial
#from subprocess import Popen, PIPE
#import time
import SM02GUI_Global      as GUI_Global
import SM02GUI_Menu        as GUI_Menu
import SM02GUI_Plot        as GUI_Plot
import SM02GUI_Paramfile   as GUI_Parameters
#import SM02GUI_Status_Text as GUI_Status_Text
import SM02GUI_Help        as GUI_Help
import SM02GUI_Data        as GUI_Data
import SM02GUI_Help_Dialog_ToolTip_Texts as TextPool
import SM02GUI_LogPopup    as GUI_Console
from SM02GUI_FitKernelDialog import FitKernelDialog

def EndBusyCursor():
    """ wxpython 3 asserts if its not set """
    try:
        wx.EndBusyCursor()
    except:
        pass

class Validator(wx.PyValidator):
    def __init__(self, base, name, value_check=None):
        wx.PyValidator.__init__(self)
        self.base = base
        self.name = name
        if value_check is not None:
            self.value_check = value_check
        else:
            self.value_check = lambda x: x

    def Clone(self):
        return self.__class__(base=self.base, name=self.name,
                              value_check=self.value_check)

    def convert(self, value):
        return value

    def Validate(self, win):
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        try:
            self.value = self.convert(text)
        except ValueError as e:
           textCtrl.SetBackgroundColour("pink")
           textCtrl.SetFocus()
           textCtrl.Refresh()
           wx.MessageBox(str(e), "Error")
           return False

        textCtrl.SetBackgroundColour(
            wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        textCtrl.Refresh()
        return True

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        value = self.value
        if value != getattr(self.base, self.name):
            setattr(self.base, self.name, value)
            GUI_Global.param_changes = True
        return True

class IntValidator(Validator):
    def convert(self, value):
        try:
            v = int(value)
        except ValueError:
            raise ValueError("Integer expected for marked control")
        return self.value_check(v)

class FloatValidator(Validator):
    def convert(self, value):
        try:
            v = float(value)
        except ValueError:
            raise ValueError("Floating point number expected for marked control")
        return self.value_check(v)


class BoolValidator(wx.PyValidator):
    """ no-op validator that sets base.name """
    def __init__(self, base, name):
        wx.PyValidator.__init__(self)
        self.base = base
        self.name = name

    def Clone(self):
        return self.__class__(base=self.base, name=self.name)

    def Validate(self, win):
        return True

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        boolCtrl = self.GetWindow()
        value = boolCtrl.IsChecked()
        if value != getattr(self.base, self.name):
            setattr(self.base, self.name, value)
            GUI_Global.param_changes = True
        return True


class ListValidator(wx.PyValidator):
    """ validator that sets base.name to value if its control is True """
    def __init__(self, base, name, value):
        wx.PyValidator.__init__(self)
        self.base = base
        self.name = name
        self.value = value

    def Clone(self):
        return self.__class__(base=self.base, name=self.name, value=self.value)

    def Validate(self, win):
        return True

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        v = self.GetWindow().GetValue()
        if v and self.value != getattr(self.base, self.name):
            setattr(self.base, self.name, self.value)
            GUI_Global.param_changes = True
        return True


def g_zero(value):
    if value <= 0:
        raise ValueError("Value must be greater zero")
    return value

def ge_zero(value):
    if value < 0:
        raise ValueError("Value must be greater or equal zero")
    return value

def in_range(value, low, high, lowincl=True, highincl=True):
    if lowincl and highincl:
        if value < low or value > high:
            raise ValueError("Value must be within [%r, %r]" % (low, high))
    elif not lowincl and highincl:
        if value <= low or value > high:
            raise ValueError("Value must be within ]%r, %r]" % (low, high))
    elif lowincl and not highincl:
        if value < low or value >= high:
            raise ValueError("Value must be within [%r, %r[" % (low, high))
    else:
        if value <= low or value >= high:
            raise ValueError("Value must be within ]%r, %r[" % (low, high))
    return value

def b_zero_one(value):
    return in_range(value, low=0, high=1, lowincl=False, highincl=True)


#===============================================================================
#
#  routines to unmap / hide the popups
#

def OnButton_help(e):
    # one feedback system for all popups is enough now
    if (GUI_Global.popup_molec.button_help == e.GetEventObject()):
        GUI_Help.refresh_help(TextPool.help_001) 
    if (GUI_Global.popup_molec_sub.button_help == e.GetEventObject()):
        GUI_Help.refresh_help(TextPool.help_002) 
    if (GUI_Global.popup_wavelength.button_help == e.GetEventObject()):
        GUI_Help.refresh_help(TextPool.help_003) 

def reset_textcontrol_colors(dialog):
    """ reset text control colors,
    validators do not reset after invalid entry followed by cancel
    """
    dflt = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
    # only get our own attributes, accessing Handle causes gdk errors
    members = set(dir(dialog)) - set(dir(wx.Dialog))
    for m in members:
        a = getattr(dialog, m)
        if isinstance(a, wx.TextCtrl):
            a.SetBackgroundColour(dflt)


#===============================================================================
#
#  POPUP for the CONTINUUM
#
def refresh_popup_cont():
    if GUI_Global.popup_continuum.fit_cont_v != 0:
        GUI_Global.popup_continuum.check_fit.SetValue(True)
    else:
        GUI_Global.popup_continuum.check_fit.SetValue(False)
    text = repr(GUI_Global.popup_continuum.cont_const_v)
    GUI_Global.popup_continuum.cont_const_t.SetValue(text)
    text = "%i" % (GUI_Global.popup_continuum.cont_n_v)
    GUI_Global.popup_continuum.cont_n_t.SetValue(text)
    reset_textcontrol_colors(GUI_Global.popup_continuum)

def show_popup_continuum(parent):
    if GUI_Global.param_available == False:
        no_parameter_message()
        return
    refresh_popup_cont()
    GUI_Global.popup_continuum.Show(True)

def create_popup_continuum(parent):
    
    GUI_Global.popup_continuum = wx.Dialog(parent, wx.NewId(), 'CONTINUUM', 
                                      size=(270, 220))
    self = GUI_Global.popup_continuum
    sizer = wx.BoxSizer( wx.VERTICAL )
    
    # TODO layout full dialog with sizer instead of using this spacer
    sizer.AddSpacer(140)

    buttons = wx.StdDialogButtonSizer()
    self.buttonsOK = wx.Button(self, wx.ID_OK)
    buttons.AddButton(self.buttonsOK)
    self.buttonsCancel = wx.Button(self, wx.ID_CANCEL)
    buttons.AddButton(self.buttonsCancel)
    buttons.Realize();
    sizer.Add(buttons, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5)
    self.SetSizer(sizer)
    self.Layout()

    # generate checkbox and set default / initial value
    GUI_Global.popup_continuum.check_fit = wx.CheckBox(
        GUI_Global.popup_continuum, -1, "", pos=(30,20), style=wx.ALIGN_RIGHT,
        validator=BoolValidator(base=GUI_Global.popup_continuum,
                                name='fit_cont_v'))
    wx.StaticText(GUI_Global.popup_continuum, -1, 
        'Fitting continuum ON', pos=(70,20))
    # text field
    wx.StaticText(GUI_Global.popup_continuum, -1, 
        'Offset constant: ', pos=(30,64))
    GUI_Global.popup_continuum.cont_const_t = \
        wx.TextCtrl(GUI_Global.popup_continuum, -1, "", 
         size=(100, -1), pos=(160,60), style=wx.TE_RIGHT,
         validator=FloatValidator(base=GUI_Global.popup_continuum,
                                  name='cont_const_v'))

    wx.StaticText(GUI_Global.popup_continuum, -1, 
        'Poly. degree: ', pos=(30,104))
    GUI_Global.popup_continuum.cont_n_t = \
        wx.TextCtrl(GUI_Global.popup_continuum, -1, "", 
         size=(100, -1), pos=(160,100), style=wx.TE_RIGHT,
         validator=IntValidator(base=GUI_Global.popup_continuum,
                                name='cont_n_v',
                                value_check=ge_zero))


#===============================================================================
#
#  POPUP for the PATH Parameters
#

def OnButton_path_done(e):
    EndBusyCursor()
    # check all parameters / files on consistency. if not 
    # dont allow exit 
    if GUI_Parameters.check_path_and_files() == False:
        return
    GUI_Global.popup_path.Show(False) 
    if len(GUI_Global.popup_path.col_dflux_t.GetValue()) == 0:
        GUI_Global.popup_path.col_dflux_v = 'NULL'
        GUI_Global.popup_path.col_dflux_t.SetValue('NULL')
    if len(GUI_Global.popup_path.col_mask_t.GetValue()) == 0:
        GUI_Global.popup_path.col_mask_v = 'NULL'
        GUI_Global.popup_path.col_mask_t.SetValue('NULL')


# def OnText_basedir_Change(e):
#     print 'basedir change is prohibited'
#     # GUI_Global.param_changes = True
#     #GUI_Global.basedir_v = GUI_Global.popup_path.basedir_t.GetValue()
#     #GUI_Global.popup_path.basedir_t.SetValue(GUI_Global.basedir_v)
        
def OnText_workdir_Change(e):
    GUI_Global.param_changes = True
    GUI_Global.workdir_v = GUI_Global.popup_path.workdir_t.GetValue()
    # GUI_Parameters.absolute_path()
        
def OnText_path_Change(e):
    target = e.GetEventObject()
    #if GUI_Global.popup_path.basedir_t == target:
    #    GUI_Global.basedir_v = target.GetValue()
    #    GUI_Global.param_changes = True
    if GUI_Global.popup_path.workdir_t == target:
        GUI_Global.workdir_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.file == target:
        GUI_Global.filename_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.w_include == target:
        GUI_Global.wrange_include_v = target.GetValue()
        GUI_Global.param_changes = True
        #if os.path.exists(GUI_Global.wrange_include_v) == True:
        #    GUI_Data.import_mask_files(
        #            GUI_Global.wrange_include_v, False, False)
    if GUI_Global.popup_path.w_exclude == target:
        GUI_Global.wrange_exclude_v = target.GetValue()
        GUI_Global.param_changes = True
        #if os.path.exists(GUI_Global.wrange_exclude_v) == True:
        #    GUI_Data.import_mask_files(
        #            GUI_Global.wrange_exclude_v, False, False)
    if GUI_Global.popup_path.p_exclude == target:
        GUI_Global.prange_exclude_v = target.GetValue()
        GUI_Global.param_changes = True
        #if os.path.exists(GUI_Global.prange_exclude_v) == True:
        #    GUI_Data.import_mask_files(
        #            GUI_Global.prange_exclude_v, False, False)
    if GUI_Global.popup_path.output_name_t == target:
        GUI_Global.output_name_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.col_lam_t == target:
        GUI_Global.col_lam_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.col_flux_t == target:
        GUI_Global.col_flux_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.col_dflux_t == target:
        GUI_Global.col_dflux_v = target.GetValue()
        GUI_Global.param_changes = True
    if GUI_Global.popup_path.col_mask_t == target:
        GUI_Global.col_mask_v = target.GetValue()
        GUI_Global.param_changes = True

def refresh_popup_path():
    text = 'Current configuration file:   %s' % GUI_Global.Current_Setup_File
    GUI_Global.popup_path.config.SetLabel(text)
    # GUI_Global.popup_path.basedir_t.SetValue(GUI_Global.basedir_v)
    #GUI_Global.popup_path.basedir_t.SetLabel(GUI_Global.basedir_v)
    GUI_Global.popup_path.workdir_t.SetValue(GUI_Global.workdir_v)
    # GUI_Parameters.absolute_path()
    if (hasattr(GUI_Global,'filename_v') == True):
        GUI_Global.popup_path.file.SetValue(GUI_Global.filename_v) 
    else:
        GUI_Global.filename_v = '** not yet given **'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.file.SetValue(GUI_Global.filename_v)

    if (hasattr(GUI_Global,'wrange_include_v') == True):
        GUI_Global.popup_path.w_include.SetValue(GUI_Global.wrange_include_v) 
    else:
        GUI_Global.wrange_include_v = 'none'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.w_include.SetValue(GUI_Global.wrange_include_v)

    if (hasattr(GUI_Global,'wrange_exclude_v') == True):
        GUI_Global.popup_path.w_exclude.SetValue(GUI_Global.wrange_exclude_v) 
    else:
        GUI_Global.wrange_exclude_v = 'none'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.w_exclude.SetValue(GUI_Global.wrange_exclude_v)

    if (hasattr(GUI_Global,'prange_exclude_v') == True):
        GUI_Global.popup_path.p_exclude.SetValue(GUI_Global.prange_exclude_v) 
    else:
        GUI_Global.prange_exclude_v = 'none'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.p_exclude.SetValue(GUI_Global.prange_exclude_v)
        
    if (hasattr(GUI_Global,'output_name_v') == True):
        GUI_Global.popup_path.output_name_t.SetValue(GUI_Global.output_name_v) 
    else:
        GUI_Global.output_name_v = 'default_result'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.output_name_t.SetValue(GUI_Global.output_name_v)
        
    if (hasattr(GUI_Global,'col_lam_v') == True):
        GUI_Global.popup_path.col_lam_t.SetValue(GUI_Global.col_lam_v) 
    else:
        GUI_Global.col_lam_v = '** not yet given **'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.col_lam_t.SetValue(GUI_Global.col_lam_v)

    if (hasattr(GUI_Global,'col_flux_v') == True):
        GUI_Global.popup_path.col_flux_t.SetValue(GUI_Global.col_flux_v) 
    else:
        GUI_Global.col_flux_v = '** not yet given **'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.col_flux_t.SetValue(GUI_Global.col_flux_v)
        
    if (hasattr(GUI_Global,'col_dflux_v') == True):
        GUI_Global.popup_path.col_dflux_t.SetValue(GUI_Global.col_dflux_v) 
    else:
        GUI_Global.col_dflux_v = 'NULL'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.col_dflux_t.SetValue(GUI_Global.col_dflux_v)
        
    if (hasattr(GUI_Global,'col_mask_v') == True):
        GUI_Global.popup_path.col_mask_t.SetValue(GUI_Global.col_mask_v) 
    else:
        GUI_Global.col_mask_v = 'NULL'
        GUI_Global.param_changes = True
        GUI_Global.popup_path.col_mask_t.SetValue(GUI_Global.col_mask_v)
        
def OnImportMasks(e):
    #print "READ MASKS"
    if os.path.exists(GUI_Global.wrange_include_v) == True:
        GUI_Data.import_mask_files(
                GUI_Global.wrange_include_v, True, True)
    if os.path.exists(GUI_Global.wrange_exclude_v) == True:
        GUI_Data.import_mask_files(
                GUI_Global.wrange_exclude_v, False, True)
    if os.path.exists(GUI_Global.prange_exclude_v) == True:
        GUI_Data.import_mask_files(
                GUI_Global.prange_exclude_v, False, False)

def OnImportFits(e):
    #print GUI_Global.col_lam_v
    #print GUI_Global.col_flux_v
    if GUI_Global.col_lam_v == '** not yet given **':
        text = 'ERROR: no definition of FITS header/column information for the'\
            ' "Wavelength:" is given yet'
        #print text
        dlg = wx.MessageDialog( GUI_Global.frame, text,
	    GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy()   # finally destroy it when finished with OK
        return
    if GUI_Global.col_flux_v == '** not yet given **':
        text = 'ERROR: no definition of FITS header/column information for the'\
            ' "Flux:" is given yet'
        dlg = wx.MessageDialog( GUI_Global.frame, text,
	    GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy()   # finally destroy it when finished with OK
        #print text
        return
        
    if os.path.exists(GUI_Global.filename_v) == True:
        #print 'fits file exists'
        retval = GUI_Data.prepeare_fits_file(GUI_Global.filename_v)
        
        if retval == 0:

            prepearefile = \
                GUI_Data.output_dir_path() + '/' + \
                GUI_Global.output_name_v + '_gui.fits'
            GUI_Data.input_data(prepearefile,False,False,False)
            if GUI_Global.X[0] < 0.28 or GUI_Global.X[0] > 8600.:
                text = 'ERROR: First wavelength value is %g.\n'\
                 '       This is outside the valid range of 0.28-8600. microns\n'\
                 'Please check in check "Wavelength to micron conversion" in\n'\
                 '"Settings->Other Parameters" and then reload the data' %\
                 GUI_Global.X[0]
                dlg = wx.MessageDialog( GUI_Global.frame, text,
	            GUI_Global.VERSION, wx.OK)
                dlg.ShowModal() 
                dlg.Destroy() 
                return
            if GUI_Global.X[len(GUI_Global.X)-1] < 0.28 or \
                GUI_Global.X[len(GUI_Global.X)-1] > 8600.:
                text = 'ERROR: Last wavelength value is %g.\n'\
                 '       This is outside the valid range of 0.28-8600. microns\n'\
                 'Please check in check "Wavelength to micron conversion" in\n'\
                 '"Settings->Other Parameters" and then reload the data'%\
                 GUI_Global.X[len(GUI_Global.X)-1]
                dlg = wx.MessageDialog( GUI_Global.frame, text,
	            GUI_Global.VERSION, wx.OK)
                dlg.ShowModal() 
                dlg.Destroy() 
                return
            GUI_Data.read_input_file_pixel_exclusion_mask(prepearefile)
            GUI_Plot.draw_figure(False)
            if os.path.exists(GUI_Global.wrange_include_v) == True:
                GUI_Data.import_mask_files(
                    GUI_Global.wrange_include_v, True, True)
            if os.path.exists(GUI_Global.wrange_exclude_v) == True:
                GUI_Data.import_mask_files(
                    GUI_Global.wrange_exclude_v, False, True)
            if os.path.exists(GUI_Global.prange_exclude_v) == True:
                GUI_Data.import_mask_files(
                    GUI_Global.prange_exclude_v, False, False)
        else:
            text = 'The command prepguitable failed '\
              '(see molecfit documentation)\n'\
              'Use CONSOLE window to read molecfit/prepguitable result.\n\n' \
              'Also check if the FITS file structure (column names) are \n'\
              'given correctly !'
            dlg = wx.MessageDialog( GUI_Global.frame, text,
	        GUI_Global.VERSION, wx.OK)
            dlg.ShowModal() # Show it
            dlg.Destroy()   # finally destroy it when finished with OK
    else:
        text = 'ERROR: FITS file\n%s\n'\
            ' not found' % GUI_Global.filename_v
        dlg = wx.MessageDialog( GUI_Global.frame, text,
	    GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy()   # finally destroy it when finished with OK
        #print text
        return


def show_popup_path(parent):
    refresh_popup_path()
    GUI_Global.popup_path.Show(True)

def  OnButton_browse(e):
    if (GUI_Global.popup_path.button_browse_file == e.GetEventObject()):
        dialog = wx.FileDialog ( None, message = "Select a FITS file",
            style = wx.DD_DEFAULT_STYLE) 
        dialog.SetDirectory(GUI_Global.workdir_v)
        if dialog.ShowModal() == wx.ID_OK:
            selected = dialog.GetPath()
            GUI_Global.filename_v = selected
            GUI_Global.param_changes = True

    if (GUI_Global.popup_path.button_browse_work == e.GetEventObject()):
        dialog = wx.DirDialog ( None, message = "Select a directory",
            style = wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST) 
        #print GUI_Global.workdir_v
        dialog.SetPath(GUI_Global.workdir_v)
        if dialog.ShowModal() == wx.ID_OK:
            selected = dialog.GetPath()
            #print selected
            GUI_Global.workdir_v = selected
            # GUI_Parameters.absolute_path()
            GUI_Global.param_changes = True

    if (GUI_Global.popup_path.button_browse_include == e.GetEventObject()):
        dialog = wx.FileDialog ( None, message = "Select fit regions file",
            style = wx.DD_DEFAULT_STYLE) 
        dialog.SetDirectory(GUI_Global.workdir_v)
        if dialog.ShowModal() == wx.ID_OK:
            selected = dialog.GetPath()
            GUI_Global.wrange_include_v = selected
            GUI_Global.param_changes = True

    if (GUI_Global.popup_path.button_browse_exclude_w == e.GetEventObject()):
        dialog = wx.FileDialog ( None, message = "Select wavelength exclusion"+
            "mask file",
            style = wx.DD_DEFAULT_STYLE) 
        #print "working dir -->%s<--" % GUI_Global.workdir_v
        dialog.SetDirectory(GUI_Global.workdir_v)
        if dialog.ShowModal() == wx.ID_OK:
            selected = dialog.GetPath()
            GUI_Global.wrange_exclude_v = selected
            GUI_Global.param_changes = True

    if (GUI_Global.popup_path.button_browse_exclude_p == e.GetEventObject()):
        dialog = wx.FileDialog ( None, message = "Select pixel range exclusion"+
            "mask file",
            style = wx.DD_DEFAULT_STYLE) 
        dialog.SetDirectory(GUI_Global.workdir_v)
        if dialog.ShowModal() == wx.ID_OK:
            selected = dialog.GetPath()
            GUI_Global.prange_exclude_v = selected
            GUI_Global.param_changes = True

    refresh_popup_path()

def create_popup_path(parent):
    GUI_Global.popup_path = wx.Dialog(parent, wx.NewId(), 
        'FILE LOCATIONS and STRUCTURE', 
        size=(950, 720), style = wx.RESIZE_BORDER | wx.CLOSE_BOX)

    GUI_Global.popup_path.config = wx.StaticText(GUI_Global.popup_path, -1, 
        'Current configuration file: -- None --', pos=(30,24))

    # wx.StaticText(GUI_Global.popup_path, -1, 
    #     'MOLECFIT base installation path (write protected)', pos=(30,64))
    # GUI_Global.popup_path.basedir_t = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
    #      size=(500, -1), pos=(380,60), style=wx.TE_LEFT | wx.TE_READONLY)

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Base WORKING directory: ', pos=(30,104))
    GUI_Global.popup_path.workdir_t = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
         size=(600, -1), pos=(250,100), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_workdir_Change, 
        GUI_Global.popup_path.workdir_t ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Parameters to override settings in the configuration file (filenames'\
        ' without absolute path are interpreted relative to WORKING dir) ', 
        pos=(30,144))
    wx.StaticLine(GUI_Global.popup_path,-1, (10,176), (930,0) )
    wx.StaticText(GUI_Global.popup_path, -1, 
        'Spectrum (FITS):', pos=(30,184))
    GUI_Global.popup_path.file = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
         size=(600, -1), pos=(250,180), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.file ) 

    wx.StaticLine(GUI_Global.popup_path,-1, (10,176+260), (930,0) )
    wx.StaticText(GUI_Global.popup_path, -1, 
        'Regions to include for fitting (in micron):', pos=(30,184+260))
    GUI_Global.popup_path.w_include = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
         size=(600, -1), pos=(250,180+260), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.w_include ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Regions to exclude from the\nfitting (in micron):', pos=(30,224+260))
    GUI_Global.popup_path.w_exclude = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
         size=(600, -1), pos=(250,220+260), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.w_exclude ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Regions to exclude from the\nfitting (in pixels, 1st px = 1):', pos=(30,264+260))
    GUI_Global.popup_path.p_exclude = wx.TextCtrl(GUI_Global.popup_path, -1, "", 
         size=(600, -1), pos=(250,260+260), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.p_exclude ) 

    wx.StaticLine(GUI_Global.popup_path,-1, (10,296+300), (930,0) )
    wx.StaticText(GUI_Global.popup_path, -1, 
        'Root name for the output files:', pos=(30,304+300))
    GUI_Global.popup_path.output_name_t = wx.TextCtrl(GUI_Global.popup_path, 
        -1, "", size=(200, -1), pos=(250,300+300), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.output_name_t ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Parameters to set structure of input FITS file(s).\n'+
        'Names of the file columns (table) or extensions (image) containing:', 
        pos=(30,179+80))
    wx.StaticText(GUI_Global.popup_path, -1, 
        '\"Wavelength\" and \"Flux\" are required. \"Flux_Error\" and '+
        '\"Mask\" can be \"NULL\" (or empty)', 
        pos=(30,484-90))
    wx.StaticText(GUI_Global.popup_path, -1, 
        'Wavelength:', pos=(30,404-90))
    GUI_Global.popup_path.col_lam_t = wx.TextCtrl(GUI_Global.popup_path, -1,"", 
         size=(150, -1), pos=(150,400-90), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.col_lam_t ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Flux:', pos=(380,404-90))
    GUI_Global.popup_path.col_flux_t = wx.TextCtrl(GUI_Global.popup_path,-1,"", 
         size=(150, -1), pos=(440,400-90), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.col_flux_t ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Flux_Error:', pos=(30,444-90))
    GUI_Global.popup_path.col_dflux_t = wx.TextCtrl(GUI_Global.popup_path,-1,"", 
         size=(150, -1), pos=(150,440-90), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.col_dflux_t ) 

    wx.StaticText(GUI_Global.popup_path, -1, 
        'Mask:', pos=(380,444-90))
    GUI_Global.popup_path.col_mask_t = wx.TextCtrl(GUI_Global.popup_path, 
        -1, "", size=(150, -1), pos=(440,440-90), style=wx.TE_LEFT)
    GUI_Global.popup_path.Bind(wx.EVT_TEXT, OnText_path_Change, 
        GUI_Global.popup_path.col_mask_t ) 

    GUI_Global.popup_path.button_browse_work = wx.Button(
        GUI_Global.popup_path, id=-1, label="browse", 
        pos=(850,99))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_browse, 
        GUI_Global.popup_path.button_browse_work ) 

    GUI_Global.popup_path.button_browse_file = wx.Button(
        GUI_Global.popup_path, id=-1, label="browse", 
        pos=(850,179))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_browse, 
        GUI_Global.popup_path.button_browse_file ) 

    GUI_Global.popup_path.button_browse_include = wx.Button(
        GUI_Global.popup_path, id=-1, label="browse", 
        pos=(850,179+260))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_browse, 
        GUI_Global.popup_path.button_browse_include ) 

    GUI_Global.popup_path.button_browse_exclude_w = wx.Button(
        GUI_Global.popup_path, id=-1, label="browse", 
        pos=(850,219+260))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_browse, 
        GUI_Global.popup_path.button_browse_exclude_w ) 

    GUI_Global.popup_path.button_browse_exclude_p = wx.Button(
        GUI_Global.popup_path, id=-1, label="browse", 
        pos=(850,259+260))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_browse, 
        GUI_Global.popup_path.button_browse_exclude_p ) 

    GUI_Global.popup_path.button_read_masks = wx.Button(
        GUI_Global.popup_path, id=-1, 
        label="read / import mask from files to GUI", 
        pos=(680,300+260))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnImportMasks, 
        GUI_Global.popup_path.button_read_masks ) 

    GUI_Global.popup_path.button_read_fits = wx.Button(
        GUI_Global.popup_path, id=-1, 
        label="read / import FITS file to GUI", 
        pos=(720,219))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnImportFits, 
        GUI_Global.popup_path.button_read_fits ) 

    # DONE BUTTON
    GUI_Global.popup_path.button_done = wx.Button(
        GUI_Global.popup_path, id=-1, label="DONE", 
        pos=(750,520+120))
    GUI_Global.popup_path.Bind(wx.EVT_BUTTON, OnButton_path_done, 
        GUI_Global.popup_path.button_done ) 
    GUI_Global.popup_path.button_done.SetToolTip(wx.ToolTip(
        "click to if you are ready with this step"))
    GUI_Global.popup_path.button_done.Show(True)

#===============================================================================
#
#  POPUP for the WAVELENGTH CORRECTION
#
def refresh_popup_wavelength():
    if GUI_Global.popup_wavelength.fit_wlc_v != 0:
        GUI_Global.popup_wavelength.check_fit.SetValue(True)
    else:
        GUI_Global.popup_wavelength.check_fit.SetValue(False)
    text = repr(GUI_Global.popup_wavelength.wlc_const_v)
    GUI_Global.popup_wavelength.wlc_const_t.SetValue(text)
    text = "%i" % (GUI_Global.popup_wavelength.wlc_n_v)
    GUI_Global.popup_wavelength.wlc_n_t.SetValue(text)
    reset_textcontrol_colors(GUI_Global.popup_wavelength)
   
def show_popup_wavelength(parent):
    if GUI_Global.param_available == False:
        no_parameter_message()
        return
    refresh_popup_wavelength()
    GUI_Global.popup_wavelength.Show(True)

def create_popup_wavelength(parent):
    GUI_Global.popup_wavelength = wx.Dialog(parent, wx.NewId(), 
        'WAVELENGTH CORRECTION', size=(270, 220))

    self = GUI_Global.popup_wavelength
    
    sizer = wx.BoxSizer( wx.VERTICAL )
    sizer.AddSpacer(140)

    buttons = wx.StdDialogButtonSizer()

    self.buttonsOK = wx.Button(self, wx.ID_OK)
    buttons.AddButton(self.buttonsOK)
    self.buttonsCancel = wx.Button(self, wx.ID_CANCEL)
    buttons.AddButton(self.buttonsCancel)
    buttons.Realize();
    
    # generate checkbox and set default / initial value
    GUI_Global.popup_wavelength.check_fit = wx.CheckBox(
        GUI_Global.popup_wavelength, -1, "", pos=(30,20), style=wx.ALIGN_RIGHT,
        validator=BoolValidator(base=self, name='fit_wlc_v'))
    wx.StaticText(GUI_Global.popup_wavelength, -1, 
        'Correcting wavelength ON', pos=(70,20))
    # text field
    wx.StaticText(GUI_Global.popup_wavelength, -1, 
        'Offset constant: ', pos=(30,64))
    GUI_Global.popup_wavelength.wlc_const_t = wx.TextCtrl(
        GUI_Global.popup_wavelength, -1, "", 
         size=(100, -1), pos=(160,60), style=wx.TE_RIGHT,
         validator=FloatValidator(base=self, name='wlc_const_v'))

    wx.StaticText(GUI_Global.popup_wavelength, -1, 
        'Poly. degree: ', pos=(30,104))
    GUI_Global.popup_wavelength.wlc_n_t = wx.TextCtrl(
        GUI_Global.popup_wavelength, -1, "", 
         size=(100, -1), pos=(160,100), style=wx.TE_RIGHT,
         validator=IntValidator(base=self, name='wlc_n_v',
                                value_check=ge_zero))

    sizer.Add( buttons, 1, wx.EXPAND | wx.ALIGN_BOTTOM, 5 )

    self.SetSizer(sizer)
    
    self.Layout()

#===============================================================================
#
#  POPUP for the FITTING KERNEL
#
def refresh_popup_kernel():
    dlg = GUI_Global.popup_kernel
    dlg.check_box.SetValue(dlg.fit_res_box_v)
    dlg.check_kernmode.SetValue(dlg.kernmode_v)
    dlg.check_gauss.SetValue(dlg.fit_res_gauss_v)
    dlg.check_lorentz.SetValue(dlg.fit_res_lorentz_v)
    dlg.check_varkern.SetValue(dlg.varkern_v)

    # text field
    dlg.relres_box_t.SetValue(repr(dlg.relres_box_v))
    dlg.res_gauss_t.SetValue(repr(dlg.res_gauss_v))
    dlg.res_lorentz_t.SetValue(repr(dlg.res_lorentz_v))
    dlg.kernfac_t.SetValue("%i" % (dlg.kernfac_v))
    # to be backward compatible to older setup files this if is required
    if hasattr(dlg,'kernel_file_v') == False:
        dlg.kernel_file_v = 'none'
    text = "%s" % (dlg.kernel_file_v)
    dlg.kernel_file_t.SetPath(text)
    reset_textcontrol_colors(dlg)

def show_popup_kernel(parent):
    if GUI_Global.param_available == False:
        no_parameter_message()
        return
    refresh_popup_kernel()
    GUI_Global.popup_kernel.Show(True)

def create_popup_kernel(parent):
    dlg = FitKernelDialog(parent)
    dlg.check_box.SetValidator(BoolValidator(base=dlg, name='fit_res_box_v'))
    dlg.relres_box_t.SetValidator(FloatValidator(base=dlg,
                                 name='relres_box_v',
                                 value_check=partial(in_range, low=0, high=2)))
    dlg.check_gauss.SetValidator(BoolValidator(base=dlg,
                                               name='fit_res_gauss_v'))
    dlg.res_gauss_t.SetValidator(FloatValidator(base=dlg, name='res_gauss_v',
                                                value_check=ge_zero))
    dlg.check_lorentz.SetValidator(BoolValidator(base=dlg,
                                                 name='fit_res_lorentz_v'))
    dlg.res_lorentz_t.SetValidator(FloatValidator(base=dlg,
                                                  name='res_lorentz_v',
                                                  value_check=ge_zero))
    dlg.check_kernmode.SetValidator(BoolValidator(base=dlg, name='kernmode_v'))
    dlg.kernfac_t.SetValidator(FloatValidator(base=dlg, name='kernfac_v',
                                              value_check=ge_zero))
    dlg.check_varkern.SetValidator(BoolValidator(base=dlg, name='varkern_v'))
    GUI_Global.popup_kernel = dlg


#===============================================================================
#
#  POPUP for the TELESCOPE BACKGROUND
#
def refresh_popup_background():
    if GUI_Global.popup_back.fit_back_v != 0:
        GUI_Global.popup_back.check_fit_back.SetValue(True)
    else:
        GUI_Global.popup_back.check_fit_back.SetValue(False)
    text = repr(GUI_Global.popup_back.telback_v)
    GUI_Global.popup_back.telback_t.SetValue(text)
    reset_textcontrol_colors(GUI_Global.popup_back)

def show_popup_back(parent):
    if GUI_Global.param_available == False:
        no_parameter_message()
        return
    refresh_popup_background()
    GUI_Global.popup_back.Show(True)

def create_popup_background(parent):
    GUI_Global.popup_back = wx.Dialog(parent, wx.NewId(), 'BACKGROUND', 
                                      size=(400, 155))

    self = GUI_Global.popup_back
    # generate checkbox and set default / initial value
    GUI_Global.popup_back.check_fit_back = wx.CheckBox(GUI_Global.popup_back, 
        -1, "", style=wx.ALIGN_RIGHT,
        validator=BoolValidator(base=self,
                                name='fit_back_v'))
    fnameLabel0 = wx.StaticText(GUI_Global.popup_back, -1, 
        'Fitting background ON')
    # text field
    fnameLabel = wx.StaticText(GUI_Global.popup_back, -1, 
        'Telescope emissivity: ')
    GUI_Global.popup_back.telback_t = wx.TextCtrl(GUI_Global.popup_back, -1, "", 
         size=(200, -1), style=wx.TE_RIGHT,
         validator=FloatValidator(base=self, name='telback_v'))
    csizer = wx.FlexGridSizer(cols = 2 , hgap = 4, vgap = 10)
    csizer.Add(GUI_Global.popup_back.check_fit_back, 0,
               wx.ALIGN_RIGHT|wx.TOP|wx.LEFT, 5)
    csizer.Add(fnameLabel0, 0, wx.RIGHT|wx.TOP|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
    csizer.Add(fnameLabel, 0, wx.ALIGN_RIGHT|wx.TOP|wx.LEFT, 5)
    csizer.Add(GUI_Global.popup_back.telback_t, 0, wx.RIGHT|wx.TOP|wx.LEFT, 5)

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(csizer, wx.EXPAND, 5)

    buttons = wx.StdDialogButtonSizer()
    self.buttonsOK = wx.Button(self, wx.ID_OK)
    buttons.AddButton(self.buttonsOK)
    self.buttonsCancel = wx.Button(self, wx.ID_CANCEL)
    buttons.AddButton(self.buttonsCancel)
    buttons.Realize();
    sizer.Add(buttons, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5)
    self.SetSizer(sizer)
    self.Layout()

#===============================================================================
#
#  POPUP for the OTHER parameters
#
def show_popup_other(parent):
    refresh_popup_other()
    GUI_Global.popup_other.Show(True)

def refresh_popup_other():
    if (hasattr(GUI_Global,'trans_v') != True):
        GUI_Global.trans_v = 1
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'vac_air_v') != True):
        GUI_Global.vac_air_v = 'vac'
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'wlgtomicron_v') != True):
        GUI_Global.wlgtomicron_v = 1e-3
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'ftol_v') != True):
        GUI_Global.ftol_v = 1e-2
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'xtol_v') != True):
        GUI_Global.xtol_v = 1e-2
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'default_error_v') != True):
        GUI_Global.default_error_v = 0.01
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'flux_unit_v') != True):
        GUI_Global.flux_unit_v = 0
        GUI_Global.param_changes = True
    if (hasattr(GUI_Global,'pixsc_v') != True):
        GUI_Global.pixsc_v = 0.16
        GUI_Global.param_changes = True

    if GUI_Global.trans_v:
        GUI_Global.popup_other.check_trans_t.SetValue(True)
    else:
        GUI_Global.popup_other.check_emiss_t.SetValue(True)

    if GUI_Global.vac_air_v == 'vac':
        GUI_Global.popup_other.check_vac_t.SetValue(True)
    else:
        GUI_Global.popup_other.check_air_t.SetValue(True)

    if GUI_Global.flux_unit_v == 0:
        GUI_Global.popup_other.check_flux_0_t.SetValue(True)
    elif GUI_Global.flux_unit_v == 1:
        GUI_Global.popup_other.check_flux_1_t.SetValue(True)
    elif GUI_Global.flux_unit_v == 2:
        GUI_Global.popup_other.check_flux_2_t.SetValue(True)
    elif GUI_Global.flux_unit_v == 3:
        GUI_Global.popup_other.check_flux_3_t.SetValue(True)
    else:
        raise ValueError("invalid flux unit %r" % GUI_Global.flux_unit_v)

    text = "%r" % (GUI_Global.wlgtomicron_v)
    GUI_Global.popup_other.wlgtomicron_t.SetValue(text)
    GUI_Global.popup_other.ftol_t.SetValue(repr(GUI_Global.ftol_v))
    GUI_Global.popup_other.xtol_t.SetValue(repr(GUI_Global.xtol_v))
    text = "%r" % (GUI_Global.default_error_v)
    GUI_Global.popup_other.default_error_t.SetValue(text)
    GUI_Global.popup_other.pixsc_t.SetValue(repr(GUI_Global.pixsc_v))
    reset_textcontrol_colors(GUI_Global.popup_other)

def create_popup_other(parent):
    GUI_Global.popup_other = wx.Dialog(parent, wx.NewId(), 
        'OTHER PARAMETERS', size=(460, 570+80), 
        style = wx.RESIZE_BORDER | wx.CLOSE_BOX)

    self = GUI_Global.popup_other
    sizer = wx.BoxSizer( wx.VERTICAL )
    
    # TODO layout full dialog with sizer instead of using this spacer
    sizer.AddSpacer(570)

    buttons = wx.StdDialogButtonSizer()
    self.button_done = wx.Button(self, wx.ID_OK)
    buttons.AddButton(self.button_done)
    self.buttonsCancel = wx.Button(self, wx.ID_CANCEL)
    buttons.AddButton(self.buttonsCancel)
    buttons.Realize();
    sizer.Add(buttons, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5)
    self.SetSizer(sizer)
    self.Layout()

    # trans parameter
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Type of input spectrum (transmission = default, emission): ', \
        pos=(30,20))
    self.check_trans_t = wx.RadioButton(self,
            -1, "transmission", pos=(40,50), style=wx.ALIGN_LEFT|wx.RB_GROUP,
            validator=ListValidator(base=GUI_Global, name='trans_v',
                                    value=1))
    self.check_emiss_t = wx.RadioButton(self,
            -1, "emission", pos=(160,50), style=wx.ALIGN_LEFT,
            validator=ListValidator(base=GUI_Global, name='trans_v',
                                    value=0))

    # toggling vacuum / air
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Type of wavelength calibration (scale in air or in vacuum): ', \
        pos=(30,90))
    self.check_vac_t = wx.RadioButton(self,
            -1, "vacuum", pos=(40,120), style=wx.ALIGN_LEFT|wx.RB_GROUP,
            validator=ListValidator(base=GUI_Global, name='vac_air_v',
                                    value='vac'))

    self.check_air_t = wx.RadioButton(self,
            -1, "air", pos=(160,120), style=wx.ALIGN_LEFT,
            validator=ListValidator(base=GUI_Global, name='vac_air_v',
                                    value='air'))

    # wavelength to microns conversion
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Wavelength to micron conversion factor: ', pos=(30,153))
    GUI_Global.popup_other.wlgtomicron_t = wx.TextCtrl(GUI_Global.popup_other, 
            -1, "", size=(100, -1), style=wx.TE_RIGHT, pos=(310,150),
            validator=FloatValidator(base=GUI_Global,
                                     name='wlgtomicron_v',
                                     value_check=g_zero))

    # toggling flux units
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Flux units: ', pos=(30,190))
    self.check_flux_0_t = wx.RadioButton(self,
            -1, "phot/(s*m^2*microns*arcsec^2) [no conversion]", pos=(40,220), 
            style=wx.ALIGN_LEFT|wx.RB_GROUP,
            validator=ListValidator(base=GUI_Global, name='flux_unit_v',
                                    value=0))

    self.check_flux_1_t = wx.RadioButton(self,
            -1, "W/(m^2*microns*arcsec^2)", pos=(40,250), style=wx.ALIGN_LEFT,
            validator=ListValidator(base=GUI_Global, name='flux_unit_v',
                                    value=1))

    self.check_flux_2_t = wx.RadioButton(self,
            -1, "erg/(s*cm^2*A*arcsec^2)", pos=(40,280), style=wx.ALIGN_LEFT,
            validator=ListValidator(base=GUI_Global, name='flux_unit_v',
                                    value=2))

    self.check_flux_3_t = wx.RadioButton(self,
            -1, "mJy/arcsec^2", pos=(40,280+30), style=wx.ALIGN_LEFT,
            validator=ListValidator(base=GUI_Global, name='flux_unit_v',
                                    value=3))
    if (hasattr(GUI_Global,'flux_unit_v') != True):
        GUI_Global.flux_unit_v = 0
    GUI_Global.popup_other.check_flux_0_t.SetValue(True)

    # FIT AND CALCULATION ERROR CRITERIA
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Fit and error criteria for calculation: ', pos=(30,320+30))
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Chi2 convergence criteria (ftol): ', pos=(30,353+30))
    GUI_Global.popup_other.ftol_t = wx.TextCtrl(GUI_Global.popup_other, 
            -1, "", size=(100, -1), style=wx.TE_RIGHT, pos=(290+50,350+30),
            validator=FloatValidator(base=GUI_Global,
                                     name='ftol_v',
                                     value_check=b_zero_one))
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Parameter convergence criteria (xtol): ', pos=(30,383+30))
    GUI_Global.popup_other.xtol_t = wx.TextCtrl(GUI_Global.popup_other, 
            -1, "", size=(100, -1), style=wx.TE_RIGHT, pos=(290+50,380+30),
            validator=FloatValidator(base=GUI_Global,
                                     name='xtol_v',
                                     value_check=b_zero_one))
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Data relative error estimate (default_error): ', pos=(30,413+30))
    GUI_Global.popup_other.default_error_t = wx.TextCtrl(GUI_Global.popup_other, 
            -1, "", size=(100, -1), style=wx.TE_RIGHT, pos=(290+50,410+30),
            validator=FloatValidator(base=GUI_Global,
                                     name='default_error_v',
                                     value_check=b_zero_one))
    wx.StaticText(GUI_Global.popup_other, -1, 
        'The data error estimate is only used if the input FITS file\n'+
        'does not give an error (or it is disabled by Flux_Error = NULL)', 
        pos=(50,440+30))
    # pixel size
    wx.StaticText(GUI_Global.popup_other, -1, 
        'Instrument pixel scale (arcsec/pix): ', pos=(30,503+30))
    GUI_Global.popup_other.pixsc_t = wx.TextCtrl(GUI_Global.popup_other, 
            -1, "", size=(100, -1), style=wx.TE_RIGHT, pos=(290,500+30),
            validator=FloatValidator(base=GUI_Global,
                                     name='pixsc_v',
                                     value_check=b_zero_one))


#===============================================================================
#
#  POPUP for the MOELCULAR data
#
def OnCheckButton(e):
    o = e.GetEventObject()
    if o == o.m.calc_t:
        # unchecking calc disables fit too
        if not o.IsChecked():
            o.m.fit_t.SetValue(o.IsChecked())
    if o == o.m.fit_t:
        # checking fit enables fit too
        if o.IsChecked():
            o.m.calc_t.SetValue(o.IsChecked())
    GUI_Global.param_changes = True

def OnPWVChange(e):
    """ give visual indication that negative means disabled """
    o = e.GetEventObject()
    try:
        o.chk.SetValue(float(o.GetValue()) >= 0.)
    except ValueError:
        # validation is done by validator later
        pass

def refresh_popup_molec():
    # refresh the data content (e.g. if changed via read param call
    # but already open from setup before)
    dflt = wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW)
    
    for m in GUI_Global.mol:
        m.calc_t.SetValue(m.calc)
        if hasattr(m, 'fit_t') and m.fit_t is not None:
            m.fit_t.SetValue(m.fit)
        m.ab_t.SetValue(repr(m.ab))
        m.ab_t.SetBackgroundColour(dflt)
    GUI_Global.pwv_par.refresh()
       
def show_popup_molec(parent):
    refresh_popup_molec()
    GUI_Global.popup_molec.Show(True)

def show_popup_molec_sub(parent):
    refresh_popup_molec()
    GUI_Global.popup_molec_sub.Show(True)

def create_popup_molec(parent):
    GUI_Global.popup_molec = wx.Dialog(parent, wx.NewId(), 
        'MOLECULAR ABUNDANCE', size=(580, 460))
    self = GUI_Global.popup_molec

    msizer = wx.BoxSizer(wx.VERTICAL)
    molec_sizer = wx.FlexGridSizer(0, 7, 0, 0)
    molec_sizer.SetFlexibleDirection(wx.HORIZONTAL)
    molec_sizer.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

    # title
    molec_sizer.AddSpacer(0)

    txt_calc = wx.StaticText(self, wx.ID_ANY, u"calc", wx.DefaultPosition,
                             wx.DefaultSize, wx.ALIGN_RIGHT)
    txt_calc.SetToolTip(wx.ToolTip(TextPool.tip_018))
    molec_sizer.Add(txt_calc, 0,
                    wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                    wx.ALIGN_CENTER_HORIZONTAL, 5)

    txt_fit = wx.StaticText(self, wx.ID_ANY, u"fit", wx.DefaultPosition,
                            wx.DefaultSize, 0)
    txt_fit.SetToolTip(wx.ToolTip(TextPool.tip_019))
    molec_sizer.Add(txt_fit, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                    wx.ALIGN_CENTER_HORIZONTAL, 5)

    txt_ab = wx.StaticText(self, wx.ID_ANY, u"relative abundance [1 = 100%]",
                           wx.DefaultPosition, wx.DefaultSize, 0 )
    txt_ab.SetToolTip(wx.ToolTip(TextPool.tip_020))
    molec_sizer.Add(txt_ab, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                    wx.ALIGN_CENTER_VERTICAL, 5)

    molec_sizer.AddSpacer(0)
    molec_sizer.AddSpacer(0)
    molec_sizer.AddSpacer(0)
    # title end

    for m in GUI_Global.mol[:GUI_Global.NR_OF_MAIN_MOLS]:
        txt = wx.StaticText(self, wx.ID_ANY, m.name, wx.DefaultPosition,
                                wx.DefaultSize, 0)
        molec_sizer.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                        wx.ALIGN_CENTER_HORIZONTAL, 5)
        
        m.calc_t = wx.CheckBox(self, wx.ID_ANY, wx.EmptyString,
                               wx.DefaultPosition, wx.DefaultSize,
                               wx.ALIGN_RIGHT,
                               validator=BoolValidator(base=m, name='calc'))
        m.calc_t.m = m
        self.Bind(wx.EVT_CHECKBOX, OnCheckButton, m.calc_t)
        m.calc_t.SetValue(m.calc)
        molec_sizer.Add(m.calc_t, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL, 5)
        
        m.fit_t = wx.CheckBox(self, wx.ID_ANY, wx.EmptyString,
                              wx.DefaultPosition, wx.DefaultSize,
                              wx.ALIGN_RIGHT,
                              validator=BoolValidator(base=m, name='fit'))
        m.fit_t.m = m
        self.Bind(wx.EVT_CHECKBOX, OnCheckButton, m.fit_t)
        m.fit_t.SetValue(m.fit)
        molec_sizer.Add(m.fit_t, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL, 5)
        
        m.ab_t = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString,
                             wx.DefaultPosition, wx.DefaultSize,
                             wx.TE_RIGHT,
                             validator=FloatValidator(base=m, name='ab',
                                  value_check=partial(in_range,
                                                      low=1e-5, high=100.)))
        m.ab_t.SetValue(repr(m.ab))
        molec_sizer.Add(m.ab_t, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL, 5)

        if m.name != "H2O":
            
            molec_sizer.AddSpacer(0)
            molec_sizer.AddSpacer(0)
            molec_sizer.AddSpacer(0)
    
        else:
            pwv_t = wx.CheckBox(self, wx.ID_ANY, wx.EmptyString,
                   wx.DefaultPosition, wx.DefaultSize,
                   wx.ALIGN_RIGHT)
            # disabled, pwv controlled by value of txtcontrol
            pwv_t.Enable(False)
            pwv_t.SetToolTip(wx.ToolTip(TextPool.tip_022))
            GUI_Global.pwv_par.chkctrl = pwv_t
            molec_sizer.Add(pwv_t, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                            wx.ALIGN_CENTER_HORIZONTAL, 5)
            
            txt_pwv = wx.StaticText(self, wx.ID_ANY, u"PWV", wx.DefaultPosition,
                                    wx.DefaultSize, 0)
            txt_pwv.SetToolTip(wx.ToolTip(TextPool.tip_022))
            molec_sizer.Add(txt_pwv, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
            
            pwv_v_t = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString,
                         wx.DefaultPosition, wx.DefaultSize,
                         wx.TE_RIGHT,
                         validator=FloatValidator(base=GUI_Global.pwv_par,
                                                  name='pwv'))
            pwv_v_t.chk = pwv_t
            self.Bind(wx.EVT_TEXT, OnPWVChange, pwv_v_t)
            GUI_Global.pwv_par.txtctrl = pwv_v_t
            molec_sizer.Add(pwv_v_t, 0, wx.ALL, 5)

    msizer.Add(molec_sizer, 1, wx.EXPAND, 5)
                
    self.sub_button = wx.Button(self, wx.ID_ANY, u"Additional trace gases",
                               wx.DefaultPosition, wx.DefaultSize, 0)
    self.Bind(wx.EVT_BUTTON, show_popup_molec_sub, self.sub_button)

    msizer.Add(self.sub_button, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5)
    
    m_sdbSizer = wx.StdDialogButtonSizer()
    m_sdbSizer.AddButton(wx.Button(self, wx.ID_OK))
    m_sdbSizer.AddButton(wx.Button(self, wx.ID_CANCEL))
    self.button_help = wx.Button(self, wx.ID_HELP)
    m_sdbSizer.AddButton(self.button_help)
    m_sdbSizer.Realize();

    self.Bind(wx.EVT_BUTTON, OnButton_help, self.button_help)
    
    msizer.Add(m_sdbSizer, 0, wx.EXPAND, 5)
    
    self.SetSizer(msizer)
    self.Layout()
    msizer.Fit(self)
    

#===============================================================================
#
#  SUB POPUP for the MOELCULAR data
#
def create_popup_molec_sub(parent):
    GUI_Global.popup_molec_sub = wx.Dialog(parent, wx.NewId(), 
        'MOLECULAR ABUNDANCE of trace gases', 
        size=(450, 750))
    self = GUI_Global.popup_molec_sub

    msizer = wx.BoxSizer(wx.VERTICAL)
    molec_sizer = wx.FlexGridSizer(0, 3, 0, 0)
    molec_sizer.SetFlexibleDirection(wx.HORIZONTAL)
    molec_sizer.SetNonFlexibleGrowMode(wx.FLEX_GROWMODE_SPECIFIED)

    # title
    molec_sizer.AddSpacer(0)

    txt_calc = wx.StaticText(self, wx.ID_ANY, u"calc", wx.DefaultPosition,
                             wx.DefaultSize, wx.ALIGN_RIGHT)
    txt_calc.SetToolTip(wx.ToolTip(TextPool.tip_018))
    molec_sizer.Add(txt_calc, 0,
                    wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                    wx.ALIGN_CENTER_HORIZONTAL, 5)

    txt_ab = wx.StaticText(self, wx.ID_ANY, u"relative abundance [1 = 100%]",
                           wx.DefaultPosition, wx.DefaultSize, 0 )
    txt_ab.SetToolTip(wx.ToolTip(TextPool.tip_020))
    molec_sizer.Add(txt_ab, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                    wx.ALIGN_CENTER_VERTICAL, 5)
    # title end

    for m in GUI_Global.mol[GUI_Global.NR_OF_MAIN_MOLS:]:
        txt = wx.StaticText(self, wx.ID_ANY, m.name, wx.DefaultPosition,
                                wx.DefaultSize, 0)
        molec_sizer.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|
                        wx.ALIGN_CENTER_HORIZONTAL, 5)
        
        m.calc_t = wx.CheckBox(self, wx.ID_ANY, wx.EmptyString,
                               wx.DefaultPosition, wx.DefaultSize,
                               wx.ALIGN_RIGHT,
                               validator=BoolValidator(base=m, name='calc'))
        m.calc_t.SetValue(m.calc)
        molec_sizer.Add(m.calc_t, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL, 5)
        
        m.ab_t = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString,
                             wx.DefaultPosition, wx.DefaultSize,
                             wx.TE_RIGHT,
                             validator=FloatValidator(base=m, name='ab',
                                  value_check=partial(in_range,
                                                      low=1e-5, high=100.)))
        m.ab_t.SetValue(repr(m.ab))
        molec_sizer.Add(m.ab_t, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL, 5)

    msizer.Add(molec_sizer, 1, wx.EXPAND, 5)
                
    m_sdbSizer = wx.StdDialogButtonSizer()
    self.button_done = wx.Button(self, wx.ID_OK)
    m_sdbSizer.AddButton(self.button_done)
    m_sdbSizer.AddButton(wx.Button(self, wx.ID_CANCEL))
    self.button_help = wx.Button(self, wx.ID_HELP)
    m_sdbSizer.AddButton(self.button_help)
    m_sdbSizer.Realize();

    self.Bind(wx.EVT_BUTTON, OnButton_help, self.button_help)
    
    msizer.Add(m_sdbSizer, 0, wx.EXPAND, 5)
    
    self.SetSizer(msizer)
    self.Layout()
    msizer.Fit(self)

#===============================================================================
#
#  POPUP for the FIT procedure
#
def OnButton_abort(e):
    EndBusyCursor()
    GUI_Global.child.kill()
    GUI_Global.user_killed = True

def create_popup_fit(parent):
    GUI_Global.popup_fit = wx.Dialog(parent, wx.NewId(), 'RADIATIVE TRANSFER', 
        size=(250, 150), style=wx.STAY_ON_TOP | wx.DEFAULT_DIALOG_STYLE)
    GUI_Global.popup_fit.SetBackgroundColour((255,255,20))
    text = wx.StaticText(GUI_Global.popup_fit, -1, 
        '      ...radiative transfer code\n                           in execution...\n            GUI blocked  ', 
        pos=(5,3))
    text.SetForegroundColour((0,0,0))
    text.Show(True)
    #---------------------------------------------------------------------------
    # ABORT BUTTON   
    GUI_Global.popup_fit.button_abort = wx.Button(
        GUI_Global.popup_fit, 
        id=-1, label="ABORT", pos=(80,70))
    GUI_Global.popup_fit.Bind(wx.EVT_BUTTON, OnButton_abort, 
        GUI_Global.popup_fit.button_abort ) 
    GUI_Global.popup_fit.button_abort.SetToolTip(wx.ToolTip(
        "click to if you are ready want to abort the call"))
    GUI_Global.popup_fit.button_abort.Show(True)

def select_file(parent, **kwargs):
    kwargs['style'] = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
    dlg = wx.FileDialog(parent, **kwargs)

    if dlg.ShowModal() == wx.ID_CANCEL:
        return None
    else:
        return dlg.GetPath()

def check_mask_names(force):
    falsereturn = False
    # check if mask files are required in name setup
    pwd = os.getcwdu()
    if hasattr(GUI_Global,'prange_exclude_v') == False:
        GUI_Global.prange_exclude_v = 'none'
    if (GUI_Global.prange_exclude_v[:4] == 'none'):
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == True) and 
                (GUI_Global.mask[ii].include == False)):
                found = 1
        if found == 1 :
            path = select_file(GUI_Global.frame,
                   message="Select filename for pixel exclusion mask data",
                   defaultFile=os.path.join(pwd, "pixel_exclusion.dat"))
            if path is None:
                dlg  = wx.MessageDialog(GUI_Global.frame,
                    "Filename for pixel exclusion mask not yet defined\n\n"\
                    "goto \n   Settings->Working Directories and Files\n"\
                    "and add filenames.", "INFO", wx.OK)
                dlg.ShowModal()
                dlg.Destroy()
                falsereturn = True
            else:
                GUI_Global.prange_exclude_v = path
                GUI_Global.param_changes = True

    if hasattr(GUI_Global,'wrange_exclude_v') == False:
        GUI_Global.wrange_exclude_v = 'none'
    if (GUI_Global.wrange_exclude_v[:4] == 'none'):
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == False) and 
                (GUI_Global.mask[ii].include == False)):
                found = 1
        if found == 1 :
            path = select_file(GUI_Global.frame,
               message="Select filename for wavelength exclusion mask data",
               defaultFile=os.path.join(pwd, "wavelength_exclusion.dat"))
            if path is None:
                dlg  = wx.MessageDialog(GUI_Global.frame,
                    "Filename for Wavelength exclusion mask not yet defined\n\n"\
                    "goto \n   Settings->Working Directories and Files\n"\
                    "and add filenames.", "INFO", wx.OK)
                dlg.ShowModal()
                dlg.Destroy()
                falsereturn = True
            else:
                GUI_Global.wrange_exclude_v = path
                GUI_Global.param_changes = True
            

    if hasattr(GUI_Global,'wrange_include_v') == False:
        GUI_Global.wrange_include_v = 'none'
    if (GUI_Global.wrange_include_v[:4] == 'none'):
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == False) and 
                (GUI_Global.mask[ii].include == True)):
                found = 1
        if found == 1 :
            path = select_file(GUI_Global.frame,
               message="Select filename for wavelength inclusion mask data",
               defaultFile=os.path.join(pwd, "wavelength_inclusion.dat"))
            if path is None:
                dlg  = wx.MessageDialog(GUI_Global.frame,
                    "Filename for Wavelength inclusion mask not yet defined\n\n"\
                    "goto \n   Settings->Working Directories and Files\n"\
                    "and add filenames.", "INFO", wx.OK)
                dlg.ShowModal()
                dlg.Destroy()
                falsereturn = True
            else:
                GUI_Global.wrange_include_v = path
                GUI_Global.param_changes = True
            
    if falsereturn == True:
        return False

    # if force == True:
    #     if GUI_Global.prange_exclude_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.prange_exclude_v) != True:
    #             if os.path.exists(GUI_Global.basedir_v + '/' + \
    #                 GUI_Global.prange_exclude_v) != True:
    #                 text = 'ERROR: pixel exclusion mask file \n%s\n does' \
    #                         ' not exist.\n     Will be created now!' % \
    #                         GUI_Global.prange_exclude_v
    #                 GUI_Parameters.my_message(text)  
    #                 GUI_Parameters.generate_and_write_mask_files()
  
    #     if GUI_Global.wrange_exclude_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.wrange_exclude_v) != True:
    #             if os.path.exists(GUI_Global.basedir_v + '/' + \
    #                 GUI_Global.wrange_exclude_v) != True:
    #                 text = ('ERROR: wavelength exclusion mask file \n%s\n does' \
    #                         ' not exist\n     Will be created now!' % \
    #                     (GUI_Global.wrange_exclude_v))
    #                 GUI_Parameters.my_message(text) 
    #                 GUI_Parameters.generate_and_write_mask_files()
  
    #     if GUI_Global.wrange_include_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.wrange_include_v) != True:
    #             if os.path.exists(GUI_Global.basedir_v + '/' + \
    #                 GUI_Global.wrange_include_v) != True:
    #                 text = ('ERROR: wavelength fit (inclusion) mask file \n%s\n ' \
    #                         'does not exist\n     Will be created now!' % \
    #                     (GUI_Global.wrange_include_v))
    #                 GUI_Parameters.my_message(text)   
    #                 GUI_Parameters.generate_and_write_mask_files()

    return True

def show_popup_fit(parent):
    if GUI_Global.data_available == False:
        no_data_message()
        return
    if GUI_Parameters.check_path_and_files() == False:
        return
    
    
    if GUI_Global.param_changes == True:
        #dlg = Three_Button_Message_Dialogue(GUI_Global.frame,
        #       "Parameters have been changed!", 
        #       "ignore changes", "save setup file", "cancel")
        #retval = dlg.ShowModal() 
        #dlg.Destroy()   
        #if (retval == 3) :
        #    
        #    return
        #if (retval == 2) :
        #    GUI_Menu.OnSave(parent)
        GUI_Parameters.my_message(
            'parameters have been changed\nuse "Save" / "Save As" '\
            'before fit.')
        return

    # normally not needed here anymore - is incuded into Save/SaveAs
    # but users do crazy combinations of actions
    # mode True forces the writing
    if check_mask_names(True) == False:
        return

  
    # if GUI_Global.filename_v[:1] != '/':
    #     GUI_Global.filename_v = GUI_Global.basedir_v + '/' + \
    #         GUI_Global.filename_v 
    text1 = GUI_Global.filename_v
    if os.path.exists(text1) != True:
        text = ('ERROR: data input FITS file \n%s\n does'
                ' not exist' % (text1))
        GUI_Parameters.my_message(text)    
        return
  
    if GUI_Parameters.generate_and_write_mask_files() == False : 
        GUI_Parameters.my_message("ERROR: mask files not ready \nFIT aborted") 
        return

    file_fit = (GUI_Data.output_dir_path() + '/' + \
            GUI_Global.output_name_v + '_fit.fits')

    if os.path.exists(file_fit) == True:
        if (os.path.getmtime(GUI_Global.Current_Setup_File) <
            os.path.getmtime(file_fit)) and not GUI_Global.mask_changes:
            # no fit needed - exists already
            dlg = wx.MessageDialog(GUI_Global.frame, TextPool.dialog_011,
                GUI_Global.VERSION, wx.YES_NO | wx.ICON_QUESTION | 
                wx.NO_DEFAULT)
            result = dlg.ShowModal()             
            dlg.Destroy()                        
            if result == wx.ID_NO:
                GUI_Data.input_data(file_fit,False,True,False)
                GUI_Global.fit_available = True
                GUI_Global.panel.cb_model.SetValue(True)
                GUI_Plot.draw_figure(False)
                return  

            
    GUI_Global.molecfit_output = str('%s/%s_molecfit_out.txt' % (
        GUI_Data.output_dir_path(),GUI_Global.output_name_v))
    # run the fit in a terminal windowq
    cmd = str("%s %s" % (GUI_Global.molecfit_exe, 
        GUI_Global.Current_Setup_File))
    # print '##### DEBUG ######'
    # print 'cmd = ' + cmd
    # print '##### DEBUG ######'
    # print 'ENV'
    # import json
    # print json.dumps(dict(**os.environ), indent=4, sort_keys=True)
    # print '##### DEBUG ######'
    #os.system("pwd")
    #print "FIT command %s with logfile %s" % (cmd,GUI_Global.molecfit_output)
    wx.BeginBusyCursor(wx.HOURGLASS_CURSOR)
    retval = run_subcommand_in_popup(cmd, GUI_Global.molecfit_output, True)
    #print "command %s retval %i " % (cmd,retval)
    EndBusyCursor()
    if GUI_Global.user_killed == True:
        GUI_Console.hide_console(None)
        return

    if retval != 0:
        dlg = wx.MessageDialog( GUI_Global.frame, 
            "ERROR in running the FIT\nmolecfit failed - consult the result log"
            , "ERROR", wx.OK)
        retval = dlg.ShowModal() 
        dlg.Destroy()   
    else:
        GUI_Global.mask_changes = False
        prepearefile = ( GUI_Data.output_dir_path() + '/' + \
            GUI_Global.output_name_v + '_fit' + '.fits')
        #print "RESULT: %s" % prepearefile
        if os.path.exists(prepearefile) == True:
            GUI_Data.input_data(prepearefile,False,True,False)
            GUI_Global.fit_available = True
            GUI_Global.panel.cb_model.SetValue(True)
            GUI_Plot.draw_figure(False)
        else:
            dlg = wx.MessageDialog( GUI_Global.frame, 
                "ERROR: Result file not found", "ERROR", wx.OK)
            retval = dlg.ShowModal() 
            dlg.Destroy()   
            #print "result file not found"
    
    ##GUI_Global.popup_fit.Show(False)


#===============================================================================
#
#  POPUP for the APPLY procedure
#
def create_popup_apply(parent):
    GUI_Global.popup_apply = wx.Dialog(parent, wx.NewId(), 'APPLY', 
        size=(250, 130), style=wx.STAY_ON_TOP | wx.DEFAULT_DIALOG_STYLE)
    GUI_Global.popup_apply.SetBackgroundColour((255,20,20))
    text = wx.StaticText(GUI_Global.popup_fit, -1, 
        '      ...radiative transfer code\n                           in execution...\n            GUI blocked  ', 
        pos=(5,3))
    text.SetForegroundColour((0,0,0))
    text.Show(True)
    GUI_Global.apply_available = False

    #---------------------------------------------------------------------------
    # ABORT BUTTON   
    GUI_Global.popup_apply.button_abort = wx.Button(
        GUI_Global.popup_apply, 
        id=-1, label="ABORT", pos=(80,70))
    GUI_Global.popup_apply.Bind(wx.EVT_BUTTON, OnButton_abort, 
        GUI_Global.popup_apply.button_abort ) 
    GUI_Global.popup_apply.button_abort.SetToolTip(wx.ToolTip(
        "click to if you are ready want to abort the call"))
    GUI_Global.popup_apply.button_abort.Show(True)

def show_popup_apply(parent):
    # hide the old window
    GUI_Global.popup_aplot.Show(False)
    if GUI_Global.data_available == False:
        no_data_message()
        return
    if GUI_Data.check_output() == False:
        return
    # bad hack for the matplotlib problem changing axis with new wavelength
    # range
    GUI_Global.popup_aplot.Destroy()
    GUI_Plot.create_popup_apply_plot()
    # end of hack

    file_fit = str(GUI_Data.output_dir_path() + '/' \
            + GUI_Global.output_name_v + 
                '_fit.fits')
    file_tac = str(GUI_Data.output_dir_path() + '/' \
            + GUI_Global.output_name_v + 
                '_tac.fits')
    if os.path.exists(file_fit) == False:
        dlg = wx.MessageDialog( GUI_Global.frame, 
            "No fit result available - run FIT first", "ERROR", wx.OK)
        retval = dlg.ShowModal() 
        dlg.Destroy()    
        return  
    
    if ((GUI_Global.param_changes == True) or
        (os.path.getmtime(GUI_Global.Current_Setup_File) >
         os.path.getmtime(file_fit))):
        #print GUI_Global.param_changes
        #print os.path.getmtime(GUI_Global.Current_Setup_File)
        #print os.path.getmtime(file_fit)
        dlg = wx.MessageDialog( GUI_Global.frame, 
            "Parameters have been changed - run FIT first", "WARNING", wx.OK)
        retval = dlg.ShowModal() 
        dlg.Destroy()    
        return  

    # check if a corrected file is available that is younger than the file of
    # the fit
    if (os.path.exists(file_fit) and os.path.exists(file_tac)):
        if ((os.path.getmtime(file_fit) < os.path.getmtime(file_tac)) and
            (os.path.getmtime(GUI_Global.Current_Setup_File) < 
            os.path.getmtime(file_tac))):
            EndBusyCursor()
            # no correction neeeded - exists already
            dlg = wx.MessageDialog(GUI_Global.frame, TextPool.dialog_012,
                GUI_Global.VERSION, wx.YES_NO | wx.ICON_QUESTION | 
                wx.NO_DEFAULT)
            result = dlg.ShowModal()             
            dlg.Destroy()                        
            if result == wx.ID_NO:
                GUI_Data.input_data(file_tac,False,False,True)
                GUI_Global.apply_available = True
                GUI_Plot.draw_apply_figure()
                GUI_Global.popup_aplot.Show(True)
                return  
            

    if GUI_Global.param_changes == True:
        dlg = Three_Button_Message_Dialogue(GUI_Global.frame,
               "Parameters have been changed!", 
               "ignore changes", "save setup file", "cancel")
        EndBusyCursor()
        retval = dlg.ShowModal() 
        dlg.Destroy()   
        if (retval == 3) :
            
            return
        if (retval == 2) :
            EndBusyCursor()
            GUI_Menu.OnSave(parent)

    GUI_Global.apply_output = str(GUI_Data.output_dir_path() + '/' + 
        GUI_Global.output_name_v + '_apply_out.txt')
    # run the fit in a terminal window
    cmd = "%s %s" % (GUI_Global.calctrans_exe, 
        GUI_Global.Current_Setup_File)
    wx.BeginBusyCursor(wx.HOURGLASS_CURSOR)
    retval = run_subcommand_in_popup(cmd, GUI_Global.apply_output, True)
    EndBusyCursor()
    if GUI_Global.user_killed == True:
        GUI_Console.hide_console(None)
        return

    if retval != 0:
        dlg = wx.MessageDialog( GUI_Global.frame, 
            "ERROR: calctrans failed", "ERROR", wx.OK)
        retval = dlg.ShowModal() 
        dlg.Destroy()   
        #print "error in prepearing the file"
    else:  
        prepearefile = str( GUI_Data.output_dir_path() + '/' + \
                GUI_Global.output_name_v + '_tac' + '.fits')
        #print "RESULT: %s" % prepearefile
        if os.path.exists(prepearefile) == True:
            GUI_Data.input_data(prepearefile,False,False,True)
            GUI_Global.apply_available = True
            GUI_Plot.draw_apply_figure()
            GUI_Global.popup_aplot.Show(True)
        else:
            dlg = wx.MessageDialog( GUI_Global.frame, 
                "ERROR: Result file not found", "ERROR", wx.OK)
            retval = dlg.ShowModal() 
            dlg.Destroy()   
            #print "result file not found"
    
    EndBusyCursor()
    #GUI_Global.popup_apply.Show(False)


#===============================================================================
#
#  variuos routines for the popups
#
def no_data_message():
    EndBusyCursor()
    dlg = wx.MessageDialog( GUI_Global.frame, 
		 TextPool.dialog_008,
	     GUI_Global.VERSION, wx.OK)
    dlg.ShowModal() # Show it
    dlg.Destroy()   # finally destroy it when finished with OK

def no_parameter_message():
    EndBusyCursor()
    dlg = wx.MessageDialog( GUI_Global.frame, 
		 TextPool.dialog_009,
	     GUI_Global.VERSION, wx.OK)
    dlg.ShowModal() # Show it
    dlg.Destroy()   # finally destroy it when finished with OK

def show_popup_mask(parent):
    EndBusyCursor()
    if GUI_Global.data_available == False:
        no_data_message()
        return
    GUI_Global.popup_mask.Show(True)




#===============================================================================
#
# Own 3 button message dialogue as it is not fully standard yet in all 
# implimentations / distributions
#
#class Three_Button_Message_Dialogue(wx.Frame):
#    
#    def __init__(self, *args, **kwargs):            
#        self.InitUI()
#        
#    def InitUI(self):    
def Three_Button_Message_Dialogue(parent, Message,
        Button_1_Text, Button_2_Text, Button_3_Text):
        new_popup = wx.Dialog(parent, wx.NewId(), 
        'Message') 
        #size=(450, 750))
        panel = wx.Panel(new_popup)
        topSizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(topSizer)

        text = "       " + Message
        st = wx.StaticText(panel,-1,text)#,pos=(30,30))
        st.SetForegroundColour((255,0,0))

        topSizer.AddSpacer(5)
        topSizer.Add(st,  proportion=wx.ALIGN_CENTER)
        topSizer.AddSpacer(5)

        #hbox = wx.BoxSizer()
        sizer = wx.GridSizer(1, 3, 5, 5)

        btn1 = wx.Button(panel, label=Button_1_Text)
        btn2 = wx.Button(panel, label=Button_2_Text)
        btn3 = wx.Button(panel, label=Button_3_Text)

        #sizer.AddSpacer(5)
        sizer.AddMany([btn1, btn2, btn3])
        #sizer.AddSpacer(5)

        topSizer.Add(sizer)
        #topSizer.SetSizer(sizer)
        
        #panel.SetSize(topSizer.GetSize())
        
        #hbox.Fit()
        #vbox.SetSizer(hbox)

        #sizer1.Add(hbox)
        #vbox.Add(sizer1)
        

        btn1.Bind(wx.EVT_BUTTON, Action1)
        btn2.Bind(wx.EVT_BUTTON, Action2)
        btn3.Bind(wx.EVT_BUTTON, Action3)

        #panel.SetSize((300, 200))
        new_popup.SetTitle('Message')
        new_popup.Centre()
        topSizer.Fit(new_popup)
        return new_popup

def Action1(event):
    a = event.GetEventObject()
    b = a.GetParent()
    c = b.GetParent()
    c.SetReturnCode(1)
    c.EndModal(1)
    #return 1

def Action2(event):
    a = event.GetEventObject()
    b = a.GetParent()
    c = b.GetParent()
    c.SetReturnCode(2)
    c.EndModal(2)
    #return 2

def Action3(event):
    a = event.GetEventObject()
    b = a.GetParent()
    c = b.GetParent()
    c.SetReturnCode(3)
    c.EndModal(3)
    #return 3

def run_subcommand_hidden(cmd, logfile):
    cmd_grp = cmd.split()
    #print cmd_grp
    if (logfile != None):
        # open the logfile in append mode and give it as parameter for the 
        # real subprocess command, stderr is linked into stdout
        log_fd = open(logfile,"w")
        p_cmd = subprocess.Popen(cmd_grp, shell=False, stderr=subprocess.STDOUT, 
                             stdout=log_fd)
        # wait for the subprocess to terminate and then close the logfile and
        # read the return status code
        p_cmd.wait()
        log_fd.close()
    else: 
        log_fd = open("/dev/null","a")
        p_cmd = subprocess.Popen(cmd_grp, shell=False, stderr=subprocess.STDOUT, 
                             stdout=log_fd)
        # wait for the subprocess to terminate and then close the logfile and
        # read the return status code
        p_cmd.wait()
        log_fd.close()
        
    retcode = p_cmd.returncode
    #print 'retcode ' + retcode
    return retcode

def my_timer_proc():  
    if p_cmd.poll() != None:
        # hide popup
        subcommand_busy = False 
        GUI_Global.popup_fit.EndModal(0)
    else:
        # respan this time
        wx.FutureCall(2000, my_timer_proc)

def run_subcommand_in_popup(cmd, logfile, popup):
    global subcommand_busy
    global p_cmd
    subcommand_busy = False
    GUI_Global.user_killed = False
    cmd_grp = cmd.split()
    os.system('echo `date` > '+logfile)
    # open the logfile in append mode and give it as parameter for the 
    # real subprocess command, stderr is linked into stdout
    log_fd = open(logfile,"a")
    p_cmd = subprocess.Popen(cmd_grp, shell=False, stderr=subprocess.STDOUT, 
                             stdout=log_fd)

    GUI_Global.child = p_cmd
    
    GUI_Console.show_console(logfile,True)
    # wait for the subprocess to terminate and then close the logfile and
    # read the return status code
    if popup == False:
        # normal call
        p_cmd.wait()
    else:
        # call with a popup window
        subcommand_busy = True
        # call a timer to check the code running
        wx.FutureCall(2000, my_timer_proc)
        # show in Modal mode a popup
        GUI_Global.popup_fit.ShowModal()

    log_fd.close()
    retcode = p_cmd.returncode
    # terminate the popup screen logfile
    #p_echo.terminate()
    GUI_Console.stop_refresh_console()
    return retcode

