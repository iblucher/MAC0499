#!/usr/bin/python -B
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
import os
import SM02GUI_Global      as GUI_Global
import SM02GUI_Data        as GUI_Data
import sys

class PWVParameter(object):
    def __init__(self, value):
        self._value = value
        self.chkctrl = None
        self.txtctrl = None

    def __repr__(self):
        if self._value < 0:
            return 'pwv: -1.'
        else:
            return 'pwv: %r' % self._value

    def refresh(self):
        if self.chkctrl is not None:
            self.chkctrl.SetValue(self._value >= 0.)
        if self.txtctrl is not None:
            self.txtctrl.SetValue(repr(self._value))

    @property
    def pwv(self):
        if self.chkctrl is None:
            return self._value
        if self.chkctrl.IsChecked():
            return self._value
        else:
            return  -1.

    @pwv.setter
    def pwv(self, value):
        if value < 0:
            self._value = -1.
        else:
            self._value = value
        self.refresh()

    @property
    def disabled(self):
        return self._value < 0.

    @disabled.setter
    def disabled(self, d):
        if d:
            self.pwv = -1.

#===============================================================================
# this procedure defines some minimum environment and searches for the 
# molecfit executable. It defines as default working directory that one 
# were the user starts this GUI
#
def is_exe(fpath):
    # combined macro for path exists + executable
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def HideMessage():
    # time controlled function to hide message after some timeout
    if GUI_Global.cur_dlg != None:
        GUI_Global.cur_dlg.EndModal(0)
        GUI_Global.cur_dlg.Destroy()
        GUI_Global.cur_dlg = None;        
  
def check_installation_basedir_executables():
    # create the names and test the existence of other important executables
    base = GUI_Global.basedir_v
    GUI_Global.molecfit_exe = os.path.join(base, 'bin', 'molecfit')
    GUI_Global.prepgui_exe = os.path.join(base, 'bin', 'prepguitable')
    GUI_Global.calctrans_exe = os.path.join(base, 'bin', 'calctrans')
    lblrtm = os.path.join(base, 'bin', 'lblrtm')
    lnfl = os.path.join(base, 'bin', 'lnfl')
    for e in (GUI_Global.molecfit_exe, GUI_Global.prepgui_exe,
              GUI_Global.calctrans_exe, lblrtm, lnfl):
        if not is_exe(e):
            text = 'Can\'t find executable %s\nwith a basedir %s\n\n' \
                   'Check installation !' %  (e, GUI_Global.basedir_v)
            raise ValueError(text)

    if hasattr(GUI_Global,"output_dir_v") == False:
        GUI_Global.output_dir_v = 'output'
        
    if hasattr(GUI_Global,"output_name_v") == False:
        GUI_Global.output_name_v = 'default_result'
    return True

def init_path_and_environment():
    # basic init procedure
    GUI_Global.cur_dlg = None
    GUI_Global.rescale = 1.0;

    GUI_Global.param_changes = True
    GUI_Global.mask_changes  = True
    GUI_Global.mask_min = -100000000
    GUI_Global.mask_max = 200000000
    GUI_Global.param_available = False # this blocks most actions until
                                       # parameters are filled with defaults
    GUI_Global.data_available  = False # this blocks most actions until
    GUI_Global.fit_available   = False # data / file has been read
    GUI_Global.apply_available = False 

    GUI_Global.corrected_available = False
    GUI_Global.mask_exist      = False
    GUI_Global.mask            = []
    GUI_Global.NR_OF_MOLS      = 23 # these two numbers define the maximum
    GUI_Global.NR_OF_MAIN_MOLS = 7  # number of molecules and those in a
                                        # possible fit
    GUI_Global.Current_Setup_File = '*** not defined yet ***'
    # default WORKING directory at start
    GUI_Global.cwd       = os.getcwd()
    GUI_Global.workdir_v = GUI_Global.cwd.strip()
    GUI_Global.pwv_par = PWVParameter(-1.)

    # default MOLECFIT basedir 
    # first use the ENVIRONMENT variable
    # if os.getenv('MOLECFIT_BASEDIR') != None:
    #     # ENVIRONMANT was set - priority 1
    #     GUI_Global.basedir_v = \
    #       os.path.expandvars(os.path.expanduser(os.getenv('MOLECFIT_BASEDIR')))
    #     GUI_Global.basedir_v.rstrip('/')
    # else:
        # assume installation tree
    GUI_Global.basedir_v = os.path.split(os.path.split(sys.argv[0])[0])[0]

# TODO move to after parameter load
    check_installation_basedir_executables()

    # Moved to Paramfile
    # if os.path.exists(GUI_Data.output_dir_path()) == False:
    #     os.system('mkdir -p ' + GUI_Data.output_dir_path())
    

    return True

def make_absolute_path(arg):
    # remove possible ../ and ./ prefixes
    if arg[:1] != '.':
        return arg
    cwd = os.getcwd()
    if arg[:2] == './':
        return cwd + arg[1:]
    return os.path.abspath(arg)
