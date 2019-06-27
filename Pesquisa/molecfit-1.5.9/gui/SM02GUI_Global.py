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
global VERSION

import os
import wx

# class to store properties of a molecule
class Molec:
    def __init__(self, name, abundance, fit, calc):
        self.name   = name
        self.ab     = abundance
        self.fit    = fit
        self.calc   = calc
        self.name_t = None
        self.ab_t   = None
        self.fit_t  = None
        self.calc_t = None

# class to store properties of a mask or fit region
class Mask:
    def __init__(self, active, include, startx, endx, validxy, startwl, endwl):
        self.active  = active       # flag for activity (True/False)
        self.include = include      # include (True) or exclude (False) mask 
        self.startx  = startx       # start/end in pixel value - valid only if
        self.endx    = endx         #   validxy == True
        self.validxy = validxy      # if True - startx/endx in pixels is valid
        self.startwl = startwl      # wavelength start/end - valid only if 
        self.endwl   = endwl        #   validxy == False
        self.box     = None         # widget entry for the mask drawing box
        self.line    = None         # widget entry for the mask drawing line

global NR_OF_MOLS                   # these two numbers define the maximum
global NR_OF_MAIN_MOLS              # number of molecules and that ones in 
                                    # a free fit mode

global frame      # main frame of the window as "parent" for most popups
global panel      # drawing panel withn main frame
global maskpanel  # drawing panel within masking tool

# Popup windows for the submenues
global popup_mask
global popup_back
global popup_molec
global popup_mask
global popup_wavelength
global popup_kernel 
global popup_continuum 
global popup_back 
global popup_path 
global popup_other 
global popup_parameters 

#-------------------------------------------------------------------------------
# parameters for the fitting code 
# they exist two times - the ending '_v' stands for the value option
# the '_t' is the text field in the graphical widget
#

