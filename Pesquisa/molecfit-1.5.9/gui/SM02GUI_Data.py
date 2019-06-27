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

try:
    import pyfits
except ImportError:
    import astropy.io.fits as pyfits
    
import numpy as NumPy

import SM02GUI_Global as GUI_Global
import SM02GUI_Plot   as GUI_Plot
import SM02GUI_Mask   as GUI_Mask
import SM02GUI_Popups as GUI_Popups


def output_dir_path():
    
    # generates a homogenious output directory system for the underlying molecfit code
    if GUI_Global.output_dir_v[:1] == '/':
        return GUI_Global.output_dir_v
    else:
        return (GUI_Global.workdir_v + '/' + GUI_Global.output_dir_v)


def check_output():
    
    # check and repair (if possible) the output directory permissions
    path = output_dir_path()
    if os.path.exists(path) == False:
        text = 'Output directory\n%s\ndoes not exits.\n'\
            'Should I create it now ?' % path
        dlg = wx.MessageDialog( GUI_Global.frame, 
		        text,
	            GUI_Global.VERSION, 
                wx.YES_NO | wx.ICON_QUESTION | wx.NO_DEFAULT)
        if dlg.ShowModal() != wx.ID_NO:
            dlg.Destroy() 
            os.mkdir(path,0600)
            return True
        else:
            dlg.Destroy() 
            return False 

    if os.path.isdir(path) == False:
        text = '\n%s\nexits - but is not a directory.\n'\
            'Check Setup !' % path
        dlg = wx.MessageDialog( GUI_Global.frame, 
		        text,
	            GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() 
        dlg.Destroy() 
        return False
    return True
  
def empty_data_define():
    # create an empty dummy data set to be able to start matplotlib
    GUI_Global.X  = NumPy.linspace(0,0.1, 5,endpoint=True)
    GUI_Global.Y  = NumPy.linspace(0,0.1, 5,endpoint=True)
    GUI_Global.XF = GUI_Global.X
    GUI_Global.YF = GUI_Global.Y
    GUI_Global.XD = GUI_Global.X
    GUI_Global.YD = GUI_Global.Y
    GUI_Global.XR = GUI_Global.X
    GUI_Global.YR = GUI_Global.Y
    GUI_Global.valid_px = len(GUI_Global.X) - 1


def min_mean_median_max_of_data(array):
    a = NumPy.asarray(array)
    # data statistics tool
    m = [0.0] * 5
    m[0] = a.min()
    m[1] = a.mean()
    m[2] = NumPy.median(a)
    m[3] = a.max()
    m[4] = a.std()
    return m

def pixel_of_data(wavelength):
    # calculate the pixel number in file from wavelength table
    if GUI_Global.data_available == False :
        return None
    pos = NumPy.where(GUI_Global.X >= wavelength)[0]
    if pos.size != 0:
        return pos[0]
    else:
        return GUI_Global.X.size - 1

def value_of_data_from_wave(wavelength):
    # derive the intensity/flux value most nearby to a given wavelength
    if GUI_Global.data_available == False :
        return None
    pos = NumPy.where(GUI_Global.X >= wavelength)[0]
    if pos.size != 0:
        return GUI_Global.Y[pos[0]]
    else:
        return GUI_Global.Y[GUI_Global.X.size - 1]

def calculate_residuum():
    # and calculate the residuum
    GUI_Global.XD = (GUI_Global.X+GUI_Global.XF)/2
    GUI_Global.YD = (GUI_Global.Y-GUI_Global.YF)
    # mask the non fitted regions
    GUI_Global.YD[GUI_Global.YF == 0] = 0

def input_data(filename, refresh, model, corrected):
    # read a FITS data file created / prepeared by prepguitable
    # boolean refresh controls if the displays are redrawn completely
    # boolean model & corrected selects if it is a raw data set, a model
    # or a file were the correction has been applied already
    GUI_Global.data_available = True 
    hdulist = pyfits.open(filename)
    tbdata = hdulist[1].data
    if ((model == False) and (corrected == False)):
        # always activate show data if required to obtain new HOME window
        # size in matplotlib
        GUI_Global.panel.cb_raw.SetValue(True)
        GUI_Global.X  = tbdata.field('lambda')[:]
        GUI_Global.Y  = tbdata.field('flux')[:]
        GUI_Global.fit_available == False
        GUI_Global.valid_px = len(GUI_Global.X) - 1
    else:
        if (model == True):
            GUI_Global.XF = tbdata.field('lambda')[:]
            GUI_Global.YF = tbdata.field('mflux')[:]
            GUI_Global.fit_available == True
            if (GUI_Global.rescale < 0.5):
                GUI_Global.YF /= GUI_Global.rescale
            calculate_residuum()
        else:
            # print "read CORRECTED data file %s" % filename
            GUI_Global.XC = tbdata.field('lambda')[:]
            GUI_Global.YC = tbdata.field('cflux')[:]
            GUI_Global.corrected_available == True
            if (GUI_Global.rescale < 0.5):
                GUI_Global.YC /= GUI_Global.rescale           

    hdulist.close()
    # these del calls should be better - but crashing pyfits on
    # old Fedora15 systems. No broblem on Fedora17-19 (pyfits > v2.5)
    #
    # del hdulist[1].data
    # del hdulist[0].data
    # del tbdata
    GUI_Global.maskpanel.first = True
    if refresh == True:
        GUI_Plot.draw_figure(False)
        GUI_Mask.draw_figure(False)
    text = "Currently displayed file:                     %s" % \
        GUI_Global.filename_v
    GUI_Global.current_plot_file_t.SetLabel(text)

def read_input_file_pixel_exclusion_mask(filename):
    # reads, if available, the pixel excusion masks from the FITS file
    hdulist = pyfits.open(filename)
    tbdata = hdulist[1].data
    XM  = tbdata.field('lambda')[:]
    YM  = tbdata.field('mask')[:]

    # Create an alternative mask for bad pixels #
    # print 'DEBUG Alternative bad pixel mask'
    # Pad the array first
    lendata = len(YM)
    badpix = NumPy.ones(lendata + 2)
    badpix[1:-1] = YM
    # Get the boundaries 1 and -1 by construction
    indexes = NumPy.roll(badpix, 1) - badpix
    starts = NumPy.where(indexes == 1)[0]
    ends = NumPy.where(indexes == -1)[0] -1
    # Un-pad indexes
    starts -= 1
    ends -= 1
    # Check for inconsistencies
    if NumPy.any(NumPy.logical_or(starts < 0, starts > lendata - 1)) or \
       NumPy.any(NumPy.logical_or(ends < 0, ends > lendata - 1)):
       print 'ERROR: start/end indices of bad pixels out of range'
       return
    # Get the wavelengths
    badpixregions = [{'start': XM[s], 'end': XM[e]} for s, e in zip(starts, ends)]
    # print 'DEBUG number of bad pixel regions:', len(badpixregions)
    # import json
    # print json.dumps(badpixregions, indent=4)
    GUI_Global.badpixregions = badpixregions
    #############################################

    # if hasattr(GUI_Global, "mask") == False:
    #     GUI_Global.mask = []
    # maskOn = False
    # added_one = False
    # for ii in range(len(XM)):
    #     if maskOn == False:
    #         if YM[ii] == 0:
    #             startX = ii
    #             startW = XM[startX]
    #             maskOn = True
    #     else:
    #         if YM[ii] == 1:
    #             endX   = ii-1
    #             endW   = XM[endX]
    #             maskOn = False
    #             added_one = True
    #             GUI_Global.mask.append(
    #                 GUI_Global.Mask(True, False, startX, endX, True, 
    #                     startW, endW))

    # if added_one == True:
    #     if hasattr(GUI_Global, "prange_exclude_v") == False:
    #         GUI_Global.prange_exclude_v = GUI_Global.workdir_v + '/' + \
    #             GUI_Global.output_name_v + '_prange_exclude.dat'
    #         os.system('echo > '+GUI_Global.prange_exclude_v)
    #     if GUI_Global.prange_exclude_v == 'none':
    #         GUI_Global.prange_exclude_v = GUI_Global.workdir_v + '/' + \
    #             GUI_Global.output_name_v + '_prange_exclude.dat'
    #         os.system('echo > '+GUI_Global.prange_exclude_v)
    #     GUI_Mask.MergeMasksOfSameType()
        

def import_mask_files(filename, include, wavelength):
    # import a mask file of command line version of molecfit
    #   include    = True ... inclusion mask (ignores wavelength parameter)
    #   wavelength = True ... file is in wl mode - else in pixel mode

    if filename == 'none':
        return
    if hasattr(GUI_Global, "mask") == False:
        GUI_Global.mask = []
    if os.path.exists(filename) != True or os.access(filename,os.R_OK) != True:
        if include == True:
            text = 'ERROR: File \n%s\n for "Fit Regions (inclusion)" does '\
                   'not\n exist or is not readable' % filename
        else:
            if wavelength == True:
                text = 'ERROR: File \n%s\n for "Wavelength Masks (exclusion)" '\
                   'does not\n exist or is not readable' % filename
            else:
                text = 'ERROR: File \n%s\n for "Pixel Masks (exclusion)" '\
                   'does not\n exist or is not readable' % filename
        dlg = wx.MessageDialog( GUI_Global.frame, 
		        text,
	            GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() 
        dlg.Destroy() 
        del text
        return
        
    fd = open( filename, 'r' )
    ii = 0
    array = []
    for line in fd:
        if len(line) > 2:
            # Skip comments
            if line.strip()[0] == '#':
                continue
            array.append( line )
            line_arr = array[ii].split()
            if include == True:
                wstart   = float(line_arr[0])
                wend     = float(line_arr[1])
                GUI_Global.mask.append(
                    GUI_Global.Mask(True, True, 0.0, 0.0, False, 
                        wstart, wend))
            else:
                if wavelength == True:
                    wstart = float(line_arr[0])
                    wend   = float(line_arr[1])
                    GUI_Global.mask.append(
                        GUI_Global.Mask(True, False, 0.0, 0.0, False, 
                        wstart, wend))
                else:
                    # FORTRAN+FITS / C+python hack - a pixel mask starting at 1  
                    # means 0 internally in the code
                    xstart = int(line_arr[0])-1
                    xend   = int(line_arr[1])-1
                    
                    if (xend > (len(GUI_Global.X) - 1)):
                        xend = (len(GUI_Global.X) - 1)
                        if (xstart > xend):
                            xstart = xend - 1
                    GUI_Global.mask.append(
                        GUI_Global.Mask(True, False, xstart, xend, True, 
                            0.0, 0.0))
            ii += 1
    fd.close()
    GUI_Mask.MergeMasksOfSameType()
    # don't trigger a save dialog even if masks could be merged
    GUI_Global.mask_changes = False
    GUI_Plot.draw_figure(True)
    GUI_Mask.draw_figure(True)

def prepeare_fits_file(data_file):
    # calling prepguitable with a temporary setup file. The latter is required
    # as molecfit (=prepguitable is the first part of molecfit) requires
    # such a file, but the user may load a fits file in the GUI without 
    # already having a (meaningful) setup. This feature had to be added due
    # to JIRA ticket PIPE-4633 without completely breaking with the 
    # 'philosophy' of the use of molecfit
    import tempfile
    if (hasattr(GUI_Global,'wlgtomicron_v') != True):
        GUI_Global.wlgtomicron_v = 1e-3


    filen = tempfile.NamedTemporaryFile(mode='w')
    # filen.writelines('basedir: '     + GUI_Global.basedir_v + '\n')
    filen.writelines('filename: '    + data_file + '\n')
    filen.writelines('columns: '     + GUI_Global.col_lam_v + ' '  \
                                     + GUI_Global.col_flux_v + ' ' \
                                     + GUI_Global.col_dflux_v + ' ' \
                                     + GUI_Global.col_mask_v + '\n')
    filen.writelines('output_dir: '  + output_dir_path() +'\n')
    filen.writelines('output_name: ' + GUI_Global.output_name_v + '\n')
    text = '%.2g' % (GUI_Global.wlgtomicron_v)
    filen.writelines('wlgtomicron: ' + text + '\n')
    filen.writelines('vac_air: '     + GUI_Global.vac_air_v + '\n')
    filen.writelines('end\n')  
    filen.flush()
    cmd = GUI_Global.prepgui_exe + ' ' + filen.name
    GUI_Global.prepguitable_output = str('%s/%s_prepguitable_out.txt' % (
        output_dir_path(),GUI_Global.output_name_v))
    ret = GUI_Popups.run_subcommand_hidden(cmd,GUI_Global.prepguitable_output)
    filen.close()
    return ret
