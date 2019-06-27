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
import SM02GUI_Plot        as GUI_Plot
import SM02GUI_Popups      as GUI_Popups
import SM02GUI_Data        as GUI_Data
import SM02GUI_Init        as GUI_Init

import SM02GUI_Help_Dialog_ToolTip_Texts as TextPool

def strncmp(s1, s2, n):
    if s1[:n] == s2[:n]:
        return True
    return False
    
def refresh_all_popups_on_parameter_change():
    GUI_Popups.refresh_popup_molec()
    GUI_Popups.refresh_popup_wavelength()
    GUI_Popups.refresh_popup_kernel()
    GUI_Popups.refresh_popup_cont()
    GUI_Popups.refresh_popup_background()
    GUI_Popups.refresh_popup_path() 
    GUI_Popups.refresh_popup_other() 

def my_message(text):
    dlg = wx.MessageDialog(GUI_Global.frame, text,
          GUI_Global.VERSION, wx.OK | wx.ICON_QUESTION)
    dlg.ShowModal()
    dlg.Destroy()        

def check_path_and_files():
    if os.path.exists(GUI_Global.workdir_v) == False:
        text = ('WORKDIR: \n\nDirectory \n%s\n not found!') % \
            (GUI_Global.workdir_v)
        my_message(text)    
        return False 
    if os.access(GUI_Global.workdir_v,os.W_OK) == False:
        text = ('WORKDIR: \n\nNo write permission in directory \n%s') % \
            (GUI_Global.workdir_v)    
        return False 
    if GUI_Global.filename_v[:6] != '** not':
        # relative path is relative to current directory, remove this crap
        # if GUI_Global.filename_v[:1] != '/':
        #     GUI_Global.filename_v = GUI_Global.basedir_v + '/' + \
        #         GUI_Global.filename_v 
        if os.path.exists(GUI_Global.filename_v) == False:
            text = ('Spectrum (FITS): \n\nFile \n%s\nnot found!') % \
                (GUI_Global.filename_v)    
            my_message(text)    
            return False 

    if GUI_Global.wrange_include_v[:4] != 'none':
        if GUI_Global.wrange_include_v[:1] != '/':
            path = GUI_Global.workdir_v + '/' + GUI_Global.wrange_include_v
        else:
            path = GUI_Global.wrange_include_v
        if os.path.exists(path) == False:
            # generate empty one if required
            open(path, 'w').close()
    if GUI_Global.wrange_exclude_v[:4] != 'none':
        if GUI_Global.wrange_exclude_v[:1] != '/':
            path = GUI_Global.workdir_v + '/' + GUI_Global.wrange_exclude_v
        else:
            path = GUI_Global.wrange_exclude_v
        if os.path.exists(path) == False:
            # generate empty one if required
            open(path, 'w').close()
    if GUI_Global.prange_exclude_v[:4] != 'none':
        if GUI_Global.prange_exclude_v[:1] != '/':
            path = GUI_Global.workdir_v + '/' + GUI_Global.prange_exclude_v
        else:
            path = GUI_Global.prange_exclude_v
        if os.path.exists(path) == False:
            # generate empty one if required
            open(path, 'w').close()
    if os.path.exists(GUI_Data.output_dir_path()) \
        == False:
        os.mkdir(GUI_Data.output_dir_path())

    if os.access(GUI_Data.output_dir_path(),os.W_OK) \
        == False:
        text = ('OUTPUT: \n\nNo write permission in output directory \n%s') % \
            (GUI_Data.output_dir_path())    
        return False 
    return True

# def absolute_path():
    # if workdir != basedir we have to use absolute pathes
    # This is total non-sense, remove this crap
    # if GUI_Global.basedir_v != GUI_Global.workdir_v:
    #     if hasattr(GUI_Global,"wrange_include_v") == False:
    #         GUI_Global.wrange_include_v = 'none'
    #     if  GUI_Global.wrange_include_v[:1] != '/' and \
    #         GUI_Global.wrange_include_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.workdir_v + '/' + \
    #                 GUI_Global.wrange_include_v) == True:
    #             GUI_Global.wrange_include_v = \
    #                 GUI_Global.workdir_v + '/' + GUI_Global.wrange_include_v
    #         else:
    #             GUI_Global.wrange_include_v = \
    #                 GUI_Global.basedir_v + '/' + GUI_Global.wrange_include_v
    #     if hasattr(GUI_Global,"wrange_exclude_v") == False:
    #         GUI_Global.wrange_exclude_v = 'none'
    #     if GUI_Global.wrange_exclude_v[:1] != '/' and \
    #        GUI_Global.wrange_exclude_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.workdir_v + '/' + \
    #                 GUI_Global.wrange_exclude_v) == True:
    #             GUI_Global.wrange_exclude_v = \
    #                 GUI_Global.workdir_v + '/' + GUI_Global.wrange_exclude_v
    #         else:
    #             GUI_Global.wrange_exclude_v = \
    #                 GUI_Global.basedir_v + '/' + GUI_Global.wrange_exclude_v
    #     if hasattr(GUI_Global,"prange_exclude_v") == False:
    #         GUI_Global.prange_exclude_v = 'none'
    #     if GUI_Global.prange_exclude_v[:1] != '/' and \
    #        GUI_Global.prange_exclude_v[:4] != 'none':
    #         if os.path.exists(GUI_Global.workdir_v + '/' + \
    #                 GUI_Global.prange_exclude_v) == True:
    #             GUI_Global.prange_exclude_v = \
    #                 GUI_Global.workdir_v + '/' + GUI_Global.prange_exclude_v
    #         else:
    #             GUI_Global.prange_exclude_v = \
    #                 GUI_Global.basedir_v + '/' + GUI_Global.prange_exclude_v

# def patch_basedir(filename):
    # patches the basedir parameter without (!) chaning the creation date/time
    # of the setup file
    # first make a copy with the original date & time 
    # and dereference content if it is a logical link
    # This is scary...
    # cmd = 'cp -fpL %s %s.bak' % (filename,filename)
    # os.system(cmd)
    # #print cmd
    # fd = open( filename, 'r' )
    # array = []
    # for line in fd:
    #     array.append( line )
    # fd.close()
    # fd = open( filename, 'w' )
    # for ii in range(len(array)):
    #     if strncmp('basedir:', array[ii], 8):
    #         basedir_v = array[ii][8:]
    #         basedir_v = basedir_v.strip()
    #         if basedir_v == '.':
    #             fd.writelines('basedir: '+ GUI_Global.basedir_v + '\n')
    #         else:
    #             fd.writelines(array[ii])
    #     elif strncmp('plot_creation:', array[ii], 14):
    #         fd.write('plot_creation: none' + '\n')
    #     #elif strncmp('output_dir:', array[ii], 11):
    #     #    GUI_Data.check_output()
    #     #    fd.write('output_dir: ' + GUI_Global.workdir_v + '/' + \
    #     #        GUI_Global.output_dir_v + '\n')
    #     elif GUI_Global.basedir_v != GUI_Global.workdir_v:
    #         if strncmp('wrange_include:', array[ii], 15):
    #             if GUI_Global.wrange_include_v != 'none' and \
    #                 GUI_Global.wrange_include_v[:1] != '/':
    #                 GUI_Global.wrange_include_v = GUI_Global.basedir_v + '/' + \
    #                                               GUI_Global.wrange_include_v
    #             fd.write('wrange_include: ' +GUI_Global.wrange_include_v + '\n')
    #         elif strncmp('wrange_exclude:', array[ii], 15):
    #             if GUI_Global.wrange_exclude_v != 'none' and \
    #                 GUI_Global.wrange_exclude_v[:1] != '/':
    #                 GUI_Global.wrange_exclude_v = GUI_Global.basedir_v + '/' + \
    #                                               GUI_Global.wrange_exclude_v
    #             fd.write('wrange_exclude: ' +GUI_Global.wrange_exclude_v + '\n')
    #         elif strncmp('prange_exclude:', array[ii], 15):
    #             if GUI_Global.prange_exclude_v != 'none' and \
    #                 GUI_Global.prange_exclude_v[:1] != '/':
    #                 GUI_Global.prange_exclude_v = GUI_Global.basedir_v + '/' + \
    #                                               GUI_Global.prange_exclude_v
    #             fd.write('prange_exclude: ' +GUI_Global.prange_exclude_v + '\n')
    #         elif strncmp('output_dir:', array[ii], 11):
    #             fd.write('output_dir: ' + GUI_Data.output_dir_path() + '\n')
    #         elif strncmp('filename:', array[ii], 9):
    #             if GUI_Global.filename_v[:1] != '/':
    #                 GUI_Global.filename_v = GUI_Global.basedir_v + '/' + \
    #                                               GUI_Global.filename_v
    #             fd.write('filename: ' + GUI_Global.filename_v + '\n')
    #         else:
    #             fd.writelines(array[ii])
    #     else:
    #         fd.writelines(array[ii])
    # fd.close()
    # cmd = 'touch -r %s.bak %s' % (filename,filename)
    # os.system(cmd)
   

def read_param_file(filename):
    if os.path.exists(filename) == False:
        text = ('Parameter file %s not found!') % (filename)
        dlg = wx.MessageDialog(GUI_Global.frame, text,
            GUI_Global.VERSION, wx.OK | wx.ICON_QUESTION)
        result = dlg.ShowModal()             
        dlg.Destroy()        
        return False

    if os.access(filename,os.W_OK) == False:
        text = ('Parameter file %s has no write permission!\n'\
            'Make it writeable or use a local copy.') % (filename)
        dlg = wx.MessageDialog(GUI_Global.frame, text,
            GUI_Global.VERSION, wx.OK | wx.ICON_QUESTION)
        result = dlg.ShowModal()             
        dlg.Destroy()        
        return False

     
    GUI_Global.param_changes   = False
    GUI_Global.param_available = True
    # mask flushing
    if len(GUI_Global.mask) > 0:
        dlg = wx.MessageDialog(GUI_Global.frame, TextPool.dialog_013,
            GUI_Global.VERSION, wx.YES_NO | wx.ICON_QUESTION)
        result = dlg.ShowModal()             
        dlg.Destroy()                         
        if result == wx.ID_NO:
            del GUI_Global.mask
            GUI_Global.mask = []
        else:
            # flush only the pixel masks
            for ii in reversed(range(len(GUI_Global.mask))):
                if (GUI_Global.mask[ii].validxy == True) :
                    del GUI_Global.mask[ii]
    # flushing fit and difference
    GUI_Global.corrected_available = False
    GUI_Global.fit_available = False
    # read a setup file line by line
    fd = open( filename, 'r' )
    array = []
    for line in fd:
        array.append( line )
    fd.close()

    GUI_Global.Current_Setup_File = filename

    text = ("Currently used parameter file:           %s" % 
        str(GUI_Global.Current_Setup_File))
    GUI_Global.current_par_file_t.SetLabel(text)

    molecs_name = None
    molecs_fit = None
    molecs_ab = None

    for ii in range(len(array)):
        #cut off (possible) newlines
        array[ii] = array[ii].strip('\n')
        # and if file was edited in MSDOS/Windows or MacOS also
        array[ii] = array[ii].strip('\r')  

        # interpretation / parser
        # Not used anymore
        # if strncmp('basedir:', array[ii], 8):
        #     old_basedir = GUI_Global.basedir_v
        #     GUI_Global.basedir_v = array[ii][8:]
        #     GUI_Global.basedir_v = GUI_Global.basedir_v.strip()
        #     if GUI_Global.basedir_v == '.':
        #         # older setup makes no sense in the way astronomers used the
        #         # code during the test workshop
        #         # GUI_Global.basedir_v = os.getcwd()
        #         GUI_Global.basedir_v = old_basedir 
        #     # check if that exists and contains a molecfit installation - 
        #     # if not it is a failure and old value is restored
        #     if GUI_Init.check_installation_basedir_executables() == False:
        #         dlg = wx.MessageDialog( GUI_Global.frame, 
		      #   'Parameter "basedir" does not point to a valid molecfit '\
        #         'installation\nOld value is restored !',
	       #      GUI_Global.VERSION, wx.OK)
        #         dlg.ShowModal() # Show it
        #         dlg.Destroy()   # finally destroy it when finished with OK
        #         GUI_Global.basedir_v = old_basedir  
        #         # rerun the check to resore also other executables 
        #         GUI_Init.check_installation_basedir_executables()         

        if strncmp('filename:', array[ii], 9):
            GUI_Global.filename_v = array[ii][9:]
            GUI_Global.filename_v = GUI_Global.filename_v.strip()
            # No
            # if GUI_Global.filename_v[:1] != '/':
            #     GUI_Global.filename_v = GUI_Global.basedir_v + '/' + \
            #         GUI_Global.filename_v 
        
        if strncmp('user_workdir:', array[ii], 13):
            GUI_Global.workdir_v = array[ii][13:]
            GUI_Global.workdir_v = GUI_Global.workdir_v.strip()
            if GUI_Global.workdir_v == '.':
                GUI_Global.workdir_v = os.getcwd()
        
                
        if strncmp('listname:', array[ii], 9):
            GUI_Global.listname_v = array[ii][9:]
            GUI_Global.listname_v = GUI_Global.listname_v.strip()

        if strncmp('trans:', array[ii], 6):
            GUI_Global.trans_v = int(array[ii][6:])

        if strncmp('columns:', array[ii], 8):
            text = array[ii][8:]
            text = text.strip()
            words = text.split()
            # print words
            if len(words) < 4 :
                dlg = wx.MessageDialog( GUI_Global.frame, 
		        'Parameter "columns:" setup file does not contain enough '+
                'entries\n\n  It should contain at least 4 entries.\n'+
                'Correct that parameter before calling the fit button',
	            GUI_Global.VERSION, wx.OK)
                dlg.ShowModal() # Show it
                dlg.Destroy()   # finally destroy it when finished with OK
            if len(words) > 0:
                GUI_Global.col_lam_v = words[0]
            else:
                GUI_Global.col_lam_v = ''
            if len(words) > 1:
                GUI_Global.col_flux_v = words[1]
            else:
                GUI_Global.col_flux_v = ''
            if len(words) > 2:
                GUI_Global.col_dflux_v = words[2]
            else:
                GUI_Global.col_dflux_v = ''
            if len(words) > 3:
                GUI_Global.col_mask_v = words[3]
            else:
                GUI_Global.col_mask_v = ''
            del words

        if strncmp('default_error:', array[ii], 14):
            GUI_Global.default_error_v = float(array[ii][14:])

        if strncmp('wlgtomicron:', array[ii], 12):
            GUI_Global.wlgtomicron_v = float(array[ii][12:])

        if strncmp('vac_air:', array[ii], 8):
            GUI_Global.vac_air_v = array[ii][8:]
            GUI_Global.vac_air_v = GUI_Global.vac_air_v.strip()

        if strncmp('wrange_include:', array[ii], 15):
            GUI_Global.wrange_include_v = array[ii][15:]
            GUI_Global.wrange_include_v = GUI_Global.wrange_include_v.strip()
            # No
            # if GUI_Global.wrange_include_v[:1] != '/' and \
            #    GUI_Global.wrange_include_v[:4] != 'none':
            #     GUI_Global.wrange_include_v = GUI_Global.basedir_v + '/' + \
            #         GUI_Global.wrange_include_v 

        if strncmp('wrange_exclude:', array[ii], 15):
            GUI_Global.wrange_exclude_v = array[ii][15:]
            GUI_Global.wrange_exclude_v = GUI_Global.wrange_exclude_v.strip()
            # No
            # if GUI_Global.wrange_exclude_v[:1] != '/'and \
            #    GUI_Global.wrange_exclude_v[:4] != 'none':
            #     GUI_Global.wrange_exclude_v = GUI_Global.basedir_v + '/' + \
            #         GUI_Global.wrange_exclude_v 

        if strncmp('prange_exclude:', array[ii], 15):
            GUI_Global.prange_exclude_v = array[ii][15:]
            GUI_Global.prange_exclude_v = GUI_Global.prange_exclude_v.strip()
            # No
            # if GUI_Global.prange_exclude_v[:1] != '/'and \
            #    GUI_Global.prange_exclude_v[:4] != 'none':
            #     GUI_Global.prange_exclude_v = GUI_Global.basedir_v + '/' + \
            #         GUI_Global.prange_exclude_v 

        if strncmp('output_dir:', array[ii], 11):
            GUI_Global.output_dir_v = array[ii][11:]
            GUI_Global.output_dir_v = GUI_Global.output_dir_v.strip()
        if strncmp('output_name:', array[ii], 12):
            GUI_Global.output_name_v = array[ii][12:]
            GUI_Global.output_name_v = GUI_Global.output_name_v.strip()
        
        # plot_creation and plot_range is not supported via GUI
        GUI_Global.plot_range    = 0
        GUI_Global.plot_creation = ''

        if strncmp('ftol:', array[ii], 5):
            GUI_Global.ftol_v = float(array[ii][5:])

        if strncmp('xtol:', array[ii], 5):
            GUI_Global.xtol_v = float(array[ii][5:])

        if strncmp('list_molec:', array[ii], 11):
            text  = array[ii][11:]
            text  = text.strip()
            text  = text.upper()
            molecs_name = text.split()

        if strncmp('fit_molec:', array[ii], 10):
            text  = array[ii][10:]
            text  = text.strip()
            molecs_fit = []
            for w in text.split():
                try:
                    molecs_fit.append(int(w))
                except ValueError:
                    raise ValueError("Could not convert %s word of fit_molec "
                                     "parameter to 0 or 1" % w)

        if strncmp('relcol:', array[ii], 7):
            text  = array[ii][7:]
            text  = text.strip()
            molecs_ab = []
            for w in text.split():
                try:
                    molecs_ab.append(float(w))
                except ValueError:
                    raise ValueError("Could not convert %s word of fit_molec "
                                     "parameter to float" % w)

        if strncmp('flux_unit:', array[ii], 10):
            GUI_Global.flux_unit_v = int(array[ii][10:])

        if strncmp('fit_back:', array[ii], 9):
            GUI_Global.popup_back.fit_back_v = int(array[ii][9:])

        if strncmp('telback:', array[ii], 8):
            GUI_Global.popup_back.telback_v = float(array[ii][8:])

        if strncmp('fit_cont:', array[ii], 9):
            GUI_Global.popup_continuum.fit_cont_v = int(array[ii][9:])

        if strncmp('cont_n:', array[ii], 7):
            GUI_Global.popup_continuum.cont_n_v = int(array[ii][7:])

        if strncmp('cont_const:', array[ii], 11):
            GUI_Global.popup_continuum.cont_const_v = float(array[ii][11:])

        if strncmp('fit_wlc:', array[ii], 8):
            GUI_Global.popup_wavelength.fit_wlc_v = int(array[ii][8:])

        if strncmp('wlc_n:', array[ii], 6):
            GUI_Global.popup_wavelength.wlc_n_v = int(array[ii][6:])

        if strncmp('wlc_const:', array[ii], 10):
            GUI_Global.popup_wavelength.wlc_const_v = float(array[ii][10:])

        if array[ii].startswith('slitw:'):
            GUI_Global.slitw_v = float(array[ii].split(':', 1)[-1].strip())

        if array[ii].startswith('slitw_key:'):
            GUI_Global.slitw_key_v = array[ii].split(':', 1)[-1].strip()

        if array[ii].startswith('pixsc:'):
            GUI_Global.pixsc_v = float(array[ii][6:])

        if array[ii].startswith('pixsc_key:'):
            GUI_Global.pixsc_key_v = array[ii].split(':', 1)[-1].strip()

        if array[ii].startswith('ref_atm:'):
            GUI_Global.ref_atm_v = array[ii].split(':', 1)[-1].strip()

        if array[ii].startswith('gdas_dir:'):
            GUI_Global.gdas_dir_v = array[ii].split(':', 1)[-1].strip()

        if array[ii].startswith('gdas_prof:'):
            GUI_Global.gdas_prof_v = array[ii].split(':', 1)[-1].strip()

        if array[ii].startswith('layers:'):
            GUI_Global.layers_v = int(array[ii].split(':', 1)[-1].strip())

        if array[ii].startswith('emix:'):
            GUI_Global.emix_v = float(array[ii].split(':', 1)[-1].strip())

        if strncmp('kernmode:', array[ii], 9):
            GUI_Global.popup_kernel.kernmode_v = int(array[ii][9:])

        if strncmp('fit_res_gauss:', array[ii], 14):
            GUI_Global.popup_kernel.fit_res_gauss_v = int(array[ii][14:])

        if strncmp('relres_box:', array[ii], 11):
            GUI_Global.popup_kernel.relres_box_v = float(array[ii][11:])

        if strncmp('fit_res_box:', array[ii], 12):
            GUI_Global.popup_kernel.fit_res_box_v = int(array[ii][12:])

        if strncmp('res_gauss:', array[ii], 10):
            GUI_Global.popup_kernel.res_gauss_v = float(array[ii][10:])

        if strncmp('fit_res_lorentz:', array[ii], 16):
            GUI_Global.popup_kernel.fit_res_lorentz_v = int(array[ii][16:])

        if strncmp('res_lorentz:', array[ii], 12):
            GUI_Global.popup_kernel.res_lorentz_v = float(array[ii][12:])

        if strncmp('kernfac:', array[ii], 8):
            GUI_Global.popup_kernel.kernfac_v = float(array[ii][8:])

        if strncmp('varkern:', array[ii], 8):
            GUI_Global.popup_kernel.varkern_v = int(array[ii][8:])

        if strncmp('kernel_file:', array[ii], 12):
            GUI_Global.popup_kernel.kernel_file_v = str(array[ii][12:]).strip()

        # hack for files NOT from ESO pipelines
        ambient_keys = ['obsdate', 'utc', 'telalt',
                        'rhum', 'pres', 'temp', 'm1temp',
                        'geoelev', 'longitude', 'latitude']
        for ak in ambient_keys:
            if array[ii].startswith('%s:' % ak):
                v = float(array[ii].split(':', 1)[-1].strip())
                setattr(GUI_Global, '%s_v' % ak, v)
            if array[ii].startswith('%s_key:' % ak):
                v = array[ii].split(':', 1)[-1].strip()
                setattr(GUI_Global, '%s_key_v' % ak, v)
        # end of hack for NON ESO pipeline products

        if strncmp('pwv:', array[ii], 4):
            #print 'found PWV'
            GUI_Global.pwv_par.pwv = float(array[ii][4:])


    del array

    if molecs_name:
        if molecs_fit is None and molecs_ab is None:
            raise ValueError("list_molec specified without relcol or "
                             "fit_molec")
        # set not specified to None
        if molecs_fit is None:
            molecs_fit = [None] * len(molecs_name)
        if molecs_ab is None:
            molecs_ab = [None] * len(molecs_name)

        if len(molecs_name) != len(molecs_fit):
            raise ValueError("length of list_molec and fit_molec do not match")
        if len(molecs_name) != len(molecs_ab):
            raise ValueError("length of list_molec and relcol do not match")
        # disable all molecuts
        selected = dict()
        for i in range(len(molecs_name)):
            selected[molecs_name[i]] = molecs_fit[i], molecs_ab[i]

        for m in GUI_Global.mol:
            if m.name not in selected:
                m.calc = False
                m.fit  = False
                m.ab   = 1.0
                continue
            m.calc = True
            fit, ab = selected[m.name]
            # override defaults if specified
            if fit is not None:
                m.fit = bool(fit)
            if ab is not None:
                m.ab = ab

    # Stupid... not created yet!
    # if os.access(GUI_Data.output_dir_path(),os.W_OK) \
    #     == False:
    #     text = 'You have no write permissions at output directory \n%s ' % \
    #         (GUI_Data.output_dir_path())
    #     dlg = wx.MessageDialog( GUI_Global.frame, text,
	   #  GUI_Global.VERSION, wx.OK)
    #     dlg.ShowModal() # Show it
    #     dlg.Destroy()   # finally destroy it when finished with OK
    #     return False

    # Create output directory if needed
    if os.path.exists(GUI_Data.output_dir_path()) == False:
        os.system('mkdir -p ' + GUI_Data.output_dir_path())

    # patch basedir if '.', output dir and plot_creation to none
    # No
    # patch_basedir(filename)  


    # prepeare the data file given in that very setup file for a plot
    if os.path.exists(GUI_Global.filename_v):
        
        #cmd = GUI_Global.prepgui_exe + ' ' + filename #+ ' 2>&1 > /dev/null'
        retval = GUI_Data.prepeare_fits_file(GUI_Global.filename_v)
        if retval == 0:
            #print "plot the data file"   
            prepearefile = \
                GUI_Data.output_dir_path() + '/' + \
                GUI_Global.output_name_v + '_gui.fits'
            GUI_Data.input_data(prepearefile,False,False,False)
            GUI_Data.read_input_file_pixel_exclusion_mask(prepearefile)
            GUI_Plot.draw_figure(False)
            #if os.path.exists(GUI_Global.wrange_include_v) == True:
            GUI_Data.import_mask_files(
                GUI_Global.wrange_include_v, True, True)
            #if os.path.exists(GUI_Global.wrange_exclude_v) == True:
            GUI_Data.import_mask_files(
                GUI_Global.wrange_exclude_v, False, True)
            #if os.path.exists(GUI_Global.prange_exclude_v) == True:
            GUI_Data.import_mask_files(
                GUI_Global.prange_exclude_v, False, False)
        else:
            text = 'The command prepguitable failed '\
              '(see molecfit documentation)\n'\
              'Use CONSOLE window to read molecfit/prepguitable result' #\
                #% GUI_Global.filename_v
            dlg = wx.MessageDialog( GUI_Global.frame, text,
	        GUI_Global.VERSION, wx.OK)
            dlg.ShowModal() # Show it
            dlg.Destroy()   # finally destroy it when finished with OK

    else:
        dlg = wx.MessageDialog( GUI_Global.frame, 
		      'Parameter "filename:" does not point to a valid FITS file ',
	    GUI_Global.VERSION, wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy()   # finally destroy it when finished with OK

    refresh_all_popups_on_parameter_change()  
    GUI_Global.param_changes = False
    return True

def set_default_parameters():
    GUI_Global.param_changes   = True
    GUI_Global.param_available = True

    # working directory - default generated via 'cwd' of code start
    ### GUI_Global.workdir_v = '.'      
    # base directory - default generated via the UNIX 'which' command on the
    # molecfit executable
    ### GUI_Global.basedir_v = '.'                                
    # name of parameter file  
    if hasattr(GUI_Global,"filename_v") == False:
        GUI_Global.filename_v = '** not yet given **'       
    GUI_Global.listname_v = 'none'       # list of files to be corrected
    GUI_Global.trans_v = 1               # type of filename: 
                                         #  1_v = transmission, 0 = radiation 

    # fit regions for the chi2 calculation (inclusion - is in wavelength space)
    GUI_Global.wrange_include_v = 'none' 
    # exclusion masks motivated by pyhsics (is in wavelength space)
    GUI_Global.wrange_exclude_v = 'none' 
    # exclusion masks motivated by instrument (is in pixel space)
    GUI_Global.prange_exclude_v = 'none'

    # output directory for the result files
    GUI_Global.output_dir_v = 'output'   
    # Name space for output files
    if hasattr(GUI_Global,"output_name_v") == False:
        GUI_Global.output_name_v = 'default_result'  
    # plotting options and individual plots  of the command line version of
    # molecfit are switched off in the GUI as the GUI does its own matplot-lib
    # graphs             
    GUI_Global.plot_creation_v = 'none'                          
    GUI_Global.plot_range_v = 0                                   

    # telescope backgound parameters
    GUI_Global.popup_back.fit_back_v = 0
    GUI_Global.popup_back.telback_v =  0.1
    # spectral continuum parameters
    GUI_Global.popup_continuum.fit_cont_v = 1
    GUI_Global.popup_continuum.cont_n_v = 3
    GUI_Global.popup_continuum.cont_const_v = 1.0
    # wavelength correction
    GUI_Global.popup_wavelength.fit_wlc_v = 1
    GUI_Global.popup_wavelength.wlc_n_v = 0
    GUI_Global.popup_wavelength.wlc_const_v = 0.0
    # line spread function and its parameters for the fit
    GUI_Global.popup_kernel.fit_res_box_v = 0
    GUI_Global.popup_kernel.relres_box_v = 0.0
    GUI_Global.popup_kernel.kernmode_v = 1
    GUI_Global.popup_kernel.fit_res_gauss_v = 1
    GUI_Global.popup_kernel.res_gauss_v = 1.0
    GUI_Global.popup_kernel.fit_res_lorentz_v = 0
    GUI_Global.popup_kernel.res_lorentz_v = 0.5
    GUI_Global.popup_kernel.kernfac_v = 30.0
    GUI_Global.popup_kernel.varkern_v = 1
    GUI_Global.popup_kernel.kernel_file_v = 'none'

    # NO real meaningfull DEFAULTS for these values
    GUI_Global.ftol_v  = 1e-2         # chi2 convergence criterion
    GUI_Global.xtol_v  = 1e-2         # parameter convergence criterion
    GUI_Global.pixsc_v = 0.086
    if hasattr(GUI_Global,"col_lam_v") == False:
        GUI_Global.col_lam_v = '** not yet given **'   
    if hasattr(GUI_Global,"col_flux_v") == False:
        GUI_Global.col_flux_v = '** not yet given **'
    if hasattr(GUI_Global,"col_dflux_v") == False:
        GUI_Global.col_dflux_v = 'NULL'
    if hasattr(GUI_Global,"col_mask_v") == False:
        GUI_Global.col_mask_v = 'NULL'
    GUI_Global.default_error_v = 0.01 # default error
    GUI_Global.wlgtomicron_v = 1e-3   # wavelength --> micron conversion factor
    GUI_Global.vac_air_v = 'air'      # air / vacuum line wavelengthes
    GUI_Global.flux_unit_v = 0
    # change related to ticket PIPE-4592 - although it is 
    # doubtfull that it is good to 'hardcode' it here (see text 
    # in JIRA
    # this keep this mark to make it removable later if required
    for jj in range(len(GUI_Global.mol)):
        #print 'name: '+GUI_Global.mol[jj].name
        if GUI_Global.mol[jj].name == 'CO2':
            #print 'CO2'
            GUI_Global.mol[jj].ab   = 1.06
        if GUI_Global.mol[jj].name == 'O3':
            #print 'O3'
            GUI_Global.mol[jj].ab   = 1.08            
    # end of PIPE-4592 changes

    refresh_all_popups_on_parameter_change()
    # list of masks used later in masking - init section 
    GUI_Global.mask = []
    GUI_Global.pwv_par.pwv = float(-1.0)

def create_and_fill_default_molecules():
    GUI_Global.mol = []
    # Main Molecules
    # H2O and CO2 are default molecules
    GUI_Global.mol.append(GUI_Global.Molec("H2O", 1.00,1,1))
    GUI_Global.mol.append(GUI_Global.Molec("CO2", 1.06,1,1))
    GUI_Global.mol.append(GUI_Global.Molec("O3",  1.08,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("N2O", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("CO",  1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("CH4", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("O2",  1.00,0,0))

    # Other Molecules
    GUI_Global.mol.append(GUI_Global.Molec("NO",  1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("SO2", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("NO2", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("NH3", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("HNO3",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("ClO", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("OCS", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("HOCl",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("N2",  1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("HCN", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("H2O2",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("C2H2",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("C2H6",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("COF2",1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("SF6", 1.00,0,0))
    GUI_Global.mol.append(GUI_Global.Molec("ClONO2", 1.00,0,0))

def write_paramfile(parent,filename,temporary):
    #if ((GUI_Global.param_available == False) or 
    if (
        (GUI_Global.Current_Setup_File[0:2] == '***')):
        dlg  = wx.MessageDialog( GUI_Global.frame, TextPool.menu_save_error,
           'ERROR', wx.OK | wx.ICON_HAND)
        dlg.ShowModal() 
        dlg.Destroy()   
        return False

    if temporary == False:
        GUI_Global.param_changes = False
        if os.path.exists(filename) == True:
            os.rename(filename,filename+'.bak') 
    filen = open(filename, 'w')
    filen.write('### Driver for MOLECFIT\n')
    filen.write('\n')
    # No
    # filen.write('## DIRECTORY STRUCTURE\n')
    # filen.write('\n')
    # filen.write('# Base directory (default: ".") for the following folder '\
    #             'structure:\n')
    # filen.write('#\n')
    # filen.write('#               |--bin/\n')
    # filen.write('#               |\n')
    # filen.write('#               |--config/\n')
    # filen.write('#  <basedir>----|\n')
    # filen.write('#               |--data/\n')
    # filen.write('#               |\n')
    # filen.write('#               |--output/\n')
    # filen.write('#\n')
    # filen.write('# A relative or an absolute path can be provided. In the '\
    #             'former case MOLECFIT\n')
    # filen.write('# has to be started in <basedir>.\n')
    # filen.writelines('basedir: '+ GUI_Global.basedir_v + '\n')
    # filen.write('\n')
    filen.write('# user working directory only important for REFLEX '\
                'workflow and GUI\n')
    filen.write('# not used by molecfit itself.\n')
    filen.writelines('user_workdir: '+ GUI_Global.workdir_v + '\n')
    filen.write('\n')
    filen.write('## INPUT DATA\n')
    filen.write('\n')
    filen.write('# Data file name (path relative to the current directory or absolute '\
                'path)\n')
    if hasattr(GUI_Global,'filename_v'):
        filen.writelines('filename: ' + GUI_Global.filename_v + '\n')
    else:
        filen.writelines('filename: none\n')
    filen.write('\n')
    filen.write('# ASCII list of files to be corrected for telluric absorption'\
                ' using the\n')
    filen.write('# transmission curve derived from the input reference file '\
                '(path of list and\n')
    filen.write('# listed files relative to the current directory or absolute path; '\
                'default: "none")\n')
    if hasattr(GUI_Global, 'listname_v'):
        filen.writelines('listname: ' + GUI_Global.listname_v + '\n')
    else:
        filen.writelines('listname: none \n')
    filen.write('\n')
    filen.write('# Type of input spectrum -- 1 = transmission (default); 0 = '\
        'emission\n')
    if hasattr(GUI_Global, 'trans_v'):
        filen.writelines('trans: %d\n' % GUI_Global.trans_v)
    else:
        filen.writelines('trans \n')
    filen.write('\n')
    filen.write('# Names of the file columns (table) or extensions (image) '\
        'containing:\n')
    filen.write('# Wavelength  Flux  Flux_Err  Mask\n')
    filen.write('# - Flux_Err and/or Mask can be avoided by writing \'NULL\'\n')
    filen.write('# - \'NULL\' is required for Wavelength if it is given by '\
        'header keywords\n')
    filen.write('# - parameter list: col_lam, col_flux, col_dflux, '\
        'and col_mask\n')
    if hasattr(GUI_Global, 'col_lam_v') and  hasattr(GUI_Global, 'col_flux_v') \
        and\
       hasattr(GUI_Global, 'col_flux_v') and hasattr(GUI_Global, 'col_dflux_v'):
       filen.writelines('columns: ' + GUI_Global.col_lam_v + ' ' + GUI_Global.col_flux_v + \
        ' ' + GUI_Global.col_dflux_v + ' ' + GUI_Global.col_mask_v + '\n')
    else:
        filen.writelines('columns: NULL NULL NULL NULL\n')
    filen.write('\n')
    filen.write('# Default error relative to mean for the case that the error '\
        'column is missing\n')
    if hasattr(GUI_Global, 'default_error_v'):
        filen.writelines('default_error: %r\n' % GUI_Global.default_error_v)
    else:
        filen.writelines('default_error \n')
    filen.write('\n')
    filen.write('# Multiplicative factor to convert wavelength to micron\n')
    filen.write('# (e.g. nm -> wlgtomicron = 1e-3)\n')
    if hasattr(GUI_Global, 'wlgtomicron_v'):
        filen.writelines('wlgtomicron: %r\n' % GUI_Global.wlgtomicron_v)
    else:
        filen.writelines('wlgtomicron: 1e-3\n')
    filen.write('\n')
    filen.write('# Wavelengths in vacuum (= vac) or air (= air)\n')
    if hasattr(GUI_Global, 'vac_air_v'):
        filen.writelines('vac_air: ' + GUI_Global.vac_air_v + '\n')
    else:
        filen.writelines('vac_air \n')
    filen.write('\n')
    filen.write('# ASCII or FITS table for wavelength ranges in micron to be'\
        ' fitted\n')
    filen.write('# (path relative to the current directory or absolute path; default: '\
        '"none")\n')
    if hasattr(GUI_Global, 'wrange_include_v'):
        if (GUI_Global.wrange_include_v[:1] != '/' ):
            if (GUI_Global.wrange_include_v[:4] != 'none'):
                GUI_Global.wrange_include_v = (
                GUI_Global.workdir_v + '/' + GUI_Global.wrange_include_v)
        filen.writelines('wrange_include: ' + GUI_Global.wrange_include_v +'\n')
    else:
        filen.writelines('wrange_include: none\n')
    filen.write('\n')
    filen.write('# ASCII or FITS table for wavelength ranges in micron to '\
        'be excluded from the\n')
    filen.write('# fit (path relative to the current directory or absolute path; default: '\
        '"none")\n')
    if hasattr(GUI_Global, 'wrange_exclude_v'):
        if (GUI_Global.wrange_exclude_v[:1] != '/' ):
            if (GUI_Global.wrange_exclude_v[:4] != 'none'):
                GUI_Global.wrange_exclude_v = (
                GUI_Global.workdir_v + '/' + GUI_Global.wrange_exclude_v)
        filen.write('wrange_exclude: ' + GUI_Global.wrange_exclude_v + '\n')
    else:
        filen.write('wrange_exclude: none\n')
    filen.write('\n')
    filen.write('# ASCII or FITS table for pixel ranges to be excluded '\
        'from the fit\n')
    filen.write('# (path relative to the current directory or absolute path; default: '\
        '"none")\n')
    if hasattr(GUI_Global, 'prange_exclude_v'):
        if (GUI_Global.prange_exclude_v[:1] != '/' ):
            if (GUI_Global.prange_exclude_v[:4] != 'none'):
                GUI_Global.prange_exclude_v = (
                GUI_Global.workdir_v + '/' + GUI_Global.prange_exclude_v)
        filen.write('prange_exclude: ' + GUI_Global.prange_exclude_v + '\n')
    else:
        filen.write('prange_exclude: none\n')
    filen.write('\n')
    filen.write('## RESULTS\n')
    filen.write('\n')
    filen.write('# Directory for output files (path relative to the current directory or '\
        'absolute path)\n')
    if temporary == False and hasattr(GUI_Global, 'output_dir_v'):
        GUI_Data.check_output()
    if hasattr(GUI_Global, 'output_dir_v'):
        filen.write('output_dir: ' + GUI_Data.output_dir_path() + '\n')
    else:
        filen.write('output_dir: output\n')
    filen.write('\n')
    filen.write('# Name for output files\n')
    filen.write('# (supplemented by "_fit" or "_tac" as well as ".asc", ".atm", ".fits",\n')
    filen.write('# ".par, ".ps", and ".res")\n')
    if hasattr(GUI_Global, 'output_name_v'):
        filen.write('output_name: ' + GUI_Global.output_name_v + '\n')
    else:
        filen.write('output_name: output\n')
    filen.write('\n')
    filen.write('# Plot creation: gnuplot is used to create control plots\n')
    filen.write('# W - screen output only (incorporating wxt terminal in '\
        'gnuplot)\n')
    filen.write('# X - screen output only (incorporating x11 terminal in '\
        'gnuplot)\n')
    filen.write('# P - postscript file labelled \'<output_name>.ps\', '\
        'stored in <output_dir>\n')
    filen.write('# combinations possible, i.e. WP, WX, XP, WXP (however, '\
        'keep the order!)\n')
    filen.write('# all other input: no plot creation is performed\n')
    filen.write('plot_creation: none' + '\n')
    filen.write('\n')
    filen.write('# Create plots for individual fit ranges? -- 1 = yes; '\
        '0 = no\n')
    filen.write('plot_range: ' + '0' + '\n')
    filen.write('\n')
    filen.write('## FIT PRECISION\n')
    filen.write('\n')
    filen.write('# Relative chi2 convergence criterion\n')
    if hasattr(GUI_Global, 'ftol_v'):
        filen.write('ftol: %r\n' % GUI_Global.ftol_v)
    else:
        filen.write('ftol\n')
    filen.write('\n')
    filen.write('# Relative parameter convergence criterion\n')
    if hasattr(GUI_Global, 'xtol_v'):
        filen.write('xtol: %r\n' % GUI_Global.xtol_v)
    else:
        filen.write('xtol\n')
    filen.write('\n')
    filen.write('## MOLECULAR COLUMNS\n')
    filen.write('\n')

    # convert molecule CLASS to individual parameter(s)
    list_molec_v = ' '.join('%s' % m.name for m in GUI_Global.mol if m.calc)
    fit_molec_v  = ' '.join('%d' % m.fit for m in GUI_Global.mol if m.calc)
    relcol_v     = ' '.join('%r' % m.ab for m in GUI_Global.mol if m.calc)

    filen.write('# List of molecules to be included in the model\n')
    filen.write('# (default: \'H2O\', N_val: nmolec)\n')
    filen.write('list_molec: ' + list_molec_v + '\n')
    filen.write('\n')
    filen.write('# Fit flags for molecules -- 1 = yes; 0 = no (N_val: '\
        'nmolec)\n')
    filen.write('fit_molec: ' + fit_molec_v + '\n')
    filen.write('\n')
    filen.write('# Values of molecular columns, expressed relatively to '\
        'the input ATM profile\n')
    filen.write('# columns (N_val: nmolec) [1 = 100%]\n')
    filen.write('relcol: ' + relcol_v + '\n')

    filen.write('\n')
    filen.write('## BACKGROUND AND CONTINUUM\n')
    filen.write('\n')
    filen.write('# Conversion of fluxes from phot/(s*m2*mum*as2) (emission '\
        'spectrum only) to\n')
    filen.write('# flux unit of observed spectrum:\n')
    filen.write('# 0: phot/(s*m^2*mum*as^2) [no conversion]\n')
    filen.write('# 1: W/(m^2*mum*as^2)\n')
    filen.write('# 2: erg/(s*cm^2*A*as^2)\n')
    filen.write('# 3: mJy/as^2\n')
    filen.write('# For other units the conversion factor has to be '\
        'considered as constant term\n')
    filen.write('# of the continuum fit.\n')
    filen.write('flux_unit: ' + '%i' % GUI_Global.flux_unit_v + '\n')
    filen.write('\n')
    filen.write('# Fit of telescope background -- 1 = yes; 0 = no (emission '\
        'spectrum only)\n')
    if hasattr(GUI_Global.popup_back, 'fit_back_v'):
        filen.write('fit_back: ' + '%i' % GUI_Global.popup_back.fit_back_v+'\n')
    else:
        filen.write('fit_back\n')
    filen.write('\n')
    filen.write('# Initial value for telescope background fit (range: [0,1])\n')
    if hasattr(GUI_Global.popup_back, 'telback_v'):
        filen.write('telback: %r\n' % GUI_Global.popup_back.telback_v)
    else:
        filen.write('telback\n')
    filen.write('\n')
    filen.write('# Polynomial fit of continuum --> degree: cont_n\n')
    if hasattr(GUI_Global.popup_continuum, 'fit_cont_v'):
        filen.write('fit_cont: ' + '%i' %GUI_Global.popup_continuum.fit_cont_v \
            + '\n')
    else:
        filen.write('fit_cont\n')
    filen.write('\n')
    filen.write('# Degree of coefficients for continuum fit\n')
    if hasattr(GUI_Global.popup_continuum, 'cont_n_v'):
        filen.write('cont_n: ' + '%i' % GUI_Global.popup_continuum.cont_n_v \
            + '\n')
    else:
        filen.write('cont_n\n')
    filen.write('\n')
    filen.write('# Initial constant term for continuum fit (valid for all '\
        'fit ranges)\n')
    filen.write('# (emission spectrum: about 1 for correct flux_unit)\n')

    if hasattr(GUI_Global.popup_continuum, 'cont_const_v'):
        filen.write('cont_const: %r\n' %
                    GUI_Global.popup_continuum.cont_const_v)
    else:
        filen.write('cont_const\n')
    filen.write('\n')
    filen.write('## WAVELENGTH SOLUTION\n')
    filen.write('\n')
    filen.write('# Refinement of wavelength solution using a polynomial of '\
        'degree wlc_n\n')
    if hasattr(GUI_Global.popup_wavelength, 'fit_wlc_v'):
        filen.write('fit_wlc: ' + '%i' % \
            GUI_Global.popup_wavelength.fit_wlc_v + '\n')
    else:
        filen.write('fit_wlc\n')
    filen.write('\n')
    filen.write('# Polynomial degree of the refined wavelength solution\n')
    if hasattr(GUI_Global.popup_wavelength, 'wlc_n_v'):
        filen.write('wlc_n: ' + '%i' %GUI_Global.popup_wavelength.wlc_n_v +'\n')
    else:
        filen.write('wlc_n\n')
    filen.write('\n')
    filen.write('# Initial constant term for wavelength correction (shift '\
        'relative to half\n')
    filen.write('# wavelength range)\n')
    if hasattr(GUI_Global.popup_wavelength, 'wlc_const_v'):
        filen.write('wlc_const: %r\n' %
            GUI_Global.popup_wavelength.wlc_const_v)
    else:
        filen.write('wlc_const\n')
    filen.write('\n')
    filen.write('## RESOLUTION\n')
    filen.write('\n')
    filen.write('# Fit resolution by boxcar -- 1 = yes; 0 = no\n')
    if hasattr(GUI_Global.popup_kernel, 'fit_res_box_v'):
        filen.write('fit_res_box: ' + '%i' % \
            GUI_Global.popup_kernel.fit_res_box_v + '\n')
    else:
        filen.write('fit_res_box\n')
    filen.write('\n')
    filen.write('# Initial value for FWHM of boxcar relative to slit '\
        'width (>= 0. and <= 2.)\n')
    if hasattr(GUI_Global.popup_kernel, 'relres_box_v'):
        filen.write('relres_box: %r\n' % GUI_Global.popup_kernel.relres_box_v)
    else:
        filen.write('relres_box\n')
    filen.write('\n')
    filen.write('# Voigt profile approximation instead of independent '\
        'Gaussian and Lorentzian\n')
    filen.write('# kernels? -- 1 = yes; 0 = no\n')
    if hasattr(GUI_Global.popup_kernel, 'kernmode_v'):
        filen.write('kernmode: ' + '%i' %GUI_Global.popup_kernel.kernmode_v  \
            + '\n')
    else:
        filen.write('kernmode\n')
    filen.write('\n')
    filen.write('# Fit resolution by Gaussian -- 1 = yes; 0 = no\n')
    if hasattr(GUI_Global.popup_kernel, 'fit_res_gauss_v'):
        filen.write('fit_res_gauss: ' + '%i' % \
            GUI_Global.popup_kernel.fit_res_gauss_v  + '\n')
    else:
        filen.write('fit_res_gauss\n')
    filen.write('\n')
    filen.write('# Initial value for FWHM of Gaussian in pixels\n')
    if hasattr(GUI_Global.popup_kernel, 'res_gauss_v'):
        filen.write('res_gauss: %r\n' % GUI_Global.popup_kernel.res_gauss_v)
    else:
        filen.write('res_gauss\n')
    filen.write('\n')
    filen.write('# Fit resolution by Lorentzian -- 1 = yes; 0 = no\n')
    if hasattr(GUI_Global.popup_kernel, 'fit_res_lorentz_v'):
        filen.write('fit_res_lorentz: ' + '%i' % \
        GUI_Global.popup_kernel.fit_res_lorentz_v  + '\n')
    else:
        filen.write('fit_res_lorentz: 0\n')
    filen.write('\n')
    filen.write('# Initial value for FWHM of Lorentzian in pixels\n')
    if hasattr(GUI_Global.popup_kernel, 'res_lorentz_v'):
        filen.write('res_lorentz: %r\n' % GUI_Global.popup_kernel.res_lorentz_v)
    else:
        filen.write('res_lorentz\n')
    filen.write('\n')
    filen.write('# Size of Gaussian/Lorentzian/Voigtian kernel in FWHM\n')
    if hasattr(GUI_Global.popup_kernel, 'kernfac_v'):
        filen.write('kernfac: %r\n' % GUI_Global.popup_kernel.kernfac_v)
    else:
        filen.write('kernfac\n')
    filen.write('\n')
    filen.write('# Variable kernel (linear increase with wavelength)? -- '\
        '1 = yes; 0 = no\n')
    if hasattr(GUI_Global.popup_kernel, 'varkern_v'):
        filen.write('varkern: ' + '%i' %GUI_Global.popup_kernel.varkern_v +'\n')
    else:
        filen.write('varkern\n')
    filen.write('\n')
    filen.write('# ASCII file for kernel elements (one per line; '\
        'normalisation not required)\n')
    filen.write('# instead of synthetic kernel consisting of boxcar, '\
        'Gaussian, and Lorentzian\n')
    filen.write('# components (path relative to the current directory or absolute path; '\
        'default: "none")\n')
    if hasattr(GUI_Global.popup_kernel, 'kernel_file_v'):
        filen.write('kernel_file: ' + '%s' % \
            GUI_Global.popup_kernel.kernel_file_v + '\n')
    else:
        filen.write('kernel_file: none\n')
    filen.write('\n')

    # hack for files NOT from ESO pipelines
    filen.write("""\
## AMBIENT PARAMETERS

# If the input data file contains a suitable FITS header, the keyword names of
# the following parameters will be read, but the corresponding values will not
# be used. The reading of parameter values from this file can be forced by
# setting keywords to NONE.

""")
    # hack for files NOT from ESO pipelines
    ambient_keys = [
        ('obsdate', '# Observing date in years or MJD in days', (-1., 'MJD-OBS')),
        ('utc', '# UTC in s', (-1., 'UTC')),
        ('telalt', '# Telescope altitude angle in deg', (90., 'ESO TEL ALT')),
        ('rhum', '# Humidity in %', (15., 'ESO TEL AMBI RHUM')),
        ('pres', '# Pressure in hPa', (750., 'ESO TEL AMBI PRES START')),
        ('temp', '# Ambient temperature in deg C', (15., 'ESO TEL AMBI TEMP')),
        ('m1temp', '# Mirror temperature in deg C', (15., 'ESO TEL TH M1 TEMP')),
        ('geoelev', '# Elevation above sea level in m (default is Paranal: 2635m)',
         (2635., 'ESO TEL GEOELEV')),
        ('longitude', '# Longitude (default is Paranal: -70.4051)',
         (-70.4051, 'ESO TEL GEOLON')),
        ('latitude', '# Latitude (default is Paranal: -24.6276)',
         (-24.6276, 'ESO TEL GEOLAT'))
        ]
    for ak, c, d in ambient_keys:
        filen.write(c + '\n')
        vval = getattr(GUI_Global, '%s_v' % ak, None)
        vkey = getattr(GUI_Global, '%s_key_v' % ak, None)
        if vval is None and vkey is not None:
            filen.write('%s\n' % ak)
            filen.write('%s_key: %s\n\n' % (ak, vkey))
        elif vval is not None and vkey is None:
            filen.write('%s: %r\n' % (ak, vval))
            filen.write('%s_key\n\n' % ak)
        elif vval is None and vkey is None:
            filen.write('%s\n' % ak)
            filen.write('%s_key\n\n' % ak)
        else:
            filen.write('%s: %r\n' % (ak, vval))
            filen.write('%s_key: %s\n\n' % (ak, vkey))

    # end of hack for NON ESO pipeline products


    filen.write('## INSTRUMENTAL PARAMETERS\n')
    filen.write('\n')

    filen.write('# Slit width in arcsec (taken from FITS header if present)\n')
    v = getattr(GUI_Global, 'slitw_v', 0.4)
    filen.write('slitw: %r\n' % v)
    v = getattr(GUI_Global, 'slitw_key_v', 'ESO INS SLIT1 WID')
    filen.write('slitw_key: %s\n' % v)

    filen.write('\n')

    filen.write('# Pixel scale in arcsec (taken from this file only)\n')
    v = getattr(GUI_Global, 'pixsc_v', 0.089)
    filen.write('pixsc: %r\n' % v)
    v = getattr(GUI_Global, 'pixsc_key_v', 'NONE')
    filen.write('pixsc_key: %s\n' % v)

    filen.write('\n')

    filen.write('## ATMOSPHERIC PROFILES\n\n')
    filen.write('# Reference atmospheric profile\n')
    v = getattr(GUI_Global, 'ref_atm_v', 'equ.atm')
    filen.write('ref_atm: %s\n\n' % v)
    filen.write("""\
# Specific GDAS-like input profile (P[hPa] HGT[m] T[K] RELHUM[%]) (path
# relative to the installation directory or absolute path). In the case of "none", no GDAS
# profiles will be considered. The default "auto" performs an automatic
# retrieval.
""")
    v = getattr(GUI_Global, 'gdas_dir_v', 'data/profiles/grib')
    filen.write('gdas_dir: %s\n' % v)
    v = getattr(GUI_Global, 'gdas_prof_v', 'auto')
    filen.write('gdas_prof: %s\n\n' % v)
    filen.write("""\
# Grid of layer heights for merging ref_atm and GDAS profile. Fixed grid = 1
# (default) and natural grid = 0.
""")
    v = getattr(GUI_Global, 'layers_v', 1)
    filen.write('layers: %d\n\n' % v)
    filen.write("""\
# Upper mixing height in km (default: 5) for considering data of a local meteo
# station. If emix is below geoelev, rhum, pres, and temp are not used for
# modifying the corresponding profiles.
""")
    v = getattr(GUI_Global, 'emix_v', 5.)
    filen.write('emix: %r\n\n' % v)

    filen.write('# PWV value in mm for the input water vapour profile. '\
        'The merged profile\n')
    filen.write('# composed of ref_atm, GDAS, and local meteo data will '\
        'be scaled to this value\n')
    filen.write('# if pwv > 0 (default: -1 -> no scaling).\n')
    if hasattr(GUI_Global,'pwv_par'):
        filen.write('%r\n\n' % GUI_Global.pwv_par)

    # THE FOLLOWING IS A HIDDEN gui ONLY PARAMETER IMPROVING PLOTS
    filen.write('# internal GUI specific parameter\n')
    filen.write('clean_mflux: 1\n\n')

    filen.write('end\n')
    filen.close()
    if temporary == False:
        if GUI_Global.param_available == True:
            refresh_all_popups_on_parameter_change()
            GUI_Global.param_changes = False
    return True

def generate_and_write_mask_files():
    # check if mask files are required in name setup
    if hasattr(GUI_Global,'prange_exclude_v') == False:
        GUI_Global.prange_exclude_v = 'none'
    if (GUI_Global.prange_exclude_v[:4] == 'none'):
        #print "prange test"
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == True) and 
                (GUI_Global.mask[ii].include == False)):
                found = 1
        if found == 1 :
            #print "FOUND PIXEL MASK"
            dlg  = wx.MessageDialog( GUI_Global.frame, 
                "Filename for pixel exclusion mask not yet defined\n\n"\
                "goto \n   Settings->Working Directories and Files\n"\
                "and add filenames.", 
                "INFO", wx.OK)

            dlg.ShowModal() 
            dlg.Destroy()
            return False

    if hasattr(GUI_Global,'wrange_exclude_v') == False:
        GUI_Global.wrange_exclude_v = 'none'
    if (GUI_Global.wrange_exclude_v[:4] == 'none'):
        #print "wrange test"
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == False) and 
                (GUI_Global.mask[ii].include == False)):
                found = 1
        if found == 1 :
            #print "FOUND WAVELENGTH MASK"
            dlg  = wx.MessageDialog( GUI_Global.frame, 
                "Filename for wavelength exclusion mask not yet defined\n\n"\
                "goto \n   Settings->Working Directories and Files\n"\
                "and add filenames.",  
                "INFO", wx.OK)
            dlg.ShowModal() 
            dlg.Destroy()
            return False

    if hasattr(GUI_Global,'wrange_include_v') == False:
        GUI_Global.wrange_include_v = 'none'
    if (GUI_Global.wrange_include_v[:4] == 'none'):
        #print "wrange test"
        found = 0
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].validxy == False) and 
                (GUI_Global.mask[ii].include == True)):
                found = 1
        if found == 1 :
            #print "FOUND INCLUDE"
            dlg  = wx.MessageDialog( GUI_Global.frame, 
                "Filename for wavelength exclusion mask not yet defined\n\n"\
                "goto \n   Settings->Working Directories and Files\n"\
                "and add filenames.",  
                "INFO", wx.OK)
            dlg.ShowModal() 
            dlg.Destroy()
            return False

    if (GUI_Global.wrange_include_v[:4] != 'none'):
        if os.path.exists(os.path.dirname(GUI_Global.wrange_include_v)) == False:
            if len(os.path.dirname(GUI_Global.wrange_include_v)) > 1:
                os.system('mkdir -p '+os.path.dirname(GUI_Global.wrange_include_v))
        if os.path.exists(GUI_Global.wrange_include_v) == True:
            os.rename(GUI_Global.wrange_include_v,GUI_Global.wrange_include_v+'.bak') 
        filen = open(GUI_Global.wrange_include_v, 'w')
        written_mask = False
        for ii in range(len(GUI_Global.mask)):
            if GUI_Global.mask[ii].include == True:
                written_mask = True
                filen.write(
                    '%r %r\n' % (GUI_Global.mask[ii].startwl,
                                 GUI_Global.mask[ii].endwl))
        # to ensure an empty mask file if no mask exists
        if written_mask == False:
            filen.write('\n')
        filen.close()

    if (GUI_Global.wrange_exclude_v[:4] != 'none'):
        if os.path.exists(os.path.dirname(GUI_Global.wrange_exclude_v)) == False:
            if len(os.path.dirname(GUI_Global.wrange_exclude_v)) > 1:
                os.system('mkdir -p '+os.path.dirname(GUI_Global.wrange_exclude_v))
        if os.path.exists(GUI_Global.wrange_exclude_v) == True:
            os.rename(GUI_Global.wrange_exclude_v,GUI_Global.wrange_exclude_v+'.bak') 
        filen = open(GUI_Global.wrange_exclude_v, 'w')
        written_mask = False
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].include == False) and
                (GUI_Global.mask[ii].validxy == False)):
                written_mask = True
                filen.write(
                    '%r %r\n' % (GUI_Global.mask[ii].startwl,
                                 GUI_Global.mask[ii].endwl))
        # to ensure an empty mask file if no mask exists
        if written_mask == False:
            filen.write('\n')
        filen.close()

    if (GUI_Global.prange_exclude_v[:4] != 'none'):
        if os.path.exists(os.path.dirname(GUI_Global.prange_exclude_v)) == False:
            if len(os.path.dirname(GUI_Global.prange_exclude_v)) > 1:
                os.system('mkdir -p '+os.path.dirname(GUI_Global.prange_exclude_v))
        if os.path.exists(GUI_Global.prange_exclude_v) == True:
            os.rename(GUI_Global.prange_exclude_v,GUI_Global.prange_exclude_v+'.bak') 
        filen = open(GUI_Global.prange_exclude_v, 'w')
        written_mask = False
        for ii in range(len(GUI_Global.mask)):
            if ((GUI_Global.mask[ii].include == False) and
                (GUI_Global.mask[ii].validxy == True)):
                written_mask = True
                filen.write(
                    '%r %r\n' % (GUI_Global.mask[ii].startx+1,
                                 GUI_Global.mask[ii].endx+1))
        # to ensure an empty mask file if no mask exists
        if written_mask == False:
            filen.write('\n')
        filen.close()

    return True
            
