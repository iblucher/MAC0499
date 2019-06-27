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
import SM02GUI_Plot        as GUI_Plot
#import SM02GUI_Status_Text as GUI_Status_Text
import SM02GUI_Help        as GUI_Help
import SM02GUI_Popups      as GUI_Popup
import SM02GUI_Data        as GUI_Data
import SM02GUI_Help_Dialog_ToolTip_Texts as TextPool
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

PIX_EXCL_COLOR = "#fe9d00" # orange
WAV_EXCL_COLOR = "#ff0000" # red
WAV_INCL_COLOR = "#00ff00" # green
SELECT_COLOR   = "#0000ff" # blue

def EndBusyCursor():
    """ wxpython 3 asserts if its not set """
    try:
        wx.EndBusyCursor()
    except:
        pass

def MergeMasksOfSameType():
    # combine overlapping masks if they are of same type (pixel exclusion, 
    # wavelength inclusion and wavelength  exclusion).
    if len(GUI_Global.mask) < 2:
        # nothing to do
        return

    #
    # REMARK - all calculations are in PIXEL space (also for the wavelength 
    #          masks) as comparison of floats never ever will be equal
    #

    # iterate all pairs until no more masks can be merged
    prev = len(GUI_Global.mask) + 1
    while prev != len(GUI_Global.mask):
        prev = len(GUI_Global.mask)
        for ii, m1 in enumerate(GUI_Global.mask[:-1]):
            #calculate wl space if required
            if m1.validxy == True:
                m1.startwl = GUI_Global.X[m1.startx]
                m1.endwl = GUI_Global.X[m1.endx]
            # calculate pix space if required
            else :
                 m1.startx = GUI_Data.pixel_of_data(m1.startwl)
                 m1.endx = GUI_Data.pixel_of_data(m1.endwl)
            for m2 in GUI_Global.mask[ii + 1:]:
                # first check if of same type
                if ((m1.validxy == m2.validxy) and
                    (m1.include == m2.include)):
                    # calculate wl space if required
                    if m2.validxy == True:
                        m2.startwl = GUI_Global.X[m2.startx]
                        m2.endwl = GUI_Global.X[m2.endx]
                    # calculate pix space if required
                    else :
                        m2.startx = GUI_Data.pixel_of_data(m2.startwl)
                        m2.endx = GUI_Data.pixel_of_data(m2.endwl)
                    # check for overlap - if yes modify first one and set second
                    # one as inactive
                    if m1.startx <= m2.startx and \
                       m1.endx >= m2.endx:
                        # m2 completeley in m1
                        m2.active = False
                    elif m1.startx >= m2.startx and \
                         m1.endx <= m2.endx:
                        # m1 completeley in m2
                        m1.startwl = m2.startwl
                        m1.startx  = m2.startx
                        m1.endwl   = m2.endwl
                        m1.endx    = m2.endx
                        m2.active = False
                    elif m1.startx <= m2.startx and \
                         m1.endx   >  m2.startx:
                        # m2 overlaps on right of m1
                        m1.endwl = m2.endwl
                        m1.endx  = m2.endx
                        m2.active = False
                    elif m1.endx   >= m2.endx and \
                         m1.startx <= m2.endx:
                        # m2 overlaps on left of m1
                        m1.startwl = m2.startwl
                        m1.startx  = m2.startx
                        m2.active = False

        # delete inactive masks
        for ii in reversed(range(0,len(GUI_Global.mask))):
            if GUI_Global.mask[ii].active == False :
                GUI_Global.mask_changes = True
                del GUI_Global.mask[ii]
                if (ii == GUI_Global.popup_mask_edit.current):
                    GUI_Global.popup_mask_edit.current = -99
                    GUI_Global.popup_mask_edit.Show(False)
                    GUI_Global.popup_mask_edit.Showed = False
    # while end

    if GUI_Global.mask_changes == True:
        draw_figure(True)
            
def MostNearbyMask(x):
    # find the most nearby mask, on overlap return exact match
    if len(GUI_Global.mask) < 1:
        return -99
    active_masks = (m for m in GUI_Global.mask if m.active)
    diff = []
    for ii, m in enumerate(active_masks):
        diff.append((ii, min(abs(m.startwl - x), abs(m.endwl - x)), m))

    def in_mask(m):
        return (m.startwl < x) and (m.endwl > x)

    diff = sorted(diff, key=lambda x: x[1])

    # if two masks overlap return the one that was clicked
    if len(diff) > 1:
        closest, second = diff[:2]
        if in_mask(second[2]) and not in_mask(closest[2]):
            return second[0]
    # else return closest
    return diff[0][0]

def OnButton_mask_done(e):
    # feedback for DONE button in mask tool
    #GUI_Parameters.generate_and_write_mask_files()
    GUI_Popup.check_mask_names(False)
    GUI_Global.popup_mask.Show(False) 
    GUI_Global.popup_mask_edit.Show(False)
    GUI_Global.popup_mask_edit.Showed = False
    if GUI_Global.mask_changes:
        GUI_Global.param_changes = True
    EndBusyCursor()
    GUI_Plot.draw_figure(True)

def onclick_release(event):
    # release is catched only in EDIT mode when dragging 
    # in this mode release means a mask is created / ready 
    GUI_Global.edit_started = False
    GUI_Global.patch_started = False
    GUI_Global.maskpanel.canvas.mpl_disconnect(GUI_Global.maskpanel.cid_rel)

    # out of axis -> drop selection
    if event.xdata is None:
        EndBusyCursor()
        draw_figure(True)
        return

    x_start = min(GUI_Global.edit_start_x, event.xdata)
    width = abs(GUI_Global.edit_start_x - event.xdata)
    is_include = bool(GUI_Global.maskpanel.cb_edit_f.GetValue())
    m = GUI_Global.Mask(True, is_include,
                        GUI_Data.pixel_of_data(x_start),
                        GUI_Data.pixel_of_data(x_start + width),
                        False, x_start, x_start + width)
    GUI_Global.mask.append(m)
    GUI_Global.mask_changes = True
    MergeMasksOfSameType()
    draw_figure(True)

def onclick_graph(event):
    # click with the left mouse button starts a new region if we are in edit
    # mode
    # if we are in EDIT mode the click causes the start of a mask - but only
    # if the user is not in zoom or pan mode
    if (((GUI_Global.maskpanel.cb_edit_f.GetValue() == True) or
         (GUI_Global.maskpanel.cb_edit_m.GetValue() == True)) and
         (event.button == 1) and 
         (GUI_Global.maskpanel.axes.get_navigate_mode() == None)):
        if GUI_Global.edit_started == False:
            # activate a callback on mouse release
            GUI_Global.maskpanel.cid_rel = \
                GUI_Global.maskpanel.canvas.mpl_connect(
                'button_release_event' , onclick_release)    
            GUI_Global.edit_started = True 
            if (GUI_Global.popup_mask_edit.Showed == True):
                OnButton_mask_edit_done(event)
            GUI_Global.edit_start_x = event.xdata
            GUI_Global.edit_start_y = event.ydata
            #print "edit starts AT x",GUI_Global.edit_start_x
            #print "edit starts AT y",GUI_Global.edit_start_y
    else:
        if ((event.button == 3) and
            (GUI_Global.maskpanel.axes.get_navigate_mode() == None)):
            # if it is the button 3 
            # this shows the data of the most nearby mask (if any) 
            # (in edit popup)
            ii = MostNearbyMask(event.xdata)
            if ii >= 0:
                #print 'mask %i from %f to %f is most nearby' % (ii, 
                #    GUI_Global.mask[ii].startwl,GUI_Global.mask[ii].endwl)
                refresh_edit_text(ii)
                EndBusyCursor()
                GUI_Global.popup_mask_edit.Show(True)
                GUI_Global.popup_mask_edit.Showed = True
            else:
                GUI_Help.refresh_error(TextPool.error_001)

def enter_graph(event):
    # activate other events when entering the pure plot drawing region
    GUI_Global.maskpanel.cid_move = GUI_Global.maskpanel.canvas.mpl_connect(
        'motion_notify_event', move_graph)
    GUI_Global.maskpanel.cid_click = GUI_Global.maskpanel.canvas.mpl_connect(
        'button_press_event' , onclick_graph)   
    text = (('wavelength: %.6g     y: %.6g\n' +
             'pixel:      %.4i    data value: %g\n') % 
             (event.xdata,event.ydata,GUI_Data.pixel_of_data(event.xdata)+1,
              GUI_Global.Y[GUI_Data.pixel_of_data(event.xdata)]))
    GUI_Global.maskpanel.position_text.SetLabel(text)  
    # change cursor to cross if we are in EDIT mode
    if GUI_Global.maskpanel.cb_edit_f.GetValue():
        wx.BeginBusyCursor(wx.CROSS_CURSOR)
    if GUI_Global.maskpanel.cb_edit_m.GetValue():
        wx.BeginBusyCursor(wx.CROSS_CURSOR)

def exit_graph(event):
    # DE-activate events when leaving the pure drawing (but still being in the 
    # canvas) area - this is required to suppress matplotlib error messages
    GUI_Global.maskpanel.canvas.mpl_disconnect(GUI_Global.maskpanel.cid_move)
    GUI_Global.maskpanel.canvas.mpl_disconnect(GUI_Global.maskpanel.cid_click)
    GUI_Global.maskpanel.position_text.SetLabel("")    
    if GUI_Global.maskpanel.cb_edit_f.GetValue():
        EndBusyCursor()
    if GUI_Global.maskpanel.cb_edit_m.GetValue():
        EndBusyCursor()

def move_graph(event):
    # display the position of the mouse in wavelength coordinates to the screen
    # this if is required if user moves mouse very fast and thus the MOVE event
    # is catched BEFORE the EXIT event (which switches off the MOVE)
    if (event.xdata != None) :
        text = (('wavelength: %.6g     y: %.6g\n' +
             'pixel:      %.4i    data value: %g\n') % 
             (event.xdata,event.ydata,GUI_Data.pixel_of_data(event.xdata)+1,
              GUI_Global.Y[GUI_Data.pixel_of_data(event.xdata)]))
        GUI_Global.maskpanel.position_text.SetLabel(text) 
        if hasattr(GUI_Global,'edit_box'):
            del GUI_Global.edit_box 
            if hasattr(GUI_Global,'patch'):
                del GUI_Global.patch
        if GUI_Global.edit_started == True:
            if GUI_Global.patch_started == False:
                GUI_Global.patch_started = True
            else:
                del GUI_Global.maskpanel.axes.patches[-1]
            # modify the line/box drawn along the mask
            x_start = min(event.xdata, GUI_Global.edit_start_x)
            GUI_Global.edit_box = matplotlib.patches.Rectangle(
                    (x_start, GUI_Global.mask_min),
                    abs(GUI_Global.edit_start_x - event.xdata),
                    GUI_Global.mask_max,
                    color=SELECT_COLOR, alpha=0.1)
            GUI_Global.patch = \
                GUI_Global.maskpanel.axes.add_patch(GUI_Global.edit_box)
            # refresh all changes
            GUI_Global.maskpanel.canvas.draw()
       
def create_popup_mask(parent):

    # TODO: add option to cancel a mask change

    GUI_Global.edit_started  = False
    GUI_Global.patch_started = False;
    
    # Create the mask tool
    GUI_Global.popup_mask = wx.Dialog(parent, wx.NewId(), 'MASK TOOL', size=(1230,800), style = wx.RESIZE_BORDER | wx.CLOSE_BOX)
    GUI_Global.popup_mask.Bind(wx.EVT_CLOSE, OnButton_mask_done)

    GUI_Global.maskpanel = wx.Panel(GUI_Global.popup_mask, size=(950,500), pos=(30,100))
    GUI_Global.maskpanel.SetBackgroundColour("#f8f8f8")
    GUI_Global.maskpanel.edit_mode = 0 # fiting = 1, masking = 2
      
    # Create the mpl Figure and FigCanvas objects. 
    GUI_Global.maskpanel.dpi = 100
    GUI_Global.maskpanel.fig = Figure((11.0, 4.5), dpi=GUI_Global.maskpanel.dpi)
    GUI_Global.maskpanel.fig.subplots_adjust(left=0.06)
    GUI_Global.maskpanel.fig.subplots_adjust(right=0.97)
    GUI_Global.maskpanel.fig.subplots_adjust(top=0.97)
    GUI_Global.maskpanel.canvas = FigCanvas(GUI_Global.maskpanel, -1, GUI_Global.maskpanel.fig)
    GUI_Global.maskpanel.canvas.mpl_connect('axes_enter_event', enter_graph)
    GUI_Global.maskpanel.canvas.mpl_connect('axes_leave_event', exit_graph)
    GUI_Global.maskpanel.axes = GUI_Global.maskpanel.fig.add_subplot(111)
     
    # Create the navigation toolbar, tied to the canvas
    GUI_Global.maskpanel.toolbar = NavigationToolbar(GUI_Global.maskpanel.canvas)

    GUI_Global.maskpanel.position_text = wx.StaticText(GUI_Global.maskpanel,-1, "", (40,110),(300,-1), wx.ALIGN_LEFT)

    # Checkboxes below the graph for selection of graph options   
    GUI_Global.maskpanel.cb_mask   = wx.CheckBox(GUI_Global.maskpanel, -1, "Show MASKS",                             style=wx.ALIGN_RIGHT)
    GUI_Global.maskpanel.cb_edit_m = wx.CheckBox(GUI_Global.maskpanel, -1, "Enter Mask (Exclusion) EDIT MODE",       style=wx.ALIGN_RIGHT)
    GUI_Global.maskpanel.cb_region = wx.CheckBox(GUI_Global.maskpanel, -1, "Show FIT REGIONS",                       style=wx.ALIGN_RIGHT)
    GUI_Global.maskpanel.cb_edit_f = wx.CheckBox(GUI_Global.maskpanel, -1, "Enter Fit (Inclusion) Region EDIT MODE", style=wx.ALIGN_RIGHT)

    # Checkboxes default values
    GUI_Global.maskpanel.cb_mask.SetValue(True)
    GUI_Global.maskpanel.cb_region.SetValue(True)

    # DONE BUTTON 
    button_mask_done = wx.Button(GUI_Global.maskpanel, id=-1, label="DONE", size=(150,30))
    button_mask_done.SetToolTip(wx.ToolTip("click to if you are ready with this step"))

    GUI_Global.maskpanel.Bind(wx.EVT_CHECKBOX, on_cbm, GUI_Global.maskpanel.cb_mask)
    GUI_Global.maskpanel.Bind(wx.EVT_CHECKBOX, on_edit_mode, GUI_Global.maskpanel.cb_edit_m)
    GUI_Global.maskpanel.Bind(wx.EVT_CHECKBOX, on_cbm, GUI_Global.maskpanel.cb_region)
    GUI_Global.maskpanel.Bind(wx.EVT_CHECKBOX, on_edit_mode, GUI_Global.maskpanel.cb_edit_f)
    GUI_Global.maskpanel.Bind(wx.EVT_BUTTON,   OnButton_mask_done, button_mask_done ) 


    ### Layout with box sizers ###
    GUI_Global.maskpanel.vbox = wx.BoxSizer(wx.VERTICAL)
    GUI_Global.maskpanel.hbox = wx.BoxSizer(wx.HORIZONTAL)
    fl = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
    
    
    # Add figure and position text to the frame
    GUI_Global.maskpanel.vbox.Add(GUI_Global.maskpanel.canvas,  1, wx.LEFT | wx.TOP | wx.GROW)
    GUI_Global.maskpanel.vbox.Add(GUI_Global.maskpanel.toolbar, 0, wx.EXPAND)
    
    GUI_Global.maskpanel.vbox.Add(GUI_Global.maskpanel.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
    GUI_Global.maskpanel.hbox.AddSpacer(30)
    GUI_Global.maskpanel.hbox.Add(GUI_Global.maskpanel.cb_mask,   0, border=3, flag=fl)
    GUI_Global.maskpanel.hbox.Add(GUI_Global.maskpanel.cb_edit_m, 0, border=3, flag=fl)
    GUI_Global.maskpanel.hbox.AddSpacer(40)
    GUI_Global.maskpanel.hbox.Add(GUI_Global.maskpanel.cb_region, 0, border=3, flag=fl)
    GUI_Global.maskpanel.hbox.Add(GUI_Global.maskpanel.cb_edit_f, 0, border=3, flag=fl)
    GUI_Global.maskpanel.hbox.AddSpacer(50)
    GUI_Global.maskpanel.hbox.Add(button_mask_done, 0, border=3, flag=fl)
    GUI_Global.maskpanel.hbox.AddSpacer(30)
    GUI_Global.maskpanel.vbox.Add(GUI_Global.maskpanel.position_text, 0, border=3, flag=fl)
    GUI_Global.maskpanel.vbox.AddSpacer(50)

    
    
    ### Create and show popup ###    
    GUI_Global.maskpanel.SetSizer(GUI_Global.maskpanel.vbox)
    GUI_Global.maskpanel.first = True
    create_popup_edit_mask(GUI_Global.maskpanel)
    

def on_cbm(event):
    # an checkbox botton causes a redraw
    draw_figure(True)
    
def on_edit_mode(event):
    # only one of the two EDIT buttons can be active
    if ((GUI_Global.maskpanel.edit_mode == 1) and
        (GUI_Global.maskpanel.cb_edit_m.GetValue() == True)):
        GUI_Global.maskpanel.cb_edit_f.SetValue(False)
        
    if ((GUI_Global.maskpanel.edit_mode == 2) and
        (GUI_Global.maskpanel.cb_edit_f.GetValue() == True)):
        GUI_Global.maskpanel.cb_edit_m.SetValue(False)

    GUI_Global.maskpanel.edit_mode = 0
    if (GUI_Global.maskpanel.cb_edit_f.GetValue() == True) :
        GUI_Global.maskpanel.edit_mode = 1      
        GUI_Global.maskpanel.cb_region.SetValue(True)

    if (GUI_Global.maskpanel.cb_edit_m.GetValue() == True) :
        GUI_Global.maskpanel.edit_mode = 2      
        GUI_Global.maskpanel.cb_mask.SetValue(True)
    
    draw_figure(True)
    

def draw_figure(reserve):
    # refreshing the window - boolean parameter controlls whether zoom level is
    # kept or not
    cur_xl = GUI_Global.maskpanel.axes.get_xlim()
    cur_yl = GUI_Global.maskpanel.axes.get_ylim()
   
    GUI_Global.x = range(len(GUI_Global.Y))
    # clear the axes and redraw the plot anew
    #
    GUI_Global.maskpanel.axes.clear()        
    GUI_Global.maskpanel.axes.set_xlabel('Wavelength (vac.) [microns]')
    GUI_Global.maskpanel.axes.set_ylabel('Flux')
    GUI_Global.maskpanel.axes.plot(GUI_Global.X, GUI_Global.Y, linewidth=0.5, 
        color=(0, 0, 0),)[0]        
    # change wish - adding if available FIT to mask display
    if (GUI_Global.panel.cb_model.IsChecked() and 
        GUI_Global.fit_available == True):
        GUI_Global.maskpanel.axes.plot(GUI_Global.XF, GUI_Global.YF,  
        linewidth=0.5, color=(1, 0, 0),)[0]        
    if (GUI_Global.panel.cb_diff.IsChecked() and 
        GUI_Global.fit_available == True):
        GUI_Global.maskpanel.axes.plot(GUI_Global.XD, GUI_Global.YD,  
        linewidth=0.5, color=(0, 0, 1),)[0]        
    
    t = GUI_Global.maskpanel.axes.get_ylim()
    GUI_Global.mask_min  = float(t[0])
    GUI_Global.mask_max  = float(t[1])-float(t[0])

    if GUI_Global.maskpanel.cb_mask.IsChecked():
        for m in GUI_Global.mask:
            # check if it is an EXCLUSION mask at all
            if (m.include == False) and (m.active  == True):
                if m.validxy == True:
                    # pixel type masks
                    m.startwl = GUI_Global.X[m.startx]
                    # for display purpose at least one pixel
                    if m.startx == m.endx :
                        if m.endx < len(GUI_Global.X) - 1:
                            m.endwl = GUI_Global.X[m.endx+1]
                        else:
                            m.endwl   = GUI_Global.X[m.endx]
                            m.startwl = GUI_Global.X[m.startx - 1]
                    else:
                        m.endwl = GUI_Global.X[m.endx]
                    color = PIX_EXCL_COLOR
                else :
                    color = WAV_EXCL_COLOR
                del m.box
                m.box = matplotlib.patches.Rectangle(
                    (m.startwl, GUI_Global.mask_min),
                    m.endwl - m.startwl,
                    GUI_Global.mask_max, color=color, alpha=0.3)
                GUI_Global.maskpanel.axes.add_patch(m.box)
          

    if GUI_Global.maskpanel.cb_region.IsChecked():
        for m in GUI_Global.mask:
            # check if it is an INCLUSION mask at all
            if (m.include == True) and (m.active  == True):
                # coordinates in "world coordinates" of plot
                # (x,y) of lower left corner and the delta_x and delta_y 
                # alpha gives transparency
                # INCLUSION masks exist only in wavelength regime
                if m.validxy == False:
                    del m.box
                    m.box  = matplotlib.patches.Rectangle(
                        (m.startwl, GUI_Global.mask_min),
                        m.endwl - m.startwl,
                        GUI_Global.mask_max, color=WAV_INCL_COLOR, alpha=0.3)
                    GUI_Global.maskpanel.axes.add_patch(m.box)

    a = GUI_Global.maskpanel.axes.get_ylim()
    if a[0] > 0:
        GUI_Global.maskpanel.axes.set_ylim((0.0, float(a[1])))

    m = GUI_Data.min_mean_median_max_of_data(GUI_Global.Y)
    if ((m[1] / m[2]) > 10) or ((m[1] / m[2]) < -3):
        GUI_Global.maskpanel.axes.set_ylim((-m[2], 5 * m[2]))
    #print "mask plot range calculated"
    #print GUI_Global.maskpanel.axes.get_ylim()
    GUI_Global.maskpanel.axes.grid(GUI_Global.panel.cb_grid.IsChecked())
    # restore size
    if ((GUI_Global.maskpanel.first == False) and 
        (reserve == True)) :
        GUI_Global.maskpanel.axes.set_xlim(cur_xl)
        GUI_Global.maskpanel.axes.set_ylim(cur_yl)
    else:
        # request by P.Ballester - if not in restore mode use same zoom as
        # that of main window (JIRA PIPE-4631)
        GUI_Global.maskpanel.axes.set_xlim(GUI_Global.panel.axes.get_xlim())
        GUI_Global.maskpanel.axes.set_ylim(GUI_Global.panel.axes.get_ylim())       
        GUI_Global.maskpanel.first = False

    #print "mask plot range at exit"
    #print GUI_Global.maskpanel.axes.get_ylim()

    # refresh all changes
    GUI_Global.maskpanel.canvas.draw()
    
#===============================================================================
#
#  routines to show the popup
#
def show_popup_mask(parent):
    if GUI_Global.data_available == False:
        GUI_Popup.no_data_message()
        return
    if (hasattr(GUI_Global,'mask') != True):
        GUI_Global.mask = []
    # at a new show mask case do not reserve the old zoom
    draw_figure(False)
    GUI_Global.maskpanel.cb_edit_f.SetValue(False)
    GUI_Global.maskpanel.cb_edit_m.SetValue(False)
    GUI_Global.maskpanel.edit_mode = 0
    GUI_Global.popup_mask.Show(True)
#===============================================================================
#
#  the EDIT sub-popup
#

def OnButton_mask_edit_done(e):
#    if GUI_Popup.check_mask_names() == False:
#        return
##    if GUI_Parameters.generate_and_write_mask_files() == False: 
##        GUI_Parameters.my_message("INFO: mask filenames not ready yet!\n\n"\
##        "goto \n   Settings->Working Directories and Files\nand add filenames.") 
##        return
    GUI_Global.popup_mask_edit.Show(False)
    GUI_Global.popup_mask_edit.Showed = False

def OnButton_mask_edit_delete(e):
    del GUI_Global.mask[GUI_Global.popup_mask_edit.current]
    GUI_Global.mask_changes = True
    GUI_Global.param_changes = True
    GUI_Global.popup_mask_edit.current = -99
    draw_figure(e)
    GUI_Global.popup_mask_edit.Show(False)
    GUI_Global.popup_mask_edit.Showed = False

def OnButton_mask_edit_save(e):
    #print "SAVE button @ %i" % (GUI_Global.popup_mask_edit.current)
    m = GUI_Global.mask[GUI_Global.popup_mask_edit.current]
    GUI_Global.mask_changes = True
    GUI_Global.param_changes = True

    if GUI_Global.popup_mask_edit.check_wavelength.GetValue():
        m.startwl = float(GUI_Global.popup_mask_edit.wls.GetValue())
        m.endwl   = float(GUI_Global.popup_mask_edit.wle.GetValue())
        m.validxy = False
        m.startx  = GUI_Data.pixel_of_data(m.startwl)
        m.endx    = GUI_Data.pixel_of_data(m.endwl)
    else:
        m.startx  = int(GUI_Global.popup_mask_edit.pxs.GetValue()) - 1
        m.endx    = int(GUI_Global.popup_mask_edit.pxe.GetValue()) - 1
        m.validxy = True
        m.startwl = GUI_Global.X[m.startx]
        m.endwl   = GUI_Global.X[m.endx]

    MergeMasksOfSameType()
    draw_figure(True)

def OnWlPixRadio(e):
    if GUI_Global.popup_mask_edit.check_pixel.GetValue():
        # toggle existing mask from wavelength to pixel
        if (GUI_Global.popup_mask_edit.current >= 0) :
            ii = GUI_Global.popup_mask_edit.current
            GUI_Global.mask[ii].validxy = True
        GUI_Global.mask[GUI_Global.popup_mask_edit.current].startx = (
            GUI_Data.pixel_of_data(
                GUI_Global.mask[GUI_Global.popup_mask_edit.current].startwl))
        GUI_Global.mask[GUI_Global.popup_mask_edit.current].endx = (
            GUI_Data.pixel_of_data(
                GUI_Global.mask[GUI_Global.popup_mask_edit.current].endwl))
    else:
        if (GUI_Global.popup_mask_edit.current >= 0) :
            ii = GUI_Global.popup_mask_edit.current
            GUI_Global.mask[ii].validxy = False
        GUI_Global.mask[GUI_Global.popup_mask_edit.current].startwl = \
            GUI_Global.X[GUI_Global.mask[\
            GUI_Global.popup_mask_edit.current].startx]
        GUI_Global.mask[GUI_Global.popup_mask_edit.current].endx = \
            GUI_Global.X[GUI_Global.mask\
            [GUI_Global.popup_mask_edit.current].endx]
    # change texts
    refresh_edit_text(-1)

def create_popup_edit_mask(parent):
    GUI_Global.popup_mask_edit = wx.Dialog(parent, wx.NewId(), 
        'MASK TOOL EDIT', size=(500, 320))

    GUI_Global.popup_mask_edit.mask_ed_nr = wx.StaticText(
        GUI_Global.popup_mask_edit, 
        -1, 'Mask Nr:  xx   Status:', pos=(20,10))

    GUI_Global.popup_mask_edit.mask_ed_fit = wx.StaticText(
        GUI_Global.popup_mask_edit, 
        -1, 'FIT REGION (inclusion)', pos=(220,10))
    GUI_Global.popup_mask_edit.mask_ed_fit.SetForegroundColour((0, 200,0)) 
    GUI_Global.popup_mask_edit.mask_ed_fit.SetBackgroundColour((0,0,255)) 

    GUI_Global.popup_mask_edit.mask_ed_mask = wx.StaticText(
        GUI_Global.popup_mask_edit, 
        -1, 'MASK (exclusion)', pos=(220,10))
    GUI_Global.popup_mask_edit.mask_ed_mask.SetForegroundColour((200,0 ,0)) 
    GUI_Global.popup_mask_edit.mask_ed_mask.SetBackgroundColour((0,0,255)) 

    GUI_Global.popup_mask_edit.mask_ed_wavelength = wx.StaticText(
        GUI_Global.popup_mask_edit, 
        -1, 'wavelength range', pos=(20,40))
    GUI_Global.popup_mask_edit.wls = wx.TextCtrl(GUI_Global.popup_mask_edit,  
         -1, "",size=(80, -1), style=wx.TE_RIGHT, pos=(150,35))
    GUI_Global.popup_mask_edit.wle = wx.TextCtrl(GUI_Global.popup_mask_edit, 
         -1, "", size=(80, -1), style=wx.TE_RIGHT, pos=(250,35))

    GUI_Global.popup_mask_edit.mask_ed_pixel = wx.StaticText(
        GUI_Global.popup_mask_edit, 
        -1, 'pixel range', pos=(20,40))
    GUI_Global.popup_mask_edit.pxs = wx.TextCtrl(GUI_Global.popup_mask_edit,  
         -1, "",size=(80, -1), style=wx.TE_RIGHT, pos=(150,35))
    GUI_Global.popup_mask_edit.pxe = wx.TextCtrl(GUI_Global.popup_mask_edit, 
         -1, "", size=(80, -1), style=wx.TE_RIGHT, pos=(250,35))
   
    GUI_Global.popup_mask_edit.check_wavelength = wx.RadioButton(
        GUI_Global.popup_mask_edit, 
        -1, "wavelength mode", pos=(70,113), style=wx.ALIGN_RIGHT|wx.RB_GROUP)
    GUI_Global.popup_mask_edit.check_wavelength.SetValue(True)
    GUI_Global.popup_mask_edit.Bind(wx.EVT_RADIOBUTTON, OnWlPixRadio,
        GUI_Global.popup_mask_edit.check_wavelength)

    GUI_Global.popup_mask_edit.check_pixel = wx.RadioButton(
        GUI_Global.popup_mask_edit, 
        -1, "chip pixel mode", pos=(270,113), style=wx.ALIGN_RIGHT)
    GUI_Global.popup_mask_edit.Bind(wx.EVT_RADIOBUTTON, OnWlPixRadio,
        GUI_Global.popup_mask_edit.check_pixel)
    
    button_mask_edit_done = wx.Button(GUI_Global.popup_mask_edit, id=-1, 
        label="EXIT / ABORT", pos=(100,200), size=(150,30))
    GUI_Global.popup_mask_edit.Bind(wx.EVT_BUTTON, OnButton_mask_edit_done, 
        button_mask_edit_done ) 

    button_mask_edit_delete = wx.Button(GUI_Global.popup_mask_edit, id=-1, 
        label="DELETE ENTRY", pos=(300,250), size=(150,30))
    GUI_Global.popup_mask_edit.Bind(wx.EVT_BUTTON, OnButton_mask_edit_delete, 
        button_mask_edit_delete ) 

    button_mask_edit_save = wx.Button(GUI_Global.popup_mask_edit, id=-1, 
        label="SAVE CHANGES", pos=(300,200), size=(150,30))
    GUI_Global.popup_mask_edit.Bind(wx.EVT_BUTTON, OnButton_mask_edit_save, 
        button_mask_edit_save ) 

    GUI_Global.popup_mask_edit.Showed  = False
    GUI_Global.popup_mask_edit.current = -99
    refresh_edit_text(-1)

def refresh_edit_text(mask_nr):
    if (mask_nr >= 0):
        GUI_Global.popup_mask_edit.current = mask_nr
    if GUI_Global.popup_mask_edit.current < 0:
        return

    m = GUI_Global.mask[GUI_Global.popup_mask_edit.current]

    # check which state this mask has
    if m.include == False :
        if m.validxy == True :
            # it is a pixel exclusion mask
            GUI_Global.popup_mask_edit.check_pixel.SetValue(True)
        else :
            # it is a wavelength exclusion mask
            GUI_Global.popup_mask_edit.check_wavelength.SetValue(True)
        GUI_Global.popup_mask_edit.check_pixel.Show(True)
    else :
        # it is a inclusion mask - always in wavelength space
        GUI_Global.popup_mask_edit.check_wavelength.SetValue(True)
        GUI_Global.popup_mask_edit.check_pixel.Show(False)
        
    if GUI_Global.popup_mask_edit.check_wavelength.GetValue():
        GUI_Global.popup_mask_edit.mask_ed_wavelength.Show(True)
        GUI_Global.popup_mask_edit.wls.Show(True)
        GUI_Global.popup_mask_edit.wle.Show(True)
        GUI_Global.popup_mask_edit.mask_ed_pixel.Show(False)
        GUI_Global.popup_mask_edit.pxs.Show(False)
        GUI_Global.popup_mask_edit.pxe.Show(False)
        GUI_Global.popup_mask_edit.wls.SetValue("%g" % (m.startwl))
        GUI_Global.popup_mask_edit.wle.SetValue("%g" % (m.endwl))
    else:
        GUI_Global.popup_mask_edit.mask_ed_wavelength.Show(False)
        GUI_Global.popup_mask_edit.wls.Show(False)
        GUI_Global.popup_mask_edit.wle.Show(False)
        GUI_Global.popup_mask_edit.mask_ed_pixel.Show(True)
        GUI_Global.popup_mask_edit.pxs.Show(True)
        GUI_Global.popup_mask_edit.pxe.Show(True)
        GUI_Global.popup_mask_edit.pxs.SetValue("%i" % (m.startx + 1))
        GUI_Global.popup_mask_edit.pxe.SetValue("%i" % (m.endx + 1))
    
    text = 'Mask Nr:  %3i   Status:' % (GUI_Global.mask.index(m))
    GUI_Global.popup_mask_edit.mask_ed_nr.SetLabel(text)
    if m.include == True:
        GUI_Global.popup_mask_edit.mask_ed_fit.Show(True)
        GUI_Global.popup_mask_edit.mask_ed_mask.Show(False)
    else:
        GUI_Global.popup_mask_edit.mask_ed_fit.Show(False)
        GUI_Global.popup_mask_edit.mask_ed_mask.Show(True)

