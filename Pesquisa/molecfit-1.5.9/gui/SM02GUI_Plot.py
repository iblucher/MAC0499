#!/usr/bin/python -B
#
# module for the main window plot
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
import SM02GUI_Help        as GUI_Help
import SM02GUI_Data        as GUI_Data
import SM02GUI_Help_Dialog_ToolTip_Texts as TextPool

import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar


def create_main_panel(parent):
    
    # create the matplotlib plotting area and the buttons around it
    GUI_Global.panel = wx.Panel(parent, size=(1230,520), pos=(30,100))
    GUI_Global.panel.SetBackgroundColour("#f8f8f8")
      
    # Create the mpl Figure and FigCanvas objects. 
    GUI_Global.panel.dpi = 100
    GUI_Global.panel.fig = Figure((12.0, 4.5), dpi=GUI_Global.panel.dpi) #, subplotpars=[0.06, 0.10, 0.97, 0.97, 0.20, 0.20])
    GUI_Global.panel.fig.subplots_adjust(left=0.06)
    GUI_Global.panel.fig.subplots_adjust(right=0.97)
    GUI_Global.panel.fig.subplots_adjust(top=0.97)
    GUI_Global.panel.fig.subplots_adjust(bottom=0.15)
    GUI_Global.panel.canvas = FigCanvas(GUI_Global.panel, -1, GUI_Global.panel.fig)
    GUI_Global.panel.axes = GUI_Global.panel.fig.add_subplot(111)

    # Create the navigation toolbar, tied to the canvas
    GUI_Global.panel.toolbar = NavigationToolbar(GUI_Global.panel.canvas)

    GUI_Global.panel.position_text = wx.StaticText(GUI_Global.panel,-1, " \n ", (40,110),(-1,-1), wx.ALIGN_LEFT)
             
    # checkboxes below the graph for selection of graph options   
    GUI_Global.panel.cb_grid   = wx.CheckBox(GUI_Global.panel, -1, "Show: Grid",        style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_legend = wx.CheckBox(GUI_Global.panel, -1, "Legend",            style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_raw    = wx.CheckBox(GUI_Global.panel, -1, "DATA",              style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_model  = wx.CheckBox(GUI_Global.panel, -1, "MODEL",             style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_diff   = wx.CheckBox(GUI_Global.panel, -1, "DIFFERENCE",        style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_mask   = wx.CheckBox(GUI_Global.panel, -1, "EXCLUSION REGIONS", style=wx.ALIGN_RIGHT)
    GUI_Global.panel.cb_region = wx.CheckBox(GUI_Global.panel, -1, "FIT REGIONS",       style=wx.ALIGN_RIGHT)
    GUI_Global.panel.bp_region = wx.CheckBox(GUI_Global.panel, -1, "BAD PIXELS",        style=wx.ALIGN_RIGHT)

    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_grid)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_legend)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_raw)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_model)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_diff)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_mask)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.cb_region)
    parent.Bind(wx.EVT_CHECKBOX, on_cb, GUI_Global.panel.bp_region)

    # Initial default values
    GUI_Global.panel.cb_grid.SetValue(True)
    GUI_Global.panel.cb_raw.SetValue(True)
    GUI_Global.panel.cb_legend.SetValue(True)
    GUI_Global.panel.cb_mask.SetValue(True)
    GUI_Global.panel.cb_region.SetValue(True)
    GUI_Global.panel.bp_region.SetValue(True)
    
       
    GUI_Global.panel.vbox = wx.BoxSizer(wx.VERTICAL)
    GUI_Global.panel.hbox = wx.BoxSizer(wx.HORIZONTAL)
    fl = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
    
    GUI_Global.panel.vbox.Add(GUI_Global.panel.canvas,  1, wx.LEFT | wx.TOP | wx.GROW)
    GUI_Global.panel.vbox.Add(GUI_Global.panel.toolbar, 0, wx.EXPAND)
    
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_grid,   0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_legend, 0, border=3, flag=fl)
    GUI_Global.panel.hbox.AddSpacer(30)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_raw,    0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_model,  0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_diff,   0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_mask,   0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.cb_region, 0, border=3, flag=fl)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.bp_region, 0, border=3, flag=fl)
    GUI_Global.panel.hbox.AddSpacer(30)
    GUI_Global.panel.hbox.Add(GUI_Global.panel.position_text, 0, border=3, flag=fl)
        
    GUI_Global.panel.vbox.Add(GUI_Global.panel.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
    GUI_Global.panel.SetSizer(GUI_Global.panel.vbox)
    GUI_Global.panel.Fit() # required in macos with wxpython3.0
    GUI_Global.panel.first = True
    
    
    create_popup_apply_plot()
        
def draw_figure(reserve):
    #   draw or refresh the plot - the boolean variable controls reservation /
    #   ot non reservation of zoom state
    cur_xl = GUI_Global.panel.axes.get_xlim()
    cur_yl = GUI_Global.panel.axes.get_ylim()
    
    GUI_Global.x = range(len(GUI_Global.Y))

    # clear the axes and redraw the plot 
    #
    GUI_Global.panel.axes.clear()        
    GUI_Global.panel.axes.set_xlabel('Wavelength (vac.) [microns]')
    GUI_Global.panel.axes.set_ylabel('Flux')
    GUI_Global.panel.axes.grid(GUI_Global.panel.cb_grid.IsChecked())
    if GUI_Global.panel.cb_raw.IsChecked():
        if (GUI_Global.panel.first == True) :
            GUI_Global.panel.axes.plot(GUI_Global.X, GUI_Global.Y,linewidth=0.5, 
                    color=(1.0, 1.0, 1.0), label='Raw')
        else:
            GUI_Global.panel.axes.plot(GUI_Global.X, GUI_Global.Y,linewidth=0.5, 
                    color=(0, 0, 0), label='Raw')
        
    if (GUI_Global.panel.cb_model.IsChecked() and 
        GUI_Global.fit_available == True):
        GUI_Global.panel.axes.plot(GUI_Global.XF, GUI_Global.YF, linewidth=0.5, 
        color=(1, 0, 0), label='Model')
    if (GUI_Global.panel.cb_diff.IsChecked() and 
        GUI_Global.fit_available == True):
        GUI_Global.panel.axes.plot(GUI_Global.XD, GUI_Global.YD, linewidth=0.5, 
        color=(0, 0, 1), label='Raw - Model')

    if (GUI_Global.panel.cb_diff.IsChecked() and 
        GUI_Global.corrected_available == True):
        GUI_Global.panel.axes.plot(GUI_Global.XC, GUI_Global.YC, linewidth=0.5, 
        color=(1, 0, 1), label='Corrected')

    inc_labeled = False
    wex_labeled = False
    pex_labeled = False
    if GUI_Global.panel.cb_mask.IsChecked():
        if hasattr(GUI_Global,"mask") == False:
            GUI_Global.mask = []
        for m in GUI_Global.mask:
            # check if it is an EXCLUSION mask at all
            if ((m.include == False) and (m.active  == True)):
                # coordinates in "world coordinates" of plot
                # (x,y) of lower left corner and the delta_x and delta_y 
                # alpha gives transparency
                # INCLUSION masks exist in wavelength regime AND in pixel regime
                if m.validxy == False:
                    del m.box
                    m.box = matplotlib.patches.Rectangle(
                        (m.startwl, GUI_Global.mask_min),
                        m.endwl - m.startwl,
                        GUI_Global.mask_max, color='#ff0000', alpha=0.3)
                    if not wex_labeled:
                        m.box.set_label('Wavelength exclusion')
                        wex_labeled = True
                    GUI_Global.panel.axes.add_patch(m.box)
                else:
                    del m.box
                    m.startx = max(m.startx, 0)
                    m.endx = min(m.endx, len(GUI_Global.X) - 1)
                    m.startwl = GUI_Global.X[m.startx]
                    m.endwl   = GUI_Global.X[m.endx]
                    m.box = matplotlib.patches.Rectangle(
                        (m.startwl,GUI_Global.mask_min),
                        m.endwl - m.startwl,
                        GUI_Global.mask_max, color='#fe9d00', alpha=0.3)
                    if not pex_labeled:
                        m.box.set_label('Pixel exclusion')
                        pex_labeled = True
                    GUI_Global.panel.axes.add_patch(m.box)
                              
    if GUI_Global.panel.cb_region.IsChecked():
        for m in GUI_Global.mask:
            # check if it is an INCLUSION mask at all
            if (m.include == True) and (m.active  == True):
                # coordinates in "world coordinates" of plot
                # (x,y) of lower left corner and the delta_x and delta_y 
                # alpha gives transparency
                # INCLUSION masks exist only in wavelength regime
                if m.validxy == False:
                    del m.box
                    m.box = matplotlib.patches.Rectangle(
                        (m.startwl, GUI_Global.mask_min),
                        m.endwl - m.startwl,
                        GUI_Global.mask_max, color='#00FF00', alpha=0.3)
                    if not inc_labeled:
                        m.box.set_label('Fit region')
                        inc_labeled = True
                    GUI_Global.panel.axes.add_patch(m.box)

    # Bad pixel regions #
    if GUI_Global.panel.bp_region.IsChecked() and \
       hasattr(GUI_Global, 'badpixregions'):
        labeled_once = False
        for badpixregion in GUI_Global.badpixregions:
            rect = matplotlib.patches.Rectangle(
                (badpixregion['start'], GUI_Global.mask_min),
                badpixregion['end'] - badpixregion['start'],
                GUI_Global.mask_max, color='#000000', alpha=0.8)
            if not labeled_once:
                labeled_once = not labeled_once
                rect.set_label('Bad pixels')
            GUI_Global.panel.axes.add_patch(rect)
    #####################

    draw_legend(ax=GUI_Global.panel.axes)
    # connect button press to subroutine onclick for measurements
    # with a self.canvas.mpl_disconnect(cid) you can disactivate at any
    # time the measurement stream
    GUI_Global.panel.canvas.mpl_connect('axes_enter_event'   , enter_graph)
    GUI_Global.panel.canvas.mpl_connect('axes_leave_event'   , exit_graph)
    # calculate meaningfull Y cuts
    m = GUI_Data.min_mean_median_max_of_data(GUI_Global.Y)
    #rescaler disabled at the moment
    #print m
    #if m[0] < 0:
    #    GUI_Global.mask_min = int(m[0] - 2 * (m[3] - m[0]))
    #else:
    #    GUI_Global.mask_min = int(- 2 * (m[3] - m[0]))
    #GUI_Global.mask_max = int(m[3] + 5 * (m[3] - m[0]))
    #print "mask Y size %i %i" % (GUI_Global.mask_min, GUI_Global.mask_max)

    a = GUI_Global.panel.axes.get_ylim()
    if a[0] > 0:
        a = tuple([0.0,float(a[1])])
        GUI_Global.panel.axes.set_ylim(a)

    if ((m[1] / m[2]) > 10) or ((m[1] / m[2]) < -3):
        GUI_Global.panel.axes.set_ylim([-m[2],+5*m[2]])
    # restore original zoom state if Parameter is not False
    if ((GUI_Global.panel.first == False) and 
        (reserve != False)):
        GUI_Global.panel.axes.set_xlim(cur_xl)
        GUI_Global.panel.axes.set_ylim(cur_yl)
    else:
        GUI_Global.panel.first = False

    # refresh all changes
    if GUI_Global.data_available == True :
        GUI_Global.panel.canvas.draw()
    
def enter_graph(event):
    # when entering the graph the motion is monitored to display the data
    GUI_Global.panel.cid_move = GUI_Global.panel.canvas.mpl_connect(
        'motion_notify_event', move_graph)

def exit_graph(event):
    # on exit the data value refresh is turned off
    GUI_Global.panel.canvas.mpl_disconnect(GUI_Global.panel.cid_move)
    GUI_Global.panel.position_text.SetLabel("")    

def move_graph(event):
    # refresh the data display itself
    if event.xdata is not None and event.ydata is not None:
        if GUI_Data.pixel_of_data(event.xdata) is not None:
            text = (('wavelength: %.5g  y: %.5g\n' +
             'pixel:      %.4i      value: %.4g\n') % 
             (event.xdata,event.ydata,GUI_Data.pixel_of_data(event.xdata)+1,
              GUI_Global.Y[GUI_Data.pixel_of_data(event.xdata)]))
        else:
            text = 'wavelength: %g     y: %g' % (event.xdata,event.ydata)
        GUI_Global.panel.position_text.SetLabel(text)    

def on_cb(e):
    if GUI_Global.panel.cb_model == e.GetEventObject():
        if GUI_Global.fit_available == False:
            GUI_Help.refresh_error(TextPool.error_002)
            GUI_Global.panel.cb_model.SetValue(False)
    # an checkbox botton causes a redraw
    draw_figure(True)

#===============================================================================
# subplot for the APPLY FIT feature
def OnButton_apply_done(e):
    # DONE 
    GUI_Global.popup_aplot.Show(False)

def create_popup_apply_plot():
    GUI_Global.popup_aplot = wx.Dialog(GUI_Global.frame, wx.NewId(), 'APPLIED CORRECTION', 
        size=(1000, 500))

    GUI_Global.edit_started = False

    GUI_Global.applypanel = wx.Panel(GUI_Global.popup_aplot, size=(1000,500), \
        pos=(30,100))
    GUI_Global.applypanel.SetBackgroundColour("#f8f8f8")
      
    GUI_Global.applypanel.dpi = 100
    GUI_Global.applypanel.fig = Figure((10.0, 3.5), \
        dpi=GUI_Global.applypanel.dpi)
    GUI_Global.applypanel.fig.subplots_adjust(left=0.08)
    GUI_Global.applypanel.fig.subplots_adjust(right=0.97)
    GUI_Global.applypanel.fig.subplots_adjust(top=0.97)
    GUI_Global.applypanel.fig.subplots_adjust(bottom=0.15)
    GUI_Global.applypanel.canvas = FigCanvas(GUI_Global.applypanel, -1, 
        GUI_Global.applypanel.fig)
        
    GUI_Global.applypanel.axes = GUI_Global.applypanel.fig.add_subplot(111)
        
    GUI_Global.applypanel.canvas.mpl_connect('axes_enter_event'   , enter_graph)
    GUI_Global.applypanel.canvas.mpl_connect('axes_leave_event'   , exit_graph)
 
    
    # Create the navigation toolbar, tied to the canvas
    #
    GUI_Global.applypanel.toolbar = NavigationToolbar(GUI_Global.applypanel.canvas)
    #
    # Layout with box sizers
    #
        
    GUI_Global.applypanel.vbox = wx.BoxSizer(wx.VERTICAL)
    GUI_Global.applypanel.vbox.Add(GUI_Global.applypanel.canvas, 1, 
        wx.LEFT | wx.TOP )# | wx.GROW)
    GUI_Global.applypanel.vbox.Add(GUI_Global.applypanel.toolbar, 0, wx.EXPAND)
    #self.vbox.AddSpacer(10)
        
    GUI_Global.applypanel.hbox = wx.BoxSizer(wx.HORIZONTAL)
    fl = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
    GUI_Global.applypanel.vbox.Add(GUI_Global.applypanel.hbox, 0, 
        flag = wx.ALIGN_LEFT | wx.TOP)

    GUI_Global.applypanel.hbox.AddSpacer(20)
    GUI_Global.applypanel.position_text = wx.StaticText(
        GUI_Global.applypanel,-1, "", (40,110),(300,-1), wx.ALIGN_LEFT)
    GUI_Global.applypanel.hbox.Add(GUI_Global.applypanel.position_text, 0, 
        border=3, flag=fl)

     
    # DONE BUTTON 
    button_apply_done = wx.Button(GUI_Global.applypanel, id=-1, label="DONE", 
        pos=(1100,660), size=(150,30))
    GUI_Global.applypanel.Bind(wx.EVT_BUTTON, OnButton_apply_done, 
        button_apply_done ) 
    button_apply_done.SetToolTip(wx.ToolTip(
        "click to if you are ready with this step"))
 
    GUI_Global.applypanel.hbox.AddSpacer(70)
    GUI_Global.applypanel.hbox.Add(button_apply_done, 0, 
        border=3, flag=fl)
    GUI_Global.applypanel.SetSizer(GUI_Global.applypanel.vbox)
    GUI_Global.applypanel.first = True

def draw_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    if GUI_Global.panel.cb_legend.IsChecked() and handles:
        l = ax.legend(handles, labels)
        plt.setp(l.get_texts(), fontsize='small')
    
def draw_apply_figure():
    ##GUI_Global.x = range(len(GUI_Global.YC))

    # clear the axes and redraw the plot anew
    #
    try:
        wx.EndBusyCursor()
    except:
        pass
    GUI_Global.applypanel.axes.clear()        
    GUI_Global.applypanel.axes.set_xlabel('Wavelength (vac.) [microns]')
    GUI_Global.applypanel.axes.set_ylabel('Flux')
    
    GUI_Global.applypanel.axes.plot(GUI_Global.XC, GUI_Global.YC, linewidth=0.5, 
        color=(0.2, 0.5, 0.2), label='Corrected')
    GUI_Global.applypanel.axes.plot(GUI_Global.X, GUI_Global.Y, linewidth=0.5, 
        color=(0, 0, 0), label='Raw')
    if GUI_Global.fit_available == True:
        GUI_Global.applypanel.axes.plot(GUI_Global.XF, GUI_Global.YF, \
            linewidth=0.5, color=(1, 0, 0), label='Model')
    GUI_Global.applypanel.axes.grid(GUI_Global.panel.cb_grid.IsChecked())

    a = GUI_Global.applypanel.axes.get_ylim()
    if a[0] > 0:
        a = tuple([0.0,float(a[1])])
        GUI_Global.applypanel.axes.set_ylim(a)
    m = GUI_Data.min_mean_median_max_of_data(GUI_Global.YC)
    if ((m[1] / m[2]) > 10) or ((m[1] / m[2]) < -3):
        GUI_Global.applypanel.axes.set_ylim([-m[2],+5*m[2]])
    GUI_Global.applypanel.axes.relim()

    draw_legend(ax=GUI_Global.applypanel.axes)
             
    # refresh all changes
    GUI_Global.applypanel.canvas.draw()


