# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Feb 26 2014)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc

###########################################################################
## Class FitKernelDialog
###########################################################################

class FitKernelDialog ( wx.Dialog ):
        
        def __init__( self, parent ):
                wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"FIT KERNEL DEFINITION", pos = wx.DefaultPosition, size = wx.DefaultSize, style = wx.DEFAULT_DIALOG_STYLE )
                
                self.SetSizeHintsSz( wx.Size( 500,-1 ), wx.DefaultSize )
                
                dsizer = wx.BoxSizer( wx.VERTICAL )
                
                boxsizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.check_box = wx.CheckBox( self, wx.ID_ANY, u"Boxcar profile", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.check_box.SetMinSize( wx.Size( 150,-1 ) )
                
                boxsizer.Add( self.check_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.relres_box_t = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_RIGHT )
                self.relres_box_t.SetMinSize( wx.Size( 100,-1 ) )
                
                boxsizer.Add( self.relres_box_t, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.txt_box = wx.StaticText( self, wx.ID_ANY, u"FWHM relative to slit width ", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_box.Wrap( -1 )
                boxsizer.Add( self.txt_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                
                dsizer.Add( boxsizer, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )
                
                voigtsizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.check_kernmode = wx.CheckBox( self, wx.ID_ANY, u"Voigt profile approximation", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_RIGHT )
                self.check_kernmode.SetMinSize( wx.Size( 100,-1 ) )
                
                voigtsizer.Add( self.check_kernmode, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5 )
                
                
                dsizer.Add( voigtsizer, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )
                
                gaussizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.check_gauss = wx.CheckBox( self, wx.ID_ANY, u"Gaussian profile", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.check_gauss.SetMinSize( wx.Size( 150,-1 ) )
                
                gaussizer.Add( self.check_gauss, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.res_gauss_t = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_RIGHT )
                self.res_gauss_t.SetMinSize( wx.Size( 100,-1 ) )
                
                gaussizer.Add( self.res_gauss_t, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.txt_gaus = wx.StaticText( self, wx.ID_ANY, u"Initial FWHM [pixel]", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_gaus.Wrap( -1 )
                gaussizer.Add( self.txt_gaus, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                
                dsizer.Add( gaussizer, 1, wx.EXPAND, 5 )
                
                lorentzsizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.check_lorentz = wx.CheckBox( self, wx.ID_ANY, u"Lorentzian profile", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.check_lorentz.SetMinSize( wx.Size( 150,-1 ) )
                
                lorentzsizer.Add( self.check_lorentz, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.res_lorentz_t = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_RIGHT )
                self.res_lorentz_t.SetMinSize( wx.Size( 100,-1 ) )
                
                lorentzsizer.Add( self.res_lorentz_t, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.txt_lorentz = wx.StaticText( self, wx.ID_ANY, u"Initial FWHM [pixel]", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_lorentz.Wrap( -1 )
                lorentzsizer.Add( self.txt_lorentz, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                
                dsizer.Add( lorentzsizer, 1, wx.EXPAND, 5 )
                
                fwhmsizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.txt_fwhm = wx.StaticText( self, wx.ID_ANY, u"Size of Gaussian/Lorentzian/Voigt kernel", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_fwhm.Wrap( -1 )
                fwhmsizer.Add( self.txt_fwhm, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.kernfac_t = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_RIGHT )
                fwhmsizer.Add( self.kernfac_t, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.txt_fwhm = wx.StaticText( self, wx.ID_ANY, u"[FWHM]", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_fwhm.Wrap( -1 )
                fwhmsizer.Add( self.txt_fwhm, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                
                dsizer.Add( fwhmsizer, 1, wx.EXPAND, 5 )
                
                varsizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.check_varkern = wx.CheckBox( self, wx.ID_ANY, u"Kernel variable with wavelength", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_RIGHT )
                varsizer.Add( self.check_varkern, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                
                dsizer.Add( varsizer, 1, wx.EXPAND, 5 )
                
                filesizer = wx.BoxSizer( wx.HORIZONTAL )
                
                self.txt_file = wx.StaticText( self, wx.ID_ANY, u"User defined kernel file", wx.DefaultPosition, wx.DefaultSize, 0 )
                self.txt_file.Wrap( -1 )
                filesizer.Add( self.txt_file, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )
                
                self.kernel_file_t = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
                filesizer.Add( self.kernel_file_t, 1, wx.ALL, 5 )
                
                
                dsizer.Add( filesizer, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )
                
                btnsizer = wx.StdDialogButtonSizer()
                self.btnsizerOK = wx.Button( self, wx.ID_OK )
                btnsizer.AddButton( self.btnsizerOK )
                self.btnsizerCancel = wx.Button( self, wx.ID_CANCEL )
                btnsizer.AddButton( self.btnsizerCancel )
                btnsizer.Realize();
                
                dsizer.Add( btnsizer, 1, wx.EXPAND, 5 )
                
                
                self.SetSizer( dsizer )
                self.Layout()
                dsizer.Fit( self )
                
                self.Centre( wx.BOTH )
        
        def __del__( self ):
                pass
        


