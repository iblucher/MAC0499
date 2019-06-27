#!/usr/bin/python -B
#
# Collection of texts for tooltips help and dialogs/popups of the GUI to make
# maintainance / modifications / updates easier
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

#===============================================================================
#
# for SM02GUI_Mask
#
# full MENU structure of main window/frame
menu_001 = 'Save (Configuration && Masks)'
tip_001  = 'Save the current set of parameters'
menu_save_error = 'No parameters available yet !\nSAVE ignored / canceled'
menu_001b= 'Save As (Configuration && Masks)'
tip_001b = 'Save the current set of parameters asking for optional new name'
menu_002 = 'Import (Configuration File)'
tip_002  = 'Read a set of parameters from a molecfit configuration file'
menu_02c = 'Set Default Parameters'
tip_02c  =('Set a set of default parameters - for some parameters there \n'+
           'are no default values')
menu_02b = 'Import Name for RAW Spectrum (to this parameter set)'
tip_02b  = 'Import a Spectrum for Fitting'
menu_003 = 'Export Fit Spectrum'
tip_003  = 'Export the fit spectrum'
menu_004 = 'Exit'
tip_004  = 'Terminate the program'

menu_005 = 'Mask Tool'
tip_005  = 'Start a masking tool'

menu_006 = 'Molecule Abundances'
tip_006  = 'Set the abundances of the molecules'
menu_007 = 'Wavelength Grid'
tip_007  = 'Set && define the fitting of the wavelength correction'
menu_008 = 'Resolution / Fitting Kernel'
tip_008  = 'Set parameters for the fitting kernel'
menu_009 = 'Continuum of the Spectrum'
tip_009  = 'Set && define the fitting of the continuum'
menu_010 = 'Telescope Background'
tip_010  = 'Define the parameters for the telescope background'
menu_011 = 'Working Directory and Files'
tip_011  = 'Define the primary files and directories'
menu_012 = 'Other Settings / Preferences'
tip_012  = 'Define all the other preferences'

menu_013 = 'Execute / Run Fit'
tip_013  = 'Run the final fit'
menu_013a= 'Show CONSOLE'
tip_013a = 'Show a text console for underlying child coldes like molecfit'
menu_013b= 'Hide CONSOLE'
tip_013b = 'Hide the text console'

menu_014 = 'Quick / Short Help'
tip_014  = 'Short Instructions'
menu_015 = 'Tutorial and GUI Help'
tip_015  = 'Extended User Help for the GUI'
menu_016 = 'MOLECFIT Manual'
tip_016  = 'Pop up the full manual of the underlying MOLECFIT code'
menu_017 = 'About'
tip_017  = 'Information about the GUI'

dialog_000 = ('molecfit executable not found.\n'+
            'Either specify it in the PATH variable before start or set it\n'+
            'via the menu of the GUI before executing fit.')

dialog_001 = ('This GUI is the frontend for the SM-02 Molecfit program.\n\n' +
         'M. Barden, W. Kausch, S. Kimeswenger (head),\n' +
         'A.M. Jones, S. Noll, C. Szyszka\n' +
         'Austrian ESO In-Kind Team Innsbruck\n\n' +
         'Institute for Astro and Particle Physics\n' +
         'University of Innsbruck, Austria\n\n' +
         '(c) ESO 2013')

dialog_002 = ('short WORKFLOW INSTRUCTIONS: \nthe fitting procedure '+
    'contains seven steps\n\n'+
    '   1) Masking wavelength regions / borders for the exclusion \n'+
    '        of the fit by using the masking tool\n'+
    '   2) Modifying list of molecules to be fitted by adding or\n'+
    '        deleting, setting start value for the column density\n'+
    '   3) Setting method / parameters for the wavelength grid fit\n'+
    '   4) Setting profile(s) for the fit of the spectral resolution\n'+
    '   5) Setting method for the continuum fit, either polynomial \n'+
    '        or spline/linear fit\n'+
    '   6) Setting value for the telescope emission (grey body \n'+
    '        emissivity [0,1]; for the emission spectra only)\n'+
    '   7) Execute the fit\n'+
    '')

dialog_003 = 'Really close the application ?\n please confirm !\n\n'

dialog_004 = 'Open / Import a configuration file'

dialog_005 = 'Export configuration to file'

dialog_006 = 'Warning: nothing was selected\n          - command ignored'

dialog_007 = 'Export fitted spectrum as 1D FITS file'

dialog_009 = 'Import a spectrum to fit molecular data'
#===============================================================================
# for SM02GUI_Plot
error_002  = 'No FIT data available up to now'
#===============================================================================
# for SM02GUI_Popups
# in popup_molec
tip_018  = ('Activate these for adding the molecule in the ' +
            'calculation without automatic variation of abundance')
tip_019  = ('Activate these for full handling of the molecule ' +
            'in a Chi2 fit')
tip_020  = ('abundances are given relative to the standard '+
            'atmosphere. \nE.g. CO2 was 6% lower a decade ago '+
            'and has also seasonal variations.')

tip_021  = ('Additional trace gases can be activated for calculation '+
            'without the possibility of a Chi2 fit.')
tip_022 = '''\
PWV value in mm for the input water vapour profile.
The merged profile will be scaled to this value if it is > 0.'''

error_001  = ('NO active mask or fit region was found!\n' +
              'Default fitting whole spectrum is used.')

help_001   =('Molecules are defined here for the radiative transfer model.\n'+
             'Their abundance are relative to the standard atmosphere model.\n'+
             'For details and references of this model refer to the MOLECFIT'+
             ' manual.\n'+
             'E.g. \n'
             'H2O and water vapour varies a lot due to weather conditions.\n'+
             'CO2 was about 6% lower a decade ago (global warming &&\n'+
             'greenhouse effects) and has seasonal variations.\n\n'
             'If the `calculation` is activated the molecule is added\n'
             'according to the given abundance.\n\n'+
             'If additionally the `fit` is actived, the code varies the\n'+
             'abundance in a Chi2 fit. The FIT REGIONS defined by the mask\n'+
             'tool are used for that. If no regions are defined, the whole\n'+
             'spectrum is taken into account.\n'+
             'We suggest not to fit all possible molecules at once. The\n'+
             'degree of freedom and a fit over the whole spectrum (= all \n'+
             'features) might lead to ambivalent and inaccurate results\n\n'+
             'H2O additionally allows to set a PWV (column) value in [mm].\n'+
             'This setting, if activated with the checkbox, overrides other\n'+
             'settings for H2O relative values. Still fitting can be activated\n'+
             'or not individually. This feature allows to insert values\n'+
             'directly given by other facilities at the observatory (e.g.\n' +
             'by a radiometer)' )

help_002   =('Molecules are defined here for the radiative transfer model.\n'+
             'Their abundance are relative to the standard atmosphere model.\n'+
             'For details and references of this model refer to the MOLECFIT'+
              ' manual.\n'+
              'If the `calculation` is activated the molecule is added\n'
              'according to the given abundance.')

help_003   =('Wavelength correction can be obtained / activated using\n'+
             'polynomials. The programme automatically calculates the '+
             'function.\n'
             'The user can only provide a guess for the 0th order parameter\n'+
             '(= offset = parallel shift). Higher order parameters are\n'+
             'too critical. Thus, to avoid unstable solutions, this is not\n'+
             'mutable manually. This method only finds / corrects for minor\n'+
             'inaccuracy in the wavelength calibration to fit lines and\n'+
             'profiles properly. It does not `replace` an appropriate\n'+
             'wavelength correction by the pipeline beforehand')

dialog_008  = ('ERROR: no data available - define input data first!')

dialog_009  = ('ERROR: no configuration / parameters available\n'+
               ' - load a molecfit parameter file or load a default set!')

dialog_010  = ('WARNING: \n'+
               '  Default (working) directory structure is that in the \n'+
               '  home of molecfit (see molecfit manual for details).\n'+
               '  To change those go to\n         \'Menu->Settings->Working '
               '  Directory and Files\'\n\n'+
               '  Due to the large variety of instruments and pipelines, and\n'+
               '  due to the fact that some information is not available\n'+
               '  as header keywords, in some pipelines a few parameters do\n'+
               '  not have default values.\n'+
               '  You always have to set them if you do not IMPORT them \n'+
               '  via the \'Menu->File->Import Configuration File\'.\n'+
               '  You find them in:\n'+
               '         \'Menu->Settings->Other Parameters\'')

dialog_011  = ('INFO: \n'+
               '  The FIT RESULT file exists AND is younger than the\n'+
               '  parameter file. Should fit be FORCED?\n'+
               '  Otherwhise the existing one is just displayed.')

dialog_012  = ('INFO: \n'+
               '  The APPLY RESULT file exists AND is younger than parameter\n'+
               '  file and fit file. Should APPLY CORRECTION be FORCED?\n'+
               '  Otherwhise the existing one is just displayed.')

dialog_013  = ('QUESTION: \n'+
               '  Old masks and fit ranges found.\n'+
               '  Should I keep the fit ranges and the wavelength masks ?\n'+
               '  They will be merged with those defined in the setup file.\n'+
               '  (info: the pixel masks come with the new input - the old\n'+
               '         ones are always flushed)')


