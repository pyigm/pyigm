"""
#;+ 
#; NAME:
#; igmguesses
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for LineIDs and Initial guesses in IGM spectra with QT
#;      Likely only useful for lowz-IGM
#;   14-Aug-2015 by JXP
#;-
#;- NT: New version using linetools' AbsComponent
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import warnings
import copy
import pdb

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy.units import Quantity
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.analysis import voigt as lav
from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu
from linetools.guis import utils as ltgu
from linetools.guis import line_widgets as ltgl
from linetools.guis import simple_widgets as ltgsm
from linetools.guis import spec_widgets as ltspw
from linetools import utils as ltu

# Global variables; defined as globals mainly to increase speed and convenience
c_kms = const.c.to('km/s').value
COLOR_MODEL = '#999966'
COLOR_EXTMODEL =  '#D4AC0D'
COLOR_RELIABLE = 'g'
COLOR_PROBABLE = 'b'
COLOR_UNCERTAIN = 'r'
COLORS = ['#0066FF','#339933','#CC3300','#660066','#FF9900','#B20047']
zero_coord = SkyCoord(ra=0.*u.deg, dec=0.*u.deg)  # Coords

# GUI for fitting LLS in a spectrum
class IGMGuessesGui(QtGui.QMainWindow):
    """ GUI to identify absorption features and provide reasonable
        first guesses of (z, logN, b) for subsequent Voigt profile
        fitting.

        v0.5
        30-Jul-2015 by JXP
    """
    def __init__(self, ispec, parent=None, previous_file=None, 
        srch_id=True, outfil=None, fwhm=None,
        plot_residuals=True,n_max_tuple=None, min_strength=None,
                 min_ew=None, vlim_disp=None, external_model=None):
        QtGui.QMainWindow.__init__(self, parent)
        """
        ispec : str
            Name of the spectrum file to load
        previous_file: str, optional
            Name of the previous IGMGuesses .json file
        smooth: float, optional
            Number of pixels to smooth on (Gaussian smoothing)
        plot_residuals : bool, optional
            Whether to plot residuals
        n_max_tuple : int, optional
            Maximum number of transitions per ion species to consider for plotting display.
        min_strength : float, optional
            Minimum strength for a transition to be considered in the analysis.
            The value should lie between (0,14.7), where 0. means 
            include everything, and 14.7 corresponds to the strength of 
            HI Lya transition assuming solar abundance.
        min_ew : float, optional
            Minimum equivalent width (in AA) of lines to be stored within a components.
            This is useful for not storing extremely weak lines.
        vlim_disp : list of Quantity, optional
            Minimum and maximum velocity limit for the display; e.g. [-500.,500.]*u.km/u.s
        external_model : str
            Name of an external model to load.

        """
        # TODO
        # 1. Fix convolve window size
        # 2. Add COS LSF (?) This seems non-practical

        self.help_message = """
Click on any white region within the velocity plots
for the following keystroke commands to work:

i,o       : zoom in/out x limits
I,O       : zoom in/out x limits (larger re-scale)
Y         : zoom out y limits
y         : guess y limits
t,b       : set y top/bottom limit
l,r       : set left/right x limit
[,]       : pan left/right
C,c       : add/remove column
K,k       : add/remove row
(         : toggle between many/few (15 or 6) panels per page
=,-       : move to next/previous page
f         : move to the first page
Space bar : set redshift from cursor position
^         : set redshift by hand
T         : update available transitions at current redshift from the
            `Strong` LineList
F         : update available transitions at current redshift from the
            full `ISM` LineList
H         : update to HI Lyman series LineList at current redshift
            (type `T` or `F` to get metals back)
U         : restrict to subset of available lines from the parent LineList
            at current redshift
A         : set limits for fitting an absorption component
            from cursor position (need to be pressed twice:
            once for each left and right limits)
a         : Add a new component with the velocity range of the current
1,2       : Modify velocity region of the component (left, right sides)
S         : select an absorption component from cursor position
D         : delete absorption component that is closest to the cursor
            (the cursor has to be in the corresponding velocity window panel
            where the component was defined in the first place)
d         : delete absorption component selected from component widget
N,n       : slightly increase/decrease column density in initial guess
V,v       : slightly increase/decrease b-value in initial guess
<,>       : slightly increase/decrease redshift in initial guess
R         : refit
X,x       : add/remove `bad pixels` (for avoiding using them in subsequent
            VP fitting; works as `A` command, i.e. need to define two limits)
L         : toggle displaying/hiding labels of currently identified lines
M         : toggle displaying/hiding the current absorption model
P         : toggle on/off "colorful" display, where components of different
            reliabilities are plotted in different colors
E         : toggle displaying/hiding the external absorption model
%         : guess a transition and redshift for a given feature at
            the cursor's position
&         : Launch an ExamineSpecWidget to view the spectrum and any absorption systems
?         : print this help message
"""

        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        self.previous_file = previous_file
        if outfil is None:
            self.outfil = 'IGM_model.json'
        else:
            self.outfil = outfil
        if n_max_tuple is None:
            n_max_tuple = 5
        if min_ew is None:
            min_ew = 0.005  # AA
        if min_strength is None:
            min_strength = 0.
        if fwhm is None:
            self.fwhm = 3.
        else:
            self.fwhm = fwhm
        self.plot_residuals = plot_residuals
        self.n_max_tuple = n_max_tuple
        self.min_strength = min_strength
        self.min_ew = min_ew * u.AA
        if vlim_disp is None:
            vlim_disp = [-500.,500.] * u.km/u.s
        self.vlim_disp = vlim_disp

        # Load spectrum
        spec, spec_fil = ltgu.read_spec(ispec)
        # Should do coordinates properly eventually
        self.coord = zero_coord
        # Normalize
        if spec.co_is_set:
            spec.normalize(co=spec.co)
        else:
            raise ValueError("Please provide a spectrum with a continuum estimation. "
                             "You can do this using linetool's `lt_continuumfit` script.")

        # Load external model spectrum
        if external_model is not None:
            self.external_model, _ = ltgu.read_spec(external_model)
        else:
            self.external_model = None

        # make sure there are no nans in uncertainty, which affects the display of residuals
        # removed because now it is handled by XSpectrum1D masked array representations
        # These commented lines should eventually be removed
        # spec.data['sig'][0] = np.where(spec.data['sig'][0].mask, 0, spec.data['sig'][0])

        # These attributes will store good/bad pixels for subsequent Voigt Profile fitting
        # spec.good_pixels = np.zeros(len(spec.wavelength),dtype=int)
        spec.bad_pixels = np.zeros(len(spec.wavelength),dtype=int)

        # Full spectrum model
        self.model = XSpectrum1D.from_tuple(
            (spec.wavelength, np.ones(len(spec.wavelength))))

        # LineList (Grab ISM, Strong and HI as defaults)
        print('Loading built in LineList: ISM, HI and Strong.')
        self.llist = ltgu.set_llist('ISM')
        self.llist['HI'] = LineList('HI')
        self.llist['Strong'] = LineList('Strong')

        # Setup available LineList; this will be the default one
        # which will be updated using a given base Linelist (e.g. 'ISM', 'Strong', 'HI')
        self.llist['available'] = LineList('ISM')

        # Sort the linelists conveniently
        self.llist['ISM'].sortdata(['abundance', 'ion_name', 'rel_strength'], reverse=True)
        self.llist['HI'].sortdata('rel_strength', reverse=True)
        self.llist['Strong'].sortdata(['abundance', 'ion_name', 'rel_strength'], reverse=True)

        # Append the new ones
        self.llist['Lists'].append('HI')
        self.llist['Lists'].append('Strong')
        self.llist['Lists'].append('available')
        if 0: # for H2 testing
            self.llist['H2'] = LineList('H2')
            self.llist['H2'].sortdata(['rel_strength'], reverse=True)
            self.llist['Lists'].append('H2')

        # Define initial redshift
        z = 0.0
        self.llist['z'] = z
        
        # Grab the pieces and tie together
        self.slines_widg = ltgl.SelectedLinesWidget(
            self.llist[self.llist['List']], parent=self, init_select='All')
        self.fiddle_widg = FiddleComponentWidget(parent=self)
        self.comps_widg = ComponentListWidget([], parent=self)
        self.velplot_widg = IGGVelPlotWidget(spec, z, 
            parent=self, llist=self.llist, fwhm=self.fwhm, plot_residuals=self.plot_residuals,
            vmnx=self.vlim_disp, external_model=self.external_model)
        self.wq_widg = ltgsm.WriteQuitWidget(parent=self)


        # Load prevoius file
        if self.previous_file is not None:
            self.read_previous()
        # Connections (buttons are above)
        #self.spec_widg.canvas.mpl_connect('key_press_event', self.on_key)
        #self.abssys_widg.abslist_widget.itemSelectionChanged.connect(
        #    self.on_list_change)

        # Layout
        anly_widg = QtGui.QWidget()
        anly_widg.setMaximumWidth(500)
        anly_widg.setMinimumWidth(250)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.fiddle_widg)
        vbox.addWidget(self.comps_widg)
        vbox.addWidget(self.slines_widg)
        vbox.addWidget(self.wq_widg)
        anly_widg.setLayout(vbox)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.velplot_widg)
        hbox.addWidget(anly_widg)

        self.main_widget.setLayout(hbox)

        # Attempt to initialize
        self.update_available_lines(linelist=self.llist[self.llist['List']])
        self.velplot_widg.init_lines()
        self.velplot_widg.on_draw(rescale=True, fig_clear=True)
        self.slines_widg.selected = self.llist['show_line']
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.slines_widg.on_list_change(self.llist[self.llist['List']])

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

        # Print help message
        print(self.help_message)

    def update_available_lines(self, linelist):
        """Grab the available LineList in the spectrum at the current
        redshift with the current parent linelist (another
        given LineList object)
        """

        print('Updating available transitions from parent LineList: {}, '
              'using n_max_tuple={} and min_strength={:.1f}'.format(linelist.list, self.n_max_tuple, self.min_strength))

        z = self.velplot_widg.z
        wvmin = self.velplot_widg.spec.wvmin
        wvmax = self.velplot_widg.spec.wvmax
        wvlims = (wvmin / (1. + z), wvmax / (1. + z))

        transitions = linelist.available_transitions(
            wvlims, n_max_tuple=self.n_max_tuple, min_strength=self.min_strength)

        if transitions is None:
            print('  There are no transitions available!')
            print('  You will probably need to modify your redshift.')
            return
        elif isinstance(transitions, dict):
            names = [transitions['name']]
        else: # transitions should be a QTable
            names = list(np.array(transitions['name']))
        if len(names) > 0:
            self.llist['available'] = linelist.subset_lines(subset=names, reset_data=True, verbose=True,
                                                            sort_by='as_given')
            self.llist['List'] = 'available'
        print('Done.')

    def on_list_change(self):
        self.update_boxes()

    def create_status_bar(self):
        self.status_text = QtGui.QLabel("IGMGuessesGUI: Please press '?' to display help message in terminal.")
        self.statusBar().addWidget(self.status_text, 1)

    def delete_component(self, component):
        '''Remove component'''
        # Component list
        self.comps_widg.remove_item(component.name)
        # Fiddle query (need to come back to it)
        if component is self.fiddle_widg.component:
            self.fiddle_widg.reset()

        # Mask
        # for line in component.lines:
        #     wvmnx = line.wrest * (1 + component.zcomp) * (1 + component.vlim.value / c_kms)
        #     gdp = np.where((self.velplot_widg.spec.wavelength > wvmnx[0])&
        #         (self.velplot_widg.spec.wavelength < wvmnx[1]))[0]
        #     self.velplot_widg.spec.good_pixels[gdp] = 0

        # Delete
        del component
        # Update
        self.velplot_widg.update_model()
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_slines(self, selected):
        self.llist['show_line'] = selected
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_component(self):
        '''Component attrib was updated. Deal with it'''
        #self.fiddle_widg.component.sync_lines()
        sync_comp_lines(self.fiddle_widg.component)
        mask_comp_lines(self.fiddle_widg.component, min_ew=self.min_ew)

        self.velplot_widg.update_model()
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_compslist(self,component):
        '''Component list was updated'''
        if component is None:
            self.fiddle_widg.reset()
        else:
            self.fiddle_widg.init_component(component)
        #self.velplot_widg.update_model()
        #self.velplot_widg.on_draw(fig_clear=True)

    def read_previous(self):
        ''' Read from a previous guesses file'''
        import json
        # Read the JSON file
        with open(self.previous_file) as data_file:    
            igmg_dict = json.load(data_file)
        # Check FWHM
        if igmg_dict['fwhm'] != self.fwhm:
            raise ValueError('Input FWHMs do not match. Please fix it!')
        # Load bad pixels
        if 'bad_pixels' in igmg_dict.keys():
            bad = igmg_dict['bad_pixels']
            if len(bad) > 0:
                self.velplot_widg.spec.bad_pixels[np.array(bad)] = 1
        # Load good pixels
        # if 'good_pixels' in igmg_dict.keys():
        #     good = igmg_dict['good_pixels']
        # elif 'mask' in igmg_dict.keys(): # old format
        #     good = igmg_dict['mask']
        # if len(good) > 0:
        #     self.velplot_widg.spec.good_pixels[np.array(good)] = 1

        # Check spectra names
        if self.velplot_widg.spec.filename != igmg_dict['spec_file']:
            warnings.warn('Spec file names do not match! Could just be path..')

        # Components
        print('Reading the components from previous file. It may take a while...')
        ntot = len(igmg_dict['cmps'].keys())
        for ii, key in enumerate(igmg_dict['cmps'].keys()):

            if 'lines' in igmg_dict['cmps'][key].keys():
                comp = AbsComponent.from_dict(igmg_dict['cmps'][key], linelist=self.llist['ISM'], coord=self.coord,
                                              chk_sep=False, chk_data=False, chk_vel=False)
                # QtCore.pyqtRemoveInputHook()
                # pdb.set_trace()
                # QtCore.pyqtRestoreInputHook()
                comp_init_attrib(comp)
                comp.init_wrest = igmg_dict['cmps'][key]['wrest']*u.AA
                try:
                    comp.mask_abslines = igmg_dict['cmps'][key]['mask_abslines']
                except KeyError:  # For compatibility
                    warnings.warn("Setting all abslines to 2")
                    comp.mask_abslines = 2*np.ones(len(comp._abslines)).astype(int)
                self.velplot_widg.add_component(comp, update_model=False)

            else:  # for compatibility, should be deprecated
                self.velplot_widg.add_component(
                        igmg_dict['cmps'][key]['wrest']*u.AA,
                        zcomp=igmg_dict['cmps'][key]['zcomp'],
                        vlim=igmg_dict['cmps'][key]['vlim']*u.km/u.s,
                        update_model=False)


            # Name
            self.velplot_widg.current_comp.name = key
            # Set N,b,z
            self.velplot_widg.current_comp.attrib['z'] = igmg_dict['cmps'][key]['zfit']
            self.velplot_widg.current_comp.attrib['b'] = igmg_dict['cmps'][key]['bfit']*u.km/u.s
            self.velplot_widg.current_comp.attrib['logN'] = igmg_dict['cmps'][key]['Nfit']
            try: # This should me removed in the future
                self.velplot_widg.current_comp.attrib['Reliability'] = igmg_dict['cmps'][key]['Reliability']
            except:
                self.velplot_widg.current_comp.attrib['Reliability'] = igmg_dict['cmps'][key]['Quality']  # old version compatibility
            self.velplot_widg.current_comp.comment = igmg_dict['cmps'][key]['Comment']
            # Sync
            sync_comp_lines(self.velplot_widg.current_comp)
            mask_comp_lines(self.velplot_widg.current_comp, min_ew=self.min_ew)

            # print message regarding progress on reading components
            import sys
            progress = int(ii * 100. / ntot)
            sys.stdout.write('Progress: {}%\r'.format(progress))
            sys.stdout.flush()

        # Updates
        self.velplot_widg.update_model()
        self.fiddle_widg.init_component(self.velplot_widg.current_comp)


    def write_out(self):
        """ Write to a JSON file"""
        import json, io
        # Create dict of the components
        out_dict = dict(cmps={},
                        spec_file=self.velplot_widg.spec.filename,
                        fwhm=self.fwhm, bad_pixels=[])

        # Write components out
        # We need a deep copy here because ._abslines will be modified before writting
        # but we want to keep the original ._abslines list in case column density
        # increases and the lines that were masked out because of EW should reappear.
        comps_aux = copy.deepcopy(self.comps_widg.all_comp)
        for kk, comp in enumerate(comps_aux):
            # get rid of masked abslines for writting out to hard drive
            abslines_aux = []
            mask_abslines_aux = []
            for ii, line in enumerate(comp._abslines):
                if comp.mask_abslines[ii] != 0:
                    abslines_aux += [line]
                    mask_abslines_aux += [comp.mask_abslines[ii]]
            comp._abslines = abslines_aux
            comp.mask_abslines = np.array(mask_abslines_aux)

            key = comp.name
            out_dict['cmps'][key] = comp.to_dict()
            out_dict['cmps'][key]['zcomp'] = comp.zcomp
            out_dict['cmps'][key]['zfit'] = comp.attrib['z']
            out_dict['cmps'][key]['Nfit'] = comp.attrib['logN']
            out_dict['cmps'][key]['bfit'] = comp.attrib['b'].value
            out_dict['cmps'][key]['wrest'] = comp.init_wrest.value
            out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['Reliability'] = str(comp.attrib['Reliability'])
            out_dict['cmps'][key]['Comment'] = str(comp.comment)
            out_dict['cmps'][key]['mask_abslines'] = comp.mask_abslines
        # Write bad/good pixels out
        # good_pixels = np.where(self.velplot_widg.spec.good_pixels == 1)[0]
        # if len(good_pixels) > 0:
        #     out_dict['good_pixels'] = list(good_pixels)
        bad_pixels = np.where(self.velplot_widg.spec.bad_pixels == 1)[0]
        if len(bad_pixels) > 0:
            out_dict['bad_pixels'] = list(bad_pixels)

        # JSONify
        gd_dict = ltu.jsonify(out_dict)

        # Write file
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(gd_dict, sort_keys=True, indent=4,
                                       separators=(',', ': '))))
        print('Wrote: {:s}'.format(self.outfil))

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.quit()

    # Quit
    def quit(self):
        self.close()

######################
class IGGVelPlotWidget(QtGui.QWidget):
    """ Widget for a velocity plot with interaction.
          Adapted from VelPlotWidget in spec_guis
        14-Aug-2015 by JXP
    """
    def __init__(self, ispec, z, parent=None, llist=None, norm=True,
                 vmnx=[-500., 500.]*u.km/u.s, fwhm=0., plot_residuals=True, external_model=None):
        '''
        spec = Spectrum1D
        Norm: Bool (False)
          Normalized spectrum?
        abs_sys: AbsSystem
          Absorption system class
        '''
        super(IGGVelPlotWidget, self).__init__(parent)

        # init help message
        self.help_message = parent.help_message

        # Initialize
        self.parent = parent
        spec, spec_fil = ltgu.read_spec(ispec)
        
        self.spec = spec
        self.spec_fil = spec_fil
        self.fwhm = fwhm
        self.z = z
        self.vmnx = vmnx
        self.norm = norm
        # Init
        self.flag_add = False
        self.flag_idlbl = False
        self.flag_mask = False
        self.flag_plotmodel = True
        self.flag_plotextmodel = False
        self.flag_colorful = False

        self.wrest = 0.
        self.avmnx = np.array([0.,0.])*u.km/u.s
        self.model = XSpectrum1D.from_tuple(
            (spec.wavelength, np.ones(len(spec.wavelength))))

        # define external model
        self.external_model = external_model

        #Define arrays for plotting residuals
        self.plot_residuals = plot_residuals
        if self.plot_residuals:
            self.residual_normalization_factor = 0.02/np.median(self.spec.sig)
            self.residual_limit = self.spec.sig * self.residual_normalization_factor
            self.residual = (self.spec.flux - self.model.flux) * self.residual_normalization_factor
            if external_model is not None:
                self.ext_res = (self.spec.flux - self.external_model.flux) * self.residual_normalization_factor

        self.psdict = {} # Dict for spectra plotting
        self.psdict['x_minmax'] = self.vmnx.value # Too much pain to use units with this
        self.psdict['y_minmax'] = [-0.1, 1.1]
        self.psdict['nav'] = ltgu.navigate(0,0,init=True)


        # Status Bar?
        #if not status is None:
        #    self.statusBar = status

        # Line List
        if llist is None:
            self.llist = ltgu.set_llist(['HI 1215', 'HI 1025'])
        else:
            self.llist = llist
        self.llist['z'] = self.z
        # QtCore.pyqtRemoveInputHook()
        # pdb.set_trace()
        # QtCore.pyqtRestoreInputHook()

        # Indexing for line plotting
        self.idx_line = 0

        self.init_lines()
        
        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 150
        self.fig = Figure((8.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        # Sub_plots
        self.sub_xy = [3,2]
        self.subxy_state = 'In'

        self.fig.subplots_adjust(hspace=0.0, wspace=0.1, left=0.04,
                                 right=0.975, top=0.9, bottom=0.07)
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        
        self.setLayout(vbox)

        # Draw on init
        self.on_draw(replot=False) # False only for init

    # Load them up for display
    def init_lines(self):
        wvmin = self.spec.wvmin
        wvmax = self.spec.wvmax
        #
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        wrest = self.llist[self.llist['List']].wrest
        wvobs = (1. + self.z) * wrest
        gdlin = np.where( (wvobs > wvmin) & (wvobs < wvmax))[0]
        self.llist['show_line'] = gdlin
        if len(gdlin) == 0:
            self.llist['show_line'] = np.array([0])

        # Update GUI
        self.parent.slines_widg.selected = self.llist['show_line']
        self.parent.slines_widg.on_list_change(
            self.llist[self.llist['List']])


    # Update model
    def update_model(self):
        if self.parent is None:
            return
        all_comp = self.parent.comps_widg.all_comp # selected_components()
        if len(all_comp) == 0:
            self.model.flux[:] = 1.
            return
        # Setup lines
        wvmin, wvmax = np.min(self.spec.wavelength), np.max(self.spec.wavelength)
        gdlin = []
        for comp in all_comp:
            for ii, line in enumerate(comp._abslines):
                if comp.mask_abslines[ii] == 0: # Do not use these absorption lines
                    continue
                wvobs = (1 + line.z) * line.wrest
                if (wvobs > wvmin) & (wvobs < wvmax):
                    line.attrib['N'] = 10.**line.attrib['logN'] / u.cm**2
                    gdlin.append(line)
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Voigt
        if len(gdlin) > 0:
            self.model = lav.voigt_from_abslines(self.spec.wavelength, gdlin, fwhm=self.fwhm)#,debug=True)
        #Define arrays for plotting residuals
        if self.plot_residuals:
            self.residual_limit = self.spec.sig * self.residual_normalization_factor
            self.residual = (self.spec.flux - self.model.flux) * self.residual_normalization_factor

    # Add a component
    def add_component(self, inp, vlim=None, zcomp=None, update_model=True):
        """Generate a component and fit with Voigt profiles

        Parameters:
        ------------
        inp : wrest or AbsComponent
        vlim : Quantity array, optional
          Velocity limits in km/s
        zcomp : float, optional
        update_model: bool, optional
          Whether to update the model. It is useful to set it to
          False when reading the previous file to increase speed.
        """
        # this is mostly from reading previous IGM_model
        if isinstance(inp, AbsComponent):

            # get rid of lines outside spectral coverage, if any
            # Todo: one may want to reconsider this, specially if more spectral coverage is obtained
            new_abslines = []
            for absline in inp._abslines:
                wobs = absline.wrest * (1 + absline.z)
                if (wobs > self.spec.wvmin) and (wobs < self.spec.wvmax):
                    new_abslines += [absline]
            inp._abslines = new_abslines
            new_comp = inp

            # compatibility with older versions
            try:
                cond = new_comp._abslines[0].analy['wvlim'] == [0,0]*u.AA
            except:
                cond = new_comp._abslines[0].limits.wvlim == [0,0]*u.AA

            if np.sum(cond)>0:
                #sync wvlims
                for aline in new_comp._abslines:
                    vlim = aline.analy['vlim']
                    aline.limits.set(vlim)

        # This is mostly for adding new components from GUI
        else:  # wrest
            # Center z and reset vmin/vmax
            if zcomp is None:
                zmin, zmax = self.z + (1 + self.z) * (self.avmnx.value / c_kms)
                zcomp = 0.5 * (zmin + zmax)
            if vlim is None:
                vlim = self.avmnx - 0.5 * (self.avmnx[1] + self.avmnx[0])

            # Create component from lines available in most complete relevant LineList,
            # in order to avoid missing something
            current_parent_linelist = self.llist[self.llist['List']].list
            if current_parent_linelist in ['ISM', 'Strong', 'EUV', 'HI']:
                linelist_aux = self.llist['ISM']
            elif current_parent_linelist in ['CO', 'H2']:
                linelist_aux = self.llist[current_parent_linelist]

            # add the component
            new_comp = create_component(zcomp, inp, linelist_aux, vlim=vlim, spec=self.spec)

        # Fit
        #print('doing fit for {:g}'.format(wrest))
        if update_model:
            self.fit_component(new_comp)

            # Masking good pixels
            # For Lyman series only mask pixels for fitting 
            # up to Ly-gamma; the rest should be done manually 
            # if wanted
            # if new_comp.lines[0].name.startswith('HI '):
            #     aux_comp_list = new_comp.lines[::-1][:3] #invert order from ISM LineList and truncate for masking
            # else:
            #     aux_comp_list = new_comp.lines

            # for line in aux_comp_list:
                # print('masking {:g}'.format(line.wrest))
                # wvmnx = line.wrest*(1+new_comp.zcomp)*(1 + vlim.value / c_kms)
                # gdp = np.where((self.spec.wavelength>wvmnx[0])&
                #     (self.spec.wavelength<wvmnx[1]))[0]
                # if len(gdp) > 0:
                #     self.spec.good_pixels[gdp] = 1

        # Add to component list and Fiddle
        if self.parent is not None:
            self.parent.comps_widg.add_component(new_comp)
            self.parent.fiddle_widg.init_component(new_comp)

        # Update model
        self.current_comp = new_comp
        if update_model:
            self.update_model()

    def fit_component(self, component):
        '''Fit the component and save values'''
        from astropy.modeling import fitting
        # Generate Fit line
        fit_line = AbsLine(component.init_wrest, linelist=self.llist[self.llist['List']], z=component.zcomp)
        fit_line.limits.set(component.vlim)
        fit_line.analy['spec'] = self.spec

        fit_line.measure_aodm(normalize=False)  # Already normalized

        # Guesses
        fmin = np.argmin(self.spec.flux[fit_line.analy['pix']])
        zguess = self.spec.wavelength[fit_line.analy['pix'][fmin]]/component.init_wrest - 1.
        bguess = 0.5 * (component.vlim[1] - component.vlim[0])
        Nguess = np.log10(fit_line.attrib['N'].to('cm**-2').value)
        # Voigt model
        fitvoigt = lav.single_voigt_model(logN=Nguess,b=bguess.value,
                                z=zguess, wrest=component.init_wrest.value,
                                gamma=fit_line.data['gamma'].value, 
                                f=fit_line.data['f'], fwhm=self.fwhm)
        # Restrict parameter space
        fitvoigt.logN.min = 12.
        fitvoigt.b.min = 10.
        fitvoigt.z.min = component.zcomp + component.vlim[0].value * (1 + component.zcomp) / c_kms
        fitvoigt.z.max = component.zcomp + component.vlim[1].value * (1 + component.zcomp) / c_kms

        # Fit
        fitter = fitting.LevMarLSQFitter()
        parm = fitter(fitvoigt,self.spec.wavelength[fit_line.analy['pix']],
            self.spec.flux[fit_line.analy['pix']].value)

        # Save and sync
        component.attrib['logN'] = parm.logN.value
        component.attrib['N'] = 10**parm.logN.value / u.cm**2
        component.attrib['z'] = parm.z.value
        component.attrib['b'] = parm.b.value * u.km/u.s
        #component.sync_lines()
        sync_comp_lines(component)
        mask_comp_lines(component, min_ew=self.parent.min_ew)

    def out_of_bounds(self,coord):
        '''Check for out of bounds
        '''
        # Default is x
        if ((coord < np.min(self.spec.wavelength))
            or (coord > np.max(self.spec.wavelength))):
            print('Out of bounds!')
            return True
        else:
            return False

    # Key stroke 
    def on_key(self,event):
        # Init
        rescale = True
        fig_clear = False
        wrest = None
        flg = 1
        sv_idx = self.idx_line

        # add/remove rows/columns
        if event.key == 'k':
            self.sub_xy[0] = max(0, self.sub_xy[0]-1)
        if event.key == 'K':
            self.sub_xy[0] = self.sub_xy[0]+1
        if event.key == 'c':
            self.sub_xy[1] = max(1, self.sub_xy[1]-1)
        if event.key == 'C':
            self.sub_xy[1] = max(1, self.sub_xy[1]+1)
        # toggle between many/few panels
        if event.key == '(':
            if self.subxy_state == 'Out':
                self.sub_xy = [3,2]
                self.subxy_state = 'In'
            else:
                self.sub_xy = [5,3]
                self.subxy_state = 'Out'

        ## NAVIGATING
        if event.key in self.psdict['nav']:
            flg = ltgu.navigate(self.psdict,event)
        if event.key == '-':  # previous page
            self.idx_line = max(0, self.idx_line-self.sub_xy[0]*self.sub_xy[1]) # Min=0
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == '=':  # next page
            self.idx_line = min(len(self.llist['show_line'])-self.sub_xy[0]*self.sub_xy[1],
                                self.idx_line + self.sub_xy[0]*self.sub_xy[1])
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == 'f':  # go to the first page
            self.idx_line = 0
            print('Edge of list')

        # Find line
        try:
            wrest = event.inaxes.get_gid()
        except AttributeError:
            return
        else:
            wvobs = wrest*(1+self.z)
            pass

        ## Fiddle with a Component
        if event.key in ['N','n','v','V','<','>','R','1','2']:
            if self.parent.fiddle_widg.component is None:
                print('Need to generate a component first!')
                return
            comp = self.parent.fiddle_widg.component
            if event.key == 'N':
                comp.attrib['logN'] += 0.05
            elif event.key == 'n':
                comp.attrib['logN'] -= 0.05
            elif event.key == 'v':
                self.parent.fiddle_widg.component.attrib['b'] -= 5*u.km/u.s
            elif event.key == 'V':
                self.parent.fiddle_widg.component.attrib['b'] += 5*u.km/u.s
            elif event.key == '<':
                self.parent.fiddle_widg.component.attrib['z'] -= 4e-5  # should be a fraction of pixel size
            elif event.key == '>':
                self.parent.fiddle_widg.component.attrib['z'] += 4e-5

            elif event.key == 'R': # Refit
                self.fit_component(self.parent.fiddle_widg.component)
            elif event.key == '1':
                # dvz_kms = c_kms * (self.z - comp.zcomp) / (1 + self.z)
                dvz_kms = -1*ltu.dv_from_z(comp.zcomp, self.z)
                comp.vlim[0] = (event.xdata*u.km/u.s + dvz_kms)
                sync_comp_lines(comp, only_lims=True)
            elif event.key == '2':
                # dvz_kms = c_kms * (self.z - comp.zcomp) / (1 + self.z)
                dvz_kms = -1*ltu.dv_from_z(comp.zcomp, self.z)
                comp.vlim[1] = (event.xdata*u.km/u.s + dvz_kms)
                sync_comp_lines(comp, only_lims=True)
            # Updates (this captures them all and redraws)
            self.parent.fiddle_widg.update_component()

        ## Grab/Delete a component
        if event.key in ['D','S','d']:
            # Delete selected component
            if event.key == 'd':
                self.parent.delete_component(self.parent.fiddle_widg.component)
                return

            components = self.parent.comps_widg.all_comp
            iwrest = np.array([comp.init_wrest.value for comp in components])*u.AA
            mtc = np.where(wrest == iwrest)[0]
            if len(mtc) == 0:
                return
            # QtCore.pyqtRemoveInputHook()
            # pdb.set_trace()
            # QtCore.pyqtRestoreInputHook()
            dvz = np.array([c_kms * (self.z - components[mt].zcomp) / (1+self.z) for mt in mtc])
            # Find minimum
            mindvz = np.argmin(np.abs(dvz+event.xdata))
            if event.key == 'S':
                self.parent.fiddle_widg.init_component(components[mtc[mindvz]])
                # update the selection in Component Widget:
                idx = mtc[mindvz]
                self.parent.comps_widg.complist_widget.setCurrentRow(idx+1) # +1 because None counts
            elif event.key == 'D': # Delete nearest component to cursor
                self.parent.delete_component(components[mtc[mindvz]])

        ## Reset z
        if event.key == ' ': # space to move redshift
            self.z = self.z + event.xdata * (1 + self.z) / c_kms
            # Drawing
            self.psdict['x_minmax'] = self.vmnx.value

        ## Launch ExamineSpecWidget
        if event.key == '&':
            # Create absorption systems
            comps = self.parent.comps_widg.all_comp
            if len(comps) > 0:
                abs_sys = ltiu.build_systems_from_components(comps)
            else:
                abs_sys = None
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            # Launch Gui
            exspecgui = ltspw.ExSpecDialog(self.spec, abs_sys=abs_sys)
            exspecgui.exec_()

        ## Set redshift
        if event.key == '^':
            zgui = ltgsm.AnsBox('Enter redshift:',float)
            zgui.exec_()
            self.z = zgui.value
            self.psdict['x_minmax'] = self.vmnx.value

        # Choose line
        if event.key == '%':
            # GUI
            aux_table = self.llist[self.llist['List']]._data  # is a table of lines
            z_aux = wvobs/aux_table['wrest'] - 1.
            aux_table['redshift'] = z_aux.value  # adding aux_z to the aux_table
            self.select_line_widg = ltgl.SelectLineWidget(aux_table)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return

            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            #
            self.z = (wvobs/wrest - 1.).value
            #self.statusBar().showMessage('z = {:f}'.format(z))
            self.init_lines()

        # Select the base LineList from keystroke
        if event.key == 'H':  # Load HI
            self.llist['List'] = 'HI'
            # self.parent.update_available_lines(linelist=self.llist['HI'])
            self.idx_line = 0
            self.init_lines()
            s = "current parent LineList set to 'HI'."
            print(s)
            self.parent.statusBar().showMessage('IGMGuessesGUI: '+ s)
        if event.key == 'T':  # Load Strong
            self.llist['List'] = 'Strong'
            # self.parent.update_available_lines(linelist=self.llist['Strong'])
            self.idx_line = 0
            self.init_lines()
            s = "current parent LineList set to 'Strong'."
            print(s)
            self.parent.statusBar().showMessage('IGMGuessesGUI: '+ s)
        if event.key == '9':  # Load H2
            self.llist['List'] = 'H2'
            # self.parent.update_available_lines(linelist=self.llist['Strong'])
            self.idx_line = 0
            self.init_lines()
            s = "current parent LineList set to 'H2'."
            print(s)
            self.parent.statusBar().showMessage('IGMGuessesGUI: '+ s)



        if event.key == 'F':  # Load Full ISM
            self.llist['List'] = 'ISM'
            # self.parent.update_available_lines(linelist=self.llist['ISM'])
            self.idx_line = 0
            self.init_lines()
            s = "current parent LineList set to 'ISM'."
            print(s)
            self.parent.statusBar().showMessage('IGMGuessesGUI: '+ s)
        if event.key == 'U':  # Update subset of the available lines
            current_linelist = self.llist[self.llist['List']]
            parent_linelist = self.llist[current_linelist.list]
            # QtCore.pyqtRemoveInputHook()
            # pdb.set_trace()
            # QtCore.pyqtRestoreInputHook()
            self.parent.update_available_lines(linelist=parent_linelist)
            self.idx_line = 0
            self.init_lines()
            s = "current LineList set to a subset of '{}' at z={:.3f}.".format(current_linelist.list, self.z)
            self.parent.statusBar().showMessage('IGMGuessesGUI: '+ s)

        ## Add component
        if event.key == 'A': # Add to lines
            if self.out_of_bounds(wvobs * (1 + event.xdata / c_kms)):
                return
            if self.flag_add is False:
                self.vtmp = event.xdata
                self.flag_add = True
                self.wrest = wrest
            else:
                self.avmnx = np.array([np.minimum(self.vtmp,event.xdata),
                    np.maximum(self.vtmp,event.xdata)])*u.km/u.s
                self.add_component(wrest)
                # Reset
                self.flag_add = False
                self.wrest = 0.
        if event.key == 'a': # Add to lines using current vlim
            if self.parent.fiddle_widg.component is None:
                print('Need to generate a component first!')
                return
            comp = self.parent.fiddle_widg.component
            dvz_kms = c_kms * (self.z - comp.zcomp) / (1 + self.z)
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.wrest = wrest
            self.avmnx = comp.vlim - dvz_kms*u.km/u.s
            self.add_component(wrest, zcomp=comp.zcomp, vlim=comp.vlim)
            self.wrest = 0.

        # Fiddle with analysis pixel mask
        if event.key in ['x','X']:
            # x = Delete mask
            # X = Add to mask
            if self.flag_mask is False:
                self.wrest = wrest
                self.wtmp = wvobs * (1 + event.xdata / c_kms)
                self.vtmp = event.xdata
                self.flag_mask = True
            else:
                wtmp2 = wvobs * (1 + event.xdata / c_kms)
                twvmnx = [np.minimum(self.wtmp,wtmp2), np.maximum(self.wtmp,wtmp2)]
                # Modify mask
                mskp = np.where((self.spec.wavelength>twvmnx[0])&
                    (self.spec.wavelength<twvmnx[1]))[0]
                #print(twvmnx,len(mskp))
                if event.key == 'x':
                    self.spec.bad_pixels[mskp] = 0
                elif event.key == 'X':
                    self.spec.bad_pixels[mskp] = 1
                # Reset
                self.flag_mask = False
                self.wrest = 0.

        # Plot component labels?
        if event.key == 'L':  # Toggle ID lines
            self.flag_idlbl = not self.flag_idlbl

        # Plot model?
        if event.key == 'M':  # Toggle plot model
            self.flag_plotmodel = not self.flag_plotmodel

        # Colorful mode?
        if event.key == 'P':  # Toggle colorful mode
            self.flag_colorful = not self.flag_colorful

        # Plot external model?
        if event.key == 'E':  # Toggle plot external model
            if self.external_model is None:
                s = "No external model loaded. Please restart IGMGuesses with --external_model input."
                print(s)
                self.parent.statusBar().showMessage('IGMGuessesGUI: ' + s)
            else:
                self.flag_plotextmodel = not self.flag_plotextmodel

        # print component info
        if event.key == '@':  #  for develop; this will be deleted in future.
            print('Computing blends between components, it may take a while...\n')
            # Write joebvp files; for developing.
            blending_info(self.parent.comps_widg.all_comp, self.spec_fil)  # this one writes them internally. Should be cleaned.

            joebvp_output = self.parent.outfil.replace('.json', '.joebvp')
            from_igmguesses_to_joebvp(self.parent.outfil, joebvp_output)
            print('Wrote: {:s}'.format(joebvp_output))

        """
        # AODM plot
        if event.key == ':':  # 
            # Grab good lines
            from xastropy.xguis import spec_guis as xsgui
            gdl = [iline.wrest for iline in self.abs_sys.lines
                if iline.analy['do_analysis'] > 0]
            # Launch AODM
            if len(gdl) > 0:
                gui = xsgui.XAODMGui(self.spec, self.z, gdl, vmnx=self.vmnx, norm=self.norm)
                gui.exec_()
            else:
                print('VelPlot.AODM: No good lines to plot')
        """

        if event.key == '?':
            print(self.help_message)

        #if wrest is not None: # Single window
        #    flg = 3
        if event.key in ['c','C','k','K','W','!', '@', '=', '-', 'X', ' ','R']: # Redraw all
            flg = 1
        if event.key in ['Y']:
            rescale = False
        if event.key in ['c','C','k','K', 'R', '(', 'H', 'F', 'T', 'U']:
            fig_clear = True

        if flg==1: # Default is to redraw
            self.on_draw(rescale=rescale, fig_clear=fig_clear)
        elif flg==2: # Layer (no clear)
            self.on_draw(replot=False, rescale=rescale)
        elif flg==3: # Layer (no clear)
            self.on_draw(in_wrest=wrest, rescale=rescale)

    # Click of main mouse button
    def on_click(self,event):
        try:
            print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
                event.button, event.x, event.y, event.xdata, event.ydata))
        except ValueError:
            return
        if event.button == 1: # Draw line
            # remove this ugly green line
            # self.ax.plot( [event.xdata,event.xdata], self.psdict['y_minmax'], ':', color='green')
            # self.on_draw(replot=False)
    
            # Print values
            try:
                self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
            except AttributeError:
                return

    def on_draw(self, replot=True, in_wrest=None, rescale=True, fig_clear=False):
        """ Redraws the figure
        """

        if replot is True:
            if fig_clear:
                self.fig.clf()
            # Title by redshift
            self.fig.suptitle('{}\nz={:.5f}'.format(self.spec.filename,self.z), fontsize='medium', ha='center')

            # Components
            components = self.parent.comps_widg.all_comp 
            iwrest = np.array([comp.init_wrest.value for comp in components])*u.AA
            # Loop on windows
            all_idx = self.llist['show_line']
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()

            # Labels and colors for reliability
            # todo: potential for speed up here... e.g. by giving more attributes to components
            if (self.flag_idlbl) or (self.flag_colorful):
                line_wvobs = []
                line_lbl = []
                line_color = []
                line_wvobs_lims = []
                for comp in components:
                    if comp.attrib['Reliability'] == 'None':
                        la = ''
                    else:
                        la = comp.attrib['Reliability']
                    for ii, line in enumerate(comp._abslines):
                        if comp.mask_abslines[ii] == 0:  # do not plot these masked out lines
                            continue
                        zline = line.z
                        wrest = line.wrest.value
                        line_wvobs.append(wrest * (1 + zline))
                        line_lbl.append(line.name+',{:.3f}{:s}'.format(line.z, la))
                        line_wvobs_lims.append(line.limits.wvlim)

                        if la == 'a':
                            line_color.append(COLOR_RELIABLE)
                        elif la == 'b':
                            line_color.append(COLOR_PROBABLE)
                        elif la == 'c':
                            line_color.append(COLOR_UNCERTAIN)
                        else:
                            line_color.append(COLOR_MODEL)
                # save them as numpy arrays
                line_wvobs = np.array(line_wvobs)*u.AA
                line_lbl = np.array(line_lbl)
                line_color = np.array(line_color)
                line_wvobs_lims = np.array(line_wvobs_lims)*u.AA
                # QtCore.pyqtRemoveInputHook()
                # pdb.set_trace()
                # QtCore.pyqtRestoreInputHook()

            # Subplots
            nplt = self.sub_xy[0]*self.sub_xy[1]
            if len(all_idx) <= nplt:
                self.idx_line = 0
            subp = np.arange(nplt) + 1
            subp_idx = np.hstack(subp.reshape(self.sub_xy[0],self.sub_xy[1]).T)
            #print('idx_l={:d}, nplt={:d}, lall={:d}'.format(self.idx_line,nplt,len(all_idx)))
            
            # try different color per ion species, and grey for model, using global
            # variables COLOR_MODEL (str) and COLORS (list of str)
            color_ind = 0

            # loop over individual velplot axes
            for jj in range(min(nplt, len(all_idx))):
                try:
                    idx = all_idx[jj+self.idx_line]
                except IndexError:
                    continue # Likely too few lines
                #print('jj={:d}, idx={:d}'.format(jj,idx))
                # Grab line
                wrest = self.llist[self.llist['List']].wrest[idx]
                kwrest = wrest.value # For the Dict
                
                #define colors for visually grouping same species together
                if jj > 0:
                    name_aux = self.llist[self.llist['List']].name[idx].split(' ')[0]
                    name_aux2 = self.llist[self.llist['List']].name[idx-1].split(' ')[0]
                    if name_aux != name_aux2:
                        color_ind += 1
                color = COLORS[color_ind % len(COLORS)]

                # Single window?
                #if in_wrest is not None:
                #    if np.abs(wrest-in_wrest) > (1e-3*u.AA):
                #        continue
                # Generate plot
                self.ax = self.fig.add_subplot(self.sub_xy[0],self.sub_xy[1], subp_idx[jj])
                self.ax.clear()        

                # GID for referencing
                self.ax.set_gid(wrest)

                # Zero velocity line
                self.ax.plot( [0., 0.], [-1e9, 1e9], ':', color='gray')
                # Unity flux level line
                self.ax.plot( [-1e9, 1e9], [1, 1], ':', color='b', lw=0.5)

                # Velocity
                wvobs = (1 + self.z) * wrest
                wvmnx = wvobs * (1 + np.array(self.psdict['x_minmax']) / c_kms)
                # velo = (self.spec.wavelength/wvobs - 1.) * c_kms * u.km/u.s
                velo = ltu.dv_from_z(self.spec.wavelength/wrest - 1., self.z)

                # Errors
                if self.plot_residuals:  # this could be always true
                    self.ax.plot(velo, self.residual_limit, 'k-', drawstyle='steps-mid', lw=0.5)
                    self.ax.plot(velo, -self.residual_limit, 'k-', drawstyle='steps-mid', lw=0.5)

                # Plot external model
                if self.flag_plotextmodel:  # can only be True if self.external_model is not None.
                    velo_external = ltu.dv_from_z(self.external_model.wavelength / wrest - 1., self.z)
                    self.ax.plot(velo_external, self.external_model.flux, '-', color=COLOR_EXTMODEL, lw=1.0)
                    if self.plot_residuals:
                        self.ax.plot(velo_external, self.ext_res, '.', color=COLOR_EXTMODEL, ms=2)

                # Plot spectrum and model
                self.ax.plot(velo, self.spec.flux, '-', color=color, drawstyle='steps-mid', lw=0.5)
                if self.flag_plotmodel:
                    self.ax.plot(velo, self.model.flux, '-', color=COLOR_MODEL, lw=0.5)
                    if self.plot_residuals:
                        self.ax.plot(velo, self.residual, '.', color=COLOR_MODEL, ms=2)


                # Labels for axes and ion names
                if (((jj+1) % self.sub_xy[0]) == 0) or ((jj+1) == len(all_idx)):
                    self.ax.set_xlabel('Relative Velocity (km/s)')
                else:
                    self.ax.get_xaxis().set_ticks([])
                lbl = self.llist[self.llist['List']].name[idx]
                self.ax.text(0.01, 0.15, lbl, color=color, transform=self.ax.transAxes,
                             size='x-small', ha='left', va='center', backgroundcolor='w',
                             bbox={'pad':0, 'edgecolor':'none', 'facecolor':'w'})

                # labels for individual components
                if (self.flag_idlbl) or (self.flag_colorful):
                    # Any lines inside?
                    mtw = np.where((line_wvobs > wvmnx[0]) & (line_wvobs < wvmnx[1]))[0]
                    for imt in mtw:
                        # v = c_kms * (line_wvobs[imt]/wvobs - 1)
                        v = ltu.dv_from_z(line_wvobs[imt]/wrest - 1., self.z)
                        if self.flag_colorful:
                            color_label = line_color[imt]
                        else:
                            color_label = COLOR_MODEL

                        if self.flag_idlbl:
                            self.ax.text(v.value, 0.5, line_lbl[imt], color=color_label, backgroundcolor='w',
                                bbox={'pad':0,'edgecolor':'none', 'facecolor':'w'}, size='xx-small',
                                    rotation=90.,ha='center',va='center')

                        if (self.flag_plotmodel) and (self.flag_colorful):
                            # plotting absline with color better done in wobs -> velo space
                            # dvmin = c_kms * ((line_wvobs_lims[imt][0] - line_wvobs[imt]) / wvobs)
                            # dvmax = c_kms * ((line_wvobs_lims[imt][1] - line_wvobs[imt]) / wvobs)
                            vmin, vmax = ltu.dv_from_z(line_wvobs_lims[imt]/wrest - 1., self.z)
                            # QtCore.pyqtRemoveInputHook()
                            # pdb.set_trace()
                            # QtCore.pyqtRestoreInputHook()
                            # cond = (velo >= v - dvmin) & (velo <= v + dvmax)
                            cond = (velo >= vmin) & (velo <= vmax)
                            self.ax.plot(velo[cond], self.model.flux[cond], '-', color=color_label, lw=1)

                # Plot good pixels
                # if np.sum(self.spec.good_pixels) > 0.:
                #     gdp = self.spec.good_pixels == 1
                #     if len(gdp) > 0:
                #         self.ax.scatter(velo[gdp],self.spec.flux[gdp],
                #             marker='o',color=color, s=3.,alpha=0.5)

                # Plot bad pixels
                if np.sum(self.spec.bad_pixels) > 0.:
                    bad = self.spec.bad_pixels == 1
                    if len(bad) > 0:
                        self.ax.scatter(velo[bad],self.spec.flux[bad],
                            marker='x',color=color, s=20., alpha=0.5, lw=0.5)

                # Reset window limits
                self.ax.set_ylim(self.psdict['y_minmax'])
                self.ax.set_xlim(self.psdict['x_minmax'])

                # Add line?
                if self.wrest == wrest:
                    self.ax.plot([self.vtmp]*2,self.psdict['y_minmax'], '--',
                        color='red')

                # Components
                mtc = np.where(wrest == iwrest)[0]
                if len(mtc) > 0:
                    #QtCore.pyqtRemoveInputHook()
                    #pdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    for mt in mtc:
                        comp = components[mt]
                        #QtCore.pyqtRemoveInputHook()
                        #pdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        # dvz_kms = c_kms * (self.z - comp.zcomp) / (1 + self.z)
                        dvz_kms = ltu.dv_from_z(comp.zcomp, self.z).to('km/s').value
                        if dvz_kms < 1.5*np.max(np.abs(self.psdict['x_minmax'])):  # 1.5 to make sure the plot is still there at the edges of the axis
                            if comp is self.parent.fiddle_widg.component:
                                lw = 1.5
                            else:
                                lw = 1.
                            # Plot
                            for vlim in comp.vlim:
                                self.ax.plot([vlim.value+dvz_kms]*2, self.psdict['y_minmax'],
                                    '--', color='r',linewidth=lw)
                            self.ax.plot([dvz_kms]*2,[1.0,1.1],
                                '-', color='grey',linewidth=lw)

                # Fonts
                set_fontsize(self.ax,6.)
        # Draw
        self.canvas.draw()
############        
class FiddleComponentWidget(QtGui.QWidget):
    ''' Widget to fiddle with a given component
    '''
    def __init__(self, component=None, parent=None):
        '''
        '''
        super(FiddleComponentWidget, self).__init__(parent)

        self.parent = parent
        #if not status is None:
        #    self.statusBar = status
        self.label = QtGui.QLabel('Component:',self)
        self.zwidget = ltgsm.EditBox(-1., 'zc=', '{:0.5f}')
        self.Nwidget = ltgsm.EditBox(-1., 'Nc=', '{:0.2f}')
        self.bwidget = ltgsm.EditBox(-1., 'bc=', '{:0.1f}')

        self.ddlbl = QtGui.QLabel('Reliability')
        self.ddlist = QtGui.QComboBox(self)
        self.ddlist.addItem('None')
        self.ddlist.addItem('a')
        self.ddlist.addItem('b')
        self.ddlist.addItem('c')
        self.Cwidget = ltgsm.EditBox('None', 'Comment=', '{:s}')

        # Init further
        if component is not None:
            self.init_component(component)
        else:
            self.component = component

        # Connect
        self.ddlist.activated[str].connect(self.setReliability)
        self.connect(self.Nwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.zwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.bwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.Cwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)

        # Layout
        zNbwidg = QtGui.QWidget()
        hbox2 = QtGui.QHBoxLayout()
        hbox2.addWidget(self.zwidget)
        hbox2.addWidget(self.Nwidget)
        hbox2.addWidget(self.bwidget)
        zNbwidg.setLayout(hbox2)

        ddwidg = QtGui.QWidget()
        vbox1 = QtGui.QVBoxLayout()
        vbox1.addWidget(self.ddlbl)
        vbox1.addWidget(self.ddlist)
        ddwidg.setLayout(vbox1)

        commwidg = QtGui.QWidget()
        hbox3 = QtGui.QHBoxLayout()
        hbox3.addWidget(ddwidg)
        hbox3.addWidget(self.Cwidget)
        commwidg.setLayout(hbox3)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.label)
        vbox.addWidget(zNbwidg)
        vbox.addWidget(commwidg)
        self.setLayout(vbox)

    def init_component(self,component):
        '''Setup Widget for the input component'''
        self.component = component
        # Values
        self.Nwidget.set_text(self.component.attrib['logN'])
        self.zwidget.set_text(self.component.attrib['z'])
        self.bwidget.set_text(self.component.attrib['b'].value)
        self.Cwidget.set_text(self.component.comment)
        # Reliability
        idx = self.ddlist.findText(self.component.attrib['Reliability'])
        self.ddlist.setCurrentIndex(idx)
        # Label
        self.set_label()

    def setReliability(self, text):
        if self.component is not None:
            self.component.attrib['Reliability'] = text

    def reset(self):
        #
        self.component = None
        #  Values
        self.Nwidget.set_text(-1.)
        self.zwidget.set_text(-1.)
        self.bwidget.set_text(-1.)
        self.Cwidget.set_text('None')
        idx = self.ddlist.findText('None')
        self.ddlist.setCurrentIndex(idx)
        # Label
        self.set_label()

    def update_component(self):
        '''Values have changed'''
        self.Nwidget.set_text(self.component.attrib['logN'])
        self.zwidget.set_text(self.component.attrib['z'])
        self.bwidget.set_text(self.component.attrib['b'].value)
        self.Cwidget.set_text(self.component.comment)
        if self.parent is not None:
            self.parent.updated_component()

    def set_label(self):
        '''Sets the label for the Widget'''
        if self.component is not None:
            self.label.setText('Component: {:s}'.format(self.component.name))            
        else:
            self.label.setText('Component:')

    def setbzN(self):
        '''Set the component column density or redshift from the boxes'''
        if self.component is None:
            print('Need to generate a component first!')
        else:
            # Grab values
            self.component.attrib['logN'] = (float(self.Nwidget.box.text()))
            self.component.attrib['z'] = (float(self.zwidget.box.text()))
            self.component.attrib['b'] = (float(self.bwidget.box.text()))*u.km/u.s
            self.component.comment = str(self.Cwidget.box.text())
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            # Update beyond
            if self.parent is not None:
                self.parent.updated_component()

# #####
class ComponentListWidget(QtGui.QWidget):
    ''' Widget to organize components on a sightline

    Parameters:
    -----------
    components: List
      List of components

    16-Dec-2014 by JXP
    '''
    def __init__(self, components, parent=None, no_buttons=False):
        '''
        only_one: bool, optional
          Restrict to one selection at a time? [False]
        no_buttons: bool, optional
          Eliminate Refine/Reload buttons?
        '''
        super(ComponentListWidget, self).__init__(parent)

        self.parent = parent

        #if not status is None:
        #    self.statusBar = status
        self.all_comp = components  # Actual components

        list_label = QtGui.QLabel('Components:')
        self.complist_widget = QtGui.QListWidget(self) 
        #self.complist_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.complist_widget.addItem('None')
        #self.abslist_widget.addItem('Test')

        # Lists
        self.items = []     # Selected
        self.all_items = [] # Names

        self.complist_widget.setCurrentRow(0)
        self.complist_widget.itemSelectionChanged.connect(self.on_list_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(list_label)
        vbox.addWidget(self.complist_widget)
        self.setLayout(vbox)

    def on_list_change(self):
        '''
        Changed an item in the list
        '''
        item = self.complist_widget.selectedItems()
        try:
            txt = item[0].text()
        except:
            # This seems to happen when the user forgets
            # to click on the white regions in the velocity
            # plots, before doing a command. Seems harmless.
            return
            # QtCore.pyqtRemoveInputHook()
            # pdb.set_trace()
            # QtCore.pyqtRestoreInputHook()
        if txt == 'None':
            if self.parent is not None:
                self.parent.updated_compslist(None)
        else:
            ii = self.all_items.index(txt)
            if self.parent is not None:
                self.parent.updated_compslist(self.all_comp[ii])

        '''
        items = self.complist_widget.selectedItems()
        # Empty the list
        #self.abs_sys = []
        if len(self.abs_sys) > 0:
            for ii in range(len(self.abs_sys)-1,-1,-1):
                self.abs_sys.pop(ii)
        # Load up abs_sys (as need be)
        new_items = []
        for item in items:
            txt = item.text()
            # Dummy
            if txt == 'None':
                continue
            print('Including {:s} in the list'.format(txt))
            # Using LLS for now.  Might change to generic
            new_items.append(txt)
            ii = self.all_items.index(txt)
            self.abs_sys.append(self.all_abssys[ii])

        # Pass back
        self.items = new_items
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        '''
    '''
    def selected_components(self):
        items = self.complist_widget.selectedItems()
        selc = []
        for item in items:
            txt = item.text()
            if txt == 'None':
                continue
            ii = self.all_items.index(txt)
            selc.append(self.all_comp[ii])
        # Return
        return selc
    '''

    def add_component(self,component):
        self.all_comp.append( component )
        self.add_item(component.name)

    def add_item(self,comp_name):
        #
        self.all_items.append(comp_name) 
        self.complist_widget.addItem(comp_name)
        self.complist_widget.item(len(self.all_items)).setSelected(True)

    def remove_item(self,comp_name):
        # Delete
        idx = self.all_items.index(comp_name)
        del self.all_items[idx]
        self.all_comp.pop(idx)

        tmp = self.complist_widget.takeItem(idx+1)  # +1 because of None counts in widget list but no in component list

        self.complist_widget.item(len(self.all_items)).setSelected(True)
        self.on_list_change()


def create_component(z, wrest, linelist, vlim=[-300.,300]*u.km/u.s,
                     spec=None):
    # Transitions
    all_trans = linelist.all_transitions(wrest)
    if isinstance(all_trans, dict):
        all_trans = [all_trans]
    abslines = []
    if spec is not None:
        wvmin, wvmax = spec.wvmin, spec.wvmax
    else:
        wvmin, wvmax = 0.*u.AA, 1e9*u.AA
    for trans in all_trans:
        # Restrict to those with spectral coverage!
        if (trans['wrest'] * (1 + z) < wvmin) or (trans['wrest'] * (1 + z) > wvmax):
            continue
        aline = AbsLine(trans['wrest'],  linelist=linelist, z=z)
        aline.limits.set(vlim)
        abslines.append(aline)
    if linelist.list != 'H2':
        if (abslines[0].data['Ej'].value > 0.):
            stars = '*'*(len(abslines[0].name.split('*'))-1)
        else:
            stars = None
    else:
        stars = None
    # AbsComponent
    comp = AbsComponent.from_abslines(abslines, stars=stars)
    # Init_wrest
    comp.init_wrest = wrest
    # Attributes
    comp_init_attrib(comp)

    # QtCore.pyqtRemoveInputHook()
    # pdb.set_trace()
    # QtCore.pyqtRestoreInputHook()

    # Mask abslines within a component
    # 0: Do not use
    # 1: Use for display only
    # 2: Use for subsequent VP fitting
    # 3: Out of wavelength range
    comp.mask_abslines = 2*np.ones(len(comp._abslines)).astype(int)

    # Component name
    comp.name = 'z{:.5f}_{:s}'.format(
            comp.zcomp, comp._abslines[0].data['name'].split(' ')[0])
    return comp


def comp_init_attrib(comp):
    # Attributes
    comp.attrib = {'N': 0./u.cm**2, 'Nsig': 0./u.cm**2, 'flagN': 0,  # Column
               'logN': 0., 'sig_logN': 0.,
               'b': 0.*u.km/u.s, 'bsig': 0.*u.km/u.s,  # Doppler
               'z': comp.zcomp, 'zsig': 0.,
               'Reliability': 'None'}


def sync_comp_lines(comp, only_lims=False):
    """Synchronize attributes of the lines and updates
    """
    if only_lims:
        # Need to update zlim *not* vlim
        zlim_comp = ltu.z_from_dv(comp.vlim, comp.zcomp, rel=False)
    for line in comp._abslines:
        if only_lims:
            line.limits.set(zlim_comp.tolist())
        else:
            line.attrib['logN'] = comp.attrib['logN']
            line.attrib['b'] = comp.attrib['b']
            line.setz(comp.attrib['z'])  # vlim is updated because _zlim is fixed


def mask_comp_lines(comp, min_ew = 0.003*u.AA, verbose=False):
    """ Mask out lines that are weaker than
    equivalent width threshold."""

    for ii, line in enumerate(comp._abslines):
        # Estimate equivalent width assuming optically thin line
        # This is ok because we want to mask out weak lines
        fosc = line.data['f']
        wrest = line.data['wrest']
        ew = fosc * wrest**2 * 10**line.attrib['logN'] / u.cm / (1.13 * 10**12 )  # eq. 9.15 Draine 2011
        if ew < min_ew:
            if verbose:
                print('Comp {}: AbsLine {} has estimated EW={:.4f} A < {} A; '
                      'masking out.'.format(comp.name, line.name, ew.to('AA').value, min_ew.to('AA').value, comp.name))
            comp.mask_abslines[ii] = 0
        else:  # line is strong enough, do not mask out
            if comp.mask_abslines[ii] == 0: # if it already masked out, unmask it
                comp.mask_abslines[ii] = 2


    # Sanity check
    if np.sum(comp.mask_abslines) == 0:
        print('Warning: Comp {} does not have any line with EW>{} A! You should consider a '
              'lower -min_ew limit to mask out lines.'.format(comp.name, min_ew.to('AA').value))
    # QtCore.pyqtRemoveInputHook()
    # pdb.set_trace()
    # QtCore.pyqtRestoreInputHook()


def set_fontsize(ax,fsz):
    """
    Generate a Table of columns and so on
    Restrict to those systems where flg_clm > 0

    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)


# Some info about the blending between components
def blending_info(components, specfile, min_vlim=100*u.km/u.s):
    """Computes blending info, and store the components in files that group them together
    depending on overlapping in wobs space. Very useful for speeding up automatic Voigt
    profile fitting

    Parameters
    ----------
    components : list
        A list of AbsComponents
    specfile : str
        Name of the spectrum file, it will be used for writing
    min_vlim : Quantity
        Minimum velocity limit for components for overlap considerations;
        should be equal or larger than the LSF core

    Returns
    -------
        Writes output files.

    """

    # create a copy of component list that has a minimum vlim incorporated
    comps_copy = copy.copy(components)
    for comp in comps_copy:
        for line in comp._abslines:
            # add minimum velocity limit consideration for blends; should be larger than the LSF core
            vlim = line.limits.vlim
            if np.fabs(vlim[0]) < min_vlim:
                vlim[0] = -1 * min_vlim
            if np.fabs(vlim[1]) < min_vlim:
                vlim[1] = min_vlim
            line.limits.set(vlim)

    # QtCore.pyqtRemoveInputHook()
    # pdb.set_trace()
    # QtCore.pyqtRestoreInputHook()

    grouped_comps = ltiu.group_coincident_compoments(comps_copy, output_type='list')
    isolated = []
    grouped = []
    for group in grouped_comps:
        if len(group) == 1:
            isolated += [group[0]]
        else:
            grouped += [group]
    # write to .joenvp output
    ltiu.joebvp_from_components(isolated, specfile, 'isolated.joebvp')
    for ii, group in enumerate(grouped):
        ltiu.joebvp_from_components(group, specfile, 'group_{}.joebvp'.format(ii+1))

    grouped = [isolated] + grouped
    for ii, group in enumerate(grouped):
        if ii == 0:
            print('The following components are isolated in wvobs space:')
            s = ' Isolated: '
        else:
            s = 'Group{}: '.format(ii)
        for comp in group:
            s += '{}, '.format(comp.name)
        s = s[:-2] + '.\n'
        if ii == 1:
            print('The following components overlap in wvobs space:')
        print(s)


# Write JOEVPFIT output from IGM_model.json
def from_igmguesses_to_joebvp(infile, outfile):
    """Reads .json file generated by IGMGuesses and translates it to
    input file(s) for joebvpfit (.joebvp).

    Parameters
    ----------
    infile : str
        Name of the .json file from IGMGuesses
    outfile : str
        Name of the .joebvp file

    Returns:
        A translated version of .json IGMGuesses output
        file that `joebvp` understands.
    """
    import json
    # Read the JSON file
    with open(infile) as data_file:
        igmg_dict = json.load(data_file)

    # Components
    comp_list = []
    for ii, key in enumerate(igmg_dict['cmps'].keys()):
        # QtCore.pyqtRemoveInputHook()
        # pdb.set_trace()
        # QtCore.pyqtRestoreInputHook()
        comp = AbsComponent.from_dict(igmg_dict['cmps'][key], chk_sep=False, chk_data=False, chk_vel=False)
        comp_list += [comp]

    ltiu.joebvp_from_components(comp_list, igmg_dict['spec_file'], outfile)



