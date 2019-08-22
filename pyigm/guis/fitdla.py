""" xfitdlagui -- based heavily off of xfitllsgui
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import warnings
import pdb

from PyQt5 import QtGui, QtCore
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget, QLabel, QPushButton, QMainWindow
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout

# Matplotlib Figure object

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.analysis import continuum as laco
from linetools.analysis import interp as laint
from linetools.isgm import utils as ltiu
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import convolve as lsc
import linetools.spectra.io as lsi
from linetools.spectralline import AbsLine
from linetools.guis import utils as ltgu
from linetools.guis import simple_widgets as ltgsm
from linetools.guis import spec_widgets as ltgsp
from linetools import utils as ltu

#from pyigm.abssys.lls import LLSSystem
from pyigm.abssys.dla import DLASystem


'''
=======
Analyzing DLA

Here is now my preferred approach to searching for DLA:

1.  Load up the spectrum.  Fiddle with the continuum
#
'''

# GUI for fitting DLA in a spectrum
class XFitDLAGUI(QMainWindow):
    """ GUI to fit DLA in a given spectrum
        v1.0
        04-Mar-2016 by JXP
    """
    def __init__(self, ispec, parent=None, dla_fit_file=None,
                 zqso=None, outfil=None, smooth=None, dw=0.1, zdla=None,
                 NHI=None,
                 skip_wveval=False, norm=True, conti_file=None):
        QMainWindow.__init__(self, parent)
        """
        ispec : Spectrum1D or specfil
        dla_fit_file: str, optional
          Name of the LLS fit file to input
        smooth : float, optional
          Number of pixels to smooth on (FWHM).  Will not smooth if 0.
        dw : float, optional
          Pixel width in Angstroms for the wavelength array used to
          generate optical depths. Default is 0.1.
        skip_wveval : bool, optional
          Skip rebinning of wavelengths in the Voigt profile generation.
          This can speed up the code considerably, but use it wisely.
        norm : bool, optional
          Whether to normalize the spectrum by dividing by the
          continuum (default True).
        conti_file : str, optional
          ASCII file containing the continuum knots (wave, flux)
        zdla : float, optional
          Create a DLA at the input redshift
        """

        self.help_message = """
Begin by left-clicking in the plot window!

Then any of the following keystroke commands will work:

i,o       : zoom in/out x limits
I,O       : zoom in/out x limits (larger re-scale)
Y         : zoom out y limits
y         : guess y limits
W         : Show original zooom
t,b       : set y top/bottom limit
l,r       : set left/right x limit
a,m,d     : Add/modify/delete continuum knot
A         : Add a new DLA
g         : Move nearest Lyman line to cursor and reset z
N/n       : Increase/decrease NHI
V/v       : Increase/decrease bvalue
Z/z       : Increase/decrease zabs
U/u       : Increase/decrease sig_NHI  [Default is 0.1 dex]
D         : Delete DLA
$         : Toggle displaying metal lines
6,7,8,9   : Add forest lines
?         : Print these help notes
P         : Save current screen in dla_plot.pdf
Q         : Quit the GUI
        """

        # Build a widget combining several others
        self.main_widget = QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        self.update = True
        self.sig_NHI = 0.1
        if outfil is None:
            self.outfil = 'DLA_fit.json'
        else:
            self.outfil = outfil
        self.count_dla = 0
        self.dla_model = None
        if smooth is None:
            self.smooth = 3.  # Pixels to smooth model by
        else:
            self.smooth = smooth
        self.base_continuum = None
        self.all_forest = []
        self.flag_write = False
        self.dw = float(dw)
        self.skip_wveval = skip_wveval
        self.zqso = zqso
        if skip_wveval:
            warnings.warn("Skipping wavelength rebinning in Voigt.")
            warnings.warn("Make sure you know what you are doing!")

        # Spectrum
        if isinstance(ispec, XSpectrum1D):
            spec = ispec
            spec_fil = spec.filename
        else:
            kwargs = {}
            kwargs['rsp_kwargs'] = {}
            kwargs['rsp_kwargs']['masking'] = 'edges'
            spec, spec_fil = ltgu.read_spec(ispec, norm=norm, **kwargs)
        self.spec = spec


        # LineList
        self.llist = ltgu.set_llist('Strong')
        self.llist['z'] = 0.
        self.plt_wv = zip(np.array([972.5367,1025.7222,1215.6700])*u.AA,
            ['Lyg','Lyb','Lya'])

        # z and N boxes
        self.zwidget = ltgsm.EditBox(-1., 'z_DLA=', '{:0.5f}')
        self.Nwidget = ltgsm.EditBox(-1., 'NHI=', '{:0.2f}')
        self.bwidget = ltgsm.EditBox(-1., 'b=', '{:0.1f}')
        self.Cwidget = ltgsm.EditBox('None', 'Comment=', '{:s}')

        # Grab the pieces and tie together
        self.abssys_widg = ltgsp.AbsSysWidget([],only_one=True,
            no_buttons=True, linelist=self.llist[self.llist['List']])

        self.spec_widg = ltgsp.ExamineSpecWidget(spec,status=self.statusBar,
                                           llist=self.llist, key_events=False,
                                           abs_sys=self.abssys_widg.abs_sys,
                                           plotzero=1, norm=norm)
        # Create a DLA?
        if zdla is not None:
            self.add_DLA(zdla, NHI=NHI, model=False)

        # Initialize continuum (and DLA if from a file)
        if dla_fit_file is not None:
            self.init_DLA(dla_fit_file,spec)
        else:
            if conti_file == 'from_spectrum':
                co = spec.co
                # knots: list of (x,y) pairs, giving the position of
                # spline knots used to generate the continuum
                #knots =
                nknots = 100
                knots = []
                for ii in range(nknots):
                    ii2 = round(ii*len(spec.wavelength)/nknots)
                    knots.append((spec.wavelength[ii2].value,spec.co[ii2]))
                # add all min and max points
                for ii in range(len(spec.wavelength)-2):
                    if ((spec.co[ii] - spec.co[ii+1])*(spec.co[ii+2] - spec.co[ii+1]) > 0):
                        iknot = (spec.wavelength[ii].value,spec.co[ii])
                        if iknot not in knots:
                            knots.append(iknot)

            else:
                if conti_file is not None:
                    # Read continuum
                    cspec = lsi.readspec(conti_file)
                    if not cspec.sig_is_set:
                        cspec.sig = 0.1*np.median(cspec.flux)
                else:
                    cspec = spec
                if zqso is not None:
                    co, knots = laco.find_continuum(cspec, redshift=self.zqso)
                else:
                    co, knots = laco.find_continuum(cspec, kind='default')

            self.conti_dict = dict(co=co, knots=knots)

        self.update_conti()

        #self.spec_widg.continuum = self.continuum

        # Full Model (LLS+continuum)
        self.full_model = XSpectrum1D.from_tuple((
            spec.wavelength,np.ones(len(spec.wavelength))))

        # Initialize as needed
        if (dla_fit_file is not None) or (zdla is not None):
            self.update_boxes()
            self.update_model()

        # Outfil
        wbtn = QPushButton(self)
        wbtn.setText('Write')
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.write_out)
        #self.out_box = QtGui.QLineEdit()
        #self.out_box.setText(self.outfil)
        #self.connect(self.out_box, QtCore.SIGNAL('editingFinished ()'), self.set_outfil)

        # Quit
        buttons = QWidget()
        wqbtn = QPushButton(self)
        wqbtn.setText('Write\n Quit')
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.write_quit)
        qbtn = QPushButton(self)
        qbtn.setText('Quit')
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.quit)

        # Connections (buttons are above)
        self.spec_widg.canvas.mpl_connect('key_press_event', self.on_key)
        self.abssys_widg.abslist_widget.itemSelectionChanged.connect(
            self.on_list_change)
        self.Nwidget.box.textChanged[str].connect(self.setbzN)
        self.zwidget.box.textChanged[str].connect(self.setbzN)
        self.bwidget.box.textChanged[str].connect(self.setbzN)
        self.Cwidget.box.textChanged[str].connect(self.setbzN)

        # Layout
        anly_widg = QWidget()
        anly_widg.setMaximumWidth(400)
        anly_widg.setMinimumWidth(250)

        # Write/Quit buttons
        hbox1 = QHBoxLayout()
        hbox1.addWidget(wbtn)
        hbox1.addWidget(wqbtn)
        hbox1.addWidget(qbtn)
        buttons.setLayout(hbox1)

        # z,N
        zNwidg = QWidget()
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.zwidget)
        hbox2.addWidget(self.Nwidget)
        hbox2.addWidget(self.bwidget)
        zNwidg.setLayout(hbox2)
        #vbox.addWidget(self.pltline_widg)

        vbox = QVBoxLayout()
        vbox.addWidget(zNwidg)
        vbox.addWidget(self.Cwidget)
        vbox.addWidget(self.abssys_widg)
        vbox.addWidget(buttons)
        anly_widg.setLayout(vbox)

        hbox = QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(anly_widg)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

        #self.spec_widg.setFixedWidth(900)
        self.spec_widg.setMinimumWidth(900)

        # Print help message
        print(self.help_message)

    @pyqtSlot()
    def on_list_change(self):
        self.update_boxes()

    def create_status_bar(self):
        self.status_text = QLabel("XFitLLS v0.5.0")
        self.statusBar().addWidget(self.status_text, 1)

    @pyqtSlot()
    def setbzN(self):
        """Set the column density or redshift from the box
        """
        if self.update is False:
            return
        idx = self.get_sngl_sel_sys()
        if idx is None:
            return
        # NHI
        try:
            self.abssys_widg.all_abssys[idx].NHI = (
                float(self.Nwidget.box.text()))
        except:
            self.abssys_widg.all_abssys[idx].NHI = 20.
        # z
        try:
            self.abssys_widg.all_abssys[idx].limits._z = (
                float(self.zwidget.box.text()))
        except:
            self.abssys_widg.all_abssys[idx].limits._z = 2.
        # b-value
        try:
            self.abssys_widg.all_abssys[idx].bval = (
                float(self.bwidget.box.text()))*u.km/u.s
        except:
            self.abssys_widg.all_abssys[idx].bval = 10 * u.km/u.s
        self.abssys_widg.all_abssys[idx].comment = (
            self.Cwidget.box.text())
        # Update the lines
        for iline in self.abssys_widg.all_abssys[idx].dla_lines:
            iline.setz(self.abssys_widg.all_abssys[idx].zabs)
            iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
            iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
        # Update the rest
        self.llist['z'] = self.abssys_widg.all_abssys[idx].zabs
        self.update_model()
        self.draw()

    def update_boxes(self):
        """Update Nbz boxes"""
        idx = self.get_sngl_sel_sys()
        if idx is None:
            return
        # z
        self.update = False  # Avoids a bit of an internal loop
        self.zwidget.box.setText(
            self.zwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].zabs))
        # N
        self.Nwidget.box.setText(
            self.Nwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].NHI))
        # b
        self.bwidget.box.setText(
            self.bwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].bval.value))
        # Comment
        self.Cwidget.box.setText(
            self.Cwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].comment))
        # Rest-frame too
        self.spec_widg.show_restframe = False
        self.spec_widg.rest_z = self.abssys_widg.all_abssys[idx].zabs
        self.update = True

    def update_conti(self):
        """Update continuum
        """
        # Update continuum from knots
        x,y = zip(*self.conti_dict['knots'])
        spl = laint.AkimaSpline(x, y)
        self.conti_dict['co'] = spl(self.spec.wavelength.value)

        if self.dla_model is not None:
            self.full_model.flux = self.dla_model * self.conti_dict['co']
        # For plotting
        self.spec_widg.spec.co = self.conti_dict['co']

    def update_model(self):
        """Update absorption model"""
        from linetools.analysis import voigt as lav

        if len(self.abssys_widg.all_abssys) == 0:
            self.dla_model = None
            self.spec_widg.model = None
            return
        # use finer wavelength array to resolve absorption features.
        wa = self.full_model.wavelength
        # Angstroms
        # should really make this a constant velocity width array instead.
        if not self.skip_wveval:
            #wa1 = np.arange(wa[0].value, wa[-1].value, self.dw) * wa.unit
            wa1 = np.arange(self.full_model.wvmin.value, self.full_model.wvmax.value, self.dw) * wa.unit
        else:
            wa1 = wa
        #all_tau_model = igmlls.tau_multi_lls(wa1,
        #   self.abssys_widg.all_abssys, skip_wveval=self.skip_wveval)
        all_lines = []
        for abssys in self.abssys_widg.all_abssys:
            for iline in abssys.dla_lines:
                all_lines.append(iline)
        #QtCore.pyqtRemoveInputHook()
        #import pdb; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        if all_lines[0].attrib['b'].value < 0.:
            QtCore.pyqtRemoveInputHook()
            import pdb; pdb.set_trace()
            QtCore.pyqtRestoreInputHook()
        tau_Lyman = lav.voigt_from_abslines(wa1, all_lines, ret='tau', skip_wveval=self.skip_wveval)
        all_tau_model = tau_Lyman
        # Loop on forest lines
        for forest in self.all_forest:
            tau_Lyman = lav.voigt_from_abslines(wa1, forest, #.lines,
                ret='tau', skip_wveval=self.skip_wveval)
            all_tau_model += tau_Lyman

        # Uncertainty
        tau_low = all_tau_model * 10**(-1*self.sig_NHI)
        tau_high = all_tau_model * 10**(self.sig_NHI)

        # Flux and smooth
        flux = np.exp(-1. * all_tau_model)
        low_flux = np.exp(-1. * tau_low)
        high_flux = np.exp(-1. * tau_high)
        if self.smooth > 0.:
            if not self.skip_wveval:
                mult = np.median(np.diff(wa.value)) / self.dw
                flux = lsc.convolve_psf(flux, self.smooth * mult)
                low_flux = lsc.convolve_psf(low_flux, self.smooth * mult)
                high_flux = lsc.convolve_psf(high_flux, self.smooth * mult)
            else:
                flux = lsc.convolve_psf(flux, self.smooth)
                low_flux = lsc.convolve_psf(low_flux, self.smooth)
                high_flux = lsc.convolve_psf(high_flux, self.smooth)
        if not self.skip_wveval:
            self.dla_model = np.interp(wa.value, wa1.value, flux)
            self.dla_low = np.interp(wa.value, wa1.value, low_flux)
            self.dla_high = np.interp(wa.value, wa1.value, high_flux)
        else:
            self.dla_model = flux
            self.dla_low = low_flux
            self.dla_high = high_flux

        # Finish
        self.full_model.flux = self.dla_model * self.conti_dict['co']
        self.dla_low *= self.conti_dict['co']
        self.dla_high *= self.conti_dict['co']
        # Over-absorbed
        #QtCore.pyqtRemoveInputHook()
        #import pdb; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.spec_widg.bad_model = np.where( (self.dla_model < 0.7) &
            (self.full_model.flux < (self.spec_widg.spec.flux-
                self.spec_widg.spec.sig*1.5)))[0]
        # Model
        self.spec_widg.model = self.full_model

    def get_sngl_sel_sys(self):
        '''Grab selected system
        '''
        items = self.abssys_widg.abslist_widget.selectedItems()
        if len(items) == 0:
            return None
        elif len(items) > 1:
            print('Need to select only 1 system!')
            return None
        #
        item = items[0]
        txt = item.text()
        if txt == 'None':
            return None
        else:
            idx = self.abssys_widg.all_items.index(item.text())
            return idx

    def on_key(self,event):
        if event.key in ['a','m','d']: # Modify knots
            if event.key == 'a':  # Add a knot
                x, y = event.xdata, event.ydata
                self.conti_dict['knots'].append([x, float(y)])
                self.conti_dict['knots'].sort()
            elif event.key in ['m','d']:
                # Find cloases
                contx,conty = zip(*self.spec_widg.ax.transData.transform(
                        self.conti_dict['knots']))
                #sep = np.hypot(event.x - np.array(contx),
                #               event.y - np.array(conty))
                sep = np.abs(event.x - np.array(contx)),
                ind = np.argmin(sep)
                #
                if event.key == 'm':
                    x, y = event.xdata, event.ydata
                    self.conti_dict['knots'][ind] = [x, float(y)]
                else: # Delete
                    self.conti_dict['knots'].pop(ind)
                    #QtCore.pyqtRemoveInputHook()
                    #pdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                self.conti_dict['knots'].sort()
            self.update_conti()
        elif event.key == 'A': # New DLA
            # Generate
            z = event.xdata/1215.670 - 1.
            self.add_DLA(z)
        elif event.key in ['u','U','a','N','n','v','V','D','$','g','z','Z']: # DLA-centric
            idx = self.get_sngl_sel_sys()
            if idx is None:
                return
            elif event.key == 'g': # Move nearest line to cursor
                wrest = event.xdata/(1+self.abssys_widg.all_abssys[idx].zabs)
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                awrest = np.array([iline.wrest.value for iline in self.abssys_widg.all_abssys[idx].dla_lines])
                imn = np.argmin(np.abs(wrest-awrest))
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/awrest[imn]-1.
            elif event.key == 'Z': #Add to redshift
                self.abssys_widg.all_abssys[idx].limits._z += 0.0002
            elif event.key == 'z': #Subtract from redshift
                self.abssys_widg.all_abssys[idx].limits._z -= 0.0002
            elif event.key == 'U': #Add to sig_NHI
                self.sig_NHI += 0.05
                print('sig_NHI={}'.format(self.sig_NHI))
            elif event.key == 'u': #Subtract from sig_NHI
                self.sig_NHI -= 0.05
                print('sig_NHI={}'.format(self.sig_NHI))
            elif event.key == 'N': #Add to NHI
                self.abssys_widg.all_abssys[idx].NHI += 0.05
            elif event.key == 'n': #Subtract from NHI
                self.abssys_widg.all_abssys[idx].NHI -= 0.05
            elif event.key == 'v': #Subtract from bval
                self.abssys_widg.all_abssys[idx].bval -= 2*u.km/u.s
            elif event.key == 'V': #Add to bval
                self.abssys_widg.all_abssys[idx].bval += 2*u.km/u.s
            elif event.key == 'D': # Delete system
                self.abssys_widg.remove_item(idx)
                idx = None
            elif event.key == '$': # Toggle metal-lines
                self.llist['Plot'] = not self.llist['Plot']
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
            else:
                raise ValueError('Not ready for this keystroke')
            # Update the lines
            if idx is not None:
                self.llist['z'] = self.abssys_widg.all_abssys[idx].zabs
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                for iline in self.abssys_widg.all_abssys[idx].dla_lines:
                    iline.setz(self.abssys_widg.all_abssys[idx].zabs)
                    iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
                    iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
            # Update the model
            self.update_model()
        elif event.key in ['6','7','8','9']: # Add forest line
            self.add_forest(event.key,event.xdata/1215.6701 - 1.)
            self.update_model()
        elif event.key == 'Q':
            # set an attribute to tell calling script to abort
            print("Setting the quit attribute, the calling script should "
                  "abort after you close the GUI")
            self.quit()
        elif event.key == '?':
            print(self.help_message)
        elif event.key == 'S':
            print('Smoothing not allowed in this GUI')
        elif event.key == 'P':  # Print
            self.draw(save=True)
            return
        else:
            self.spec_widg.on_key(event)

        # Draw by default
        self.update_boxes()
        self.draw()

    def draw(self, save=False):
        self.spec_widg.on_draw(no_draw=True)
        # Add text?
        for kk,dla in enumerate(self.abssys_widg.all_abssys):
            # Label
            ipos = self.abssys_widg.all_items[kk].rfind('_')
            ilbl = self.abssys_widg.all_items[kk][ipos+1:]
            # Add text
            for wv,lbl in self.plt_wv:
                idx = np.argmin(np.abs(self.spec.wavelength-wv*(1+dla.zabs)))
                self.spec_widg.ax.text(wv.value*(1+dla.zabs),
                    self.conti_dict['co'][idx],
                    '{:s}_{:s}'.format(ilbl,lbl), ha='center',
                    color='blue', size='small', rotation=90.)
        # Ticks for selected DLA
        idxl = self.get_sngl_sel_sys()
        if idxl is not None:
            dla = self.abssys_widg.all_abssys[idxl]
            # Label
            ipos = self.abssys_widg.all_items[idxl].rfind('_')
            ilbl = self.abssys_widg.all_items[idxl][ipos+1:]
            for line in dla.dla_lines:
                if line.wrest < 915.*u.AA:
                    continue
                idx = np.argmin(np.abs(self.spec.wavelength-
                    line.wrest*(1+dla.zabs)))
                self.spec_widg.ax.text(line.wrest.value*(1+dla.zabs),
                    self.conti_dict['co'][idx],
                    '-{:s}'.format(ilbl), ha='center',
                    color='red', size='small', rotation=90.)
            # Error interval
            #QtCore.pyqtRemoveInputHook()
            #import pdb; pdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.spec_widg.ax.fill_between(self.full_model.wavelength.value, self.dla_low, self.dla_high, color='gray', alpha=0.5)
        # Continuum knots
        xknot = np.array([knot[0] for knot in self.conti_dict['knots']])
        yknot = np.array([knot[1] for knot in self.conti_dict['knots']])
        self.spec_widg.ax.scatter(xknot, yknot, marker='o', color='pink',
                                  zorder=10)
        # Draw
        self.spec_widg.canvas.draw()
        if save:
            print("Printing the window to dla_plot.pdf")
            self.spec_widg.canvas.print_figure('dla_plot.pdf')
        #QtCore.pyqtRemoveInputHook()
        #import pdb; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

    def add_forest(self,inp,z):
        """Add a Lya/Lyb forest line
        """
        forest = []
        # NHI
        NHI_dict = {'6':12.,'7':13.,'8':14.,'9':15.}
        forest_NHI=NHI_dict[inp]
        # Lines
        for name in ['HI 1215','HI 1025', 'HI 972']:
            aline = AbsLine(name, linelist=self.llist[self.llist['List']], z=z)
            # Attributes
            aline.attrib['N'] = 10**forest_NHI * u.cm**-2
            aline.attrib['b'] = 20.*u.km/u.s
            # Append
            forest.append(aline)
        # Append to forest lines
        self.all_forest.append(forest)

    def add_DLA(self,z, NHI=None, bval=30.*u.km/u.s, comment='None', model=True):
        """Generate a new DLA
        """
        if NHI is None:
            NHI = 20.3
        # Lya, Lyb
        dla_lines = []  # For convenience
        for trans in ['HI 1025', 'HI 1215']:
            iline = AbsLine(trans, z=z)
            iline.attrib['flag_N'] = 1
            iline.attrib['N'] = 10**NHI / u.cm**2
            iline.attrib['sig_N'] = 1 / u.cm**2  # Avoid nan
            iline.attrib['b'] = bval
            iline.attrib['coord'] = SkyCoord(ra=0*u.deg,dec=0*u.deg)
            iline.analy['do_analysis'] = 0  # Avoids plotting
            dla_lines.append(iline)
        # Generate system
        HIcomponent = ltiu.build_components_from_abslines(dla_lines)[0]
        HIcomponent.synthesize_colm(overwrite=True)
        new_sys = DLASystem.from_components([HIcomponent]) #(0*u.deg,0*u.deg),z,None,NHI)
        new_sys.bval = bval # This is not standard, but for convenience
        new_sys.comment = comment
        new_sys.dla_lines = dla_lines  # Also for convenience
        # Name
        self.count_dla += 1
        new_sys.label = 'DLA_Sys_{:d}'.format(self.count_dla)
        # Add
        self.abssys_widg.add_fil(new_sys.label)
        self.abssys_widg.all_abssys.append(new_sys)
        #QtCore.pyqtRemoveInputHook()
        #import pdb; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.abssys_widg.abslist_widget.item(len(
            self.abssys_widg.all_abssys)).setSelected(True)

        # Update
        self.llist['Plot'] = False # Turn off metal-lines
        if model:  # For dealing with initialization
            self.update_model()

    #def refine_abesys(self):
    #    item = self.abssys_widg.abslist_widget.selectedItems()
    #    if len(item) != 1:
    #        self.statusBar().showMessage('AbsSys: Must select only 1 system!')
    #        print('AbsSys: Must select only 1 system!')
    #    txt = item[0].text()
    #    ii = self.abssys_widg.all_items.index(txt)
    #    iabs_sys = self.abssys_widg.all_abssys[ii]
    #    # Launch
    #    gui = xspg.XVelPltGui(self.spec_widg.spec, outfil=iabs_sys.absid_file,
    #                           abs_sys=iabs_sys, norm=self.spec_widg.norm)
    #    gui.exec_()


    def init_DLA(self,fit_file,spec):
        """ Read from a JSON file
        Parameters
        ----------
        fit_file
        spec

        Returns
        -------

        """
        import json
        # Read the JSON file
        with open(fit_file) as data_file:
            dla_dict = json.load(data_file)
        # Init continuum
        self.conti_dict = dla_dict['conti_model']
        self.update_conti()
        # Check spectra names
        if spec.filename != dla_dict['spec_file']:
            warnings.warn('Spec file names do not match!')
        # LLS
        for key in dla_dict['DLA'].keys():
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_DLA(dla_dict['DLA'][key]['z'],
                NHI=dla_dict['DLA'][key]['NHI'],
                bval=dla_dict['DLA'][key]['bval']*u.km/u.s,
                comment=dla_dict['DLA'][key]['comment'], model=False)
        self.smooth = dla_dict['smooth']
        try:
            self.zqso = dla_dict['zqso']
        except KeyError:
            self.zqso = None
        # Updates
        #self.update_boxes()
        #self.update_model()
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

    # Write
    @pyqtSlot()
    def write_out(self):
        import json, io
        # Create dict
        out_dict = dict(DLA={},conti_model=self.conti_dict,
            spec_file=self.spec_widg.spec.filename,smooth=self.smooth)
        if self.zqso is not None:
            out_dict['zqso'] = self.zqso
        # Load
        for kk,dla in enumerate(self.abssys_widg.all_abssys):
            key = '{:d}'.format(kk+1)
            out_dict['DLA'][key] = {}
            out_dict['DLA'][key]['z'] = dla.zabs
            out_dict['DLA'][key]['NHI'] = dla.NHI
            out_dict['DLA'][key]['bval'] = dla.dla_lines[0].attrib['b'].value
            out_dict['DLA'][key]['comment'] = str(dla.comment).strip()
        # Write
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        clean_dict = ltu.jsonify(out_dict)
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(json.dumps(clean_dict, sort_keys=True, indent=4,
                separators=(',', ': ')))
        print("Wrote fit and continuum to {:s}".format(self.outfil))
        self.flag_write = True

    # Write + Quit
    @pyqtSlot()
    def write_quit(self):
        self.write_out()
        self.quit()

    # Quit
    @pyqtSlot()
    def quit(self):
        self.close()



