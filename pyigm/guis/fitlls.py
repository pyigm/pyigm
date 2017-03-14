"""
#;+
#; NAME:
#; spec_guis
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Spectroscopy Guis with QT
#;      These call pieces from spec_widgets
#;   12-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
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

from astropy.units import Quantity
from astropy import units as u

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import convolve as lsc
import linetools.spectra.io as lsi
from linetools.spectralline import AbsLine
from linetools.guis import utils as ltgu
from linetools import utils as ltu

from pyigm.abssys.lls import LLSSystem
from pyigm.abssys import lls as igmlls
from pyigm.continuum import core as pycc
from pyigm.continuum import quasar as pycq

from linetools.guis import simple_widgets as ltgsm
from linetools.guis import spec_widgets as ltgsp

try:
    ustr = unicode
except NameError:  # For Python 3
    ustr = str

'''
=======
Analyzing spectra with auto_plls

Here is now my preferred approach to searching for
LLS with auto_plls:

1.  Load up the spectrum.  Fiddle with the continuum
normalization (and tilt, if necessary).

2.  Eyeball search for a putative break that one might
associate with a PLLS.

3.  Put the cursor a bit to the left (blueward) of the break
and at the approximate flux level of the data, post-break.
The former sets the starting redshift to search and the latter
sets a guess for NHI.

4. Hit "F" and wait for magic to happen.

5. If an LLS satisfying a rather simple criterion is found, it
should appear.  Otherwise, a message prints to the terminal
stating none found.  You should then inspect the model,
including at Lyb and Lya and modify it (or even delete it).

6. Once you are happy with what you got (hopefully you are),
move on to the next putative break and try again, i.e.
repeat steps 2-5.
'''

# GUI for fitting LLS in a spectrum
class XFitLLSGUI(QMainWindow):
    """ GUI to fit LLS in a given spectrum
        v1.2
        30-Jul-2015 by JXP
    """
    def __init__(self, ispec, zqso, parent=None, lls_fit_file=None,
        outfil=None, smooth=3., fN_gamma=None, template=None,
        dw=0.1, skip_wveval=False, norm=True):
        QMainWindow.__init__(self, parent)
        '''
        ispec : Spectrum1D or specfil
        lls_fit_file: str, optional
          Name of the LLS fit file to input
        smooth : float, optional
          Number of pixels to smooth on (FWHM)
        zqso : float, optional
          Redshift of the quasar.  If input, a Telfer continuum is used
        fN_gamma : float, optional
          Redshift evolution of f(N) or IGM fiddled continuum
        template : str, optional
          Filename of a QSO template to use instead of the Telfer
          continuum. Only used if zqso is also given.
        dw : float, optional
          Pixel width in Angstroms for the wavelength array used to
          generate optical depths. Default is 0.1.
        skip_wveval : bool, optional
          Skip rebinning of wavelengths in the Voigt profile generation.
          This can speed up the code considerably, but use it wisely.
        norm : bool, optional
          Whether to normalize the spectrum by dividing by the
          continuum (default True).
        '''

        # Build a widget combining several others
        self.main_widget = QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        self.update = True
        if outfil is None:
            self.outfil = 'LLS_fit.json'
        else:
            self.outfil = outfil
        self.count_lls = 0
        self.lls_model = None
        self.smooth = None
        self.base_continuum = None
        self.all_forest = []
        self.flag_write = False
        self.dw = float(dw)
        self.skip_wveval = skip_wveval
        if skip_wveval:
            warnings.warn("Skipping wavelength rebinning in Voigt.")
            warnings.warn("Make sure you know what you are doing!")

        # Spectrum
        if isinstance(ispec, XSpectrum1D):
            spec = ispec
            spec_fil = spec.filename
        else:
            # this is broken
            spec, spec_fil = ltgu.read_spec(ispec)


        # LineList
        self.llist = ltgu.set_llist('Strong')
        self.llist['z'] = 0.
        self.plt_wv = zip(np.array([911.7, 949.743, 972.5367,1025.7222,1215.6700])*u.AA,
            ['LL','Lyd', 'Lyg','Lyb','Lya'])

        # z and N boxes
        self.zwidget = ltgsm.EditBox(-1., 'z_LLS=', '{:0.5f}')
        self.Nwidget = ltgsm.EditBox(-1., 'NHI=', '{:0.2f}')
        self.bwidget = ltgsm.EditBox(-1., 'b=', '{:0.1f}')
        self.Cwidget = ltgsm.EditBox('None', 'Comment=', '{:s}')

        # Grab the pieces and tie together
        self.abssys_widg = ltgsp.AbsSysWidget([],only_one=True,
            no_buttons=True, linelist=self.llist[self.llist['List']])

        vlines = [912 * (1 + zqso)]
        self.spec_widg = ltgsp.ExamineSpecWidget(spec,status=self.statusBar,
                                           llist=self.llist, key_events=False,
                                           abs_sys=self.abssys_widg.abs_sys,
                                           vlines=vlines, plotzero=1,
                                                 norm=norm)
        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)
        # Initialize continuum (and LLS if from a file)
        if lls_fit_file is not None:
            self.init_LLS(lls_fit_file,spec)
        else:
            self.conti_dict = pycc.init_conti_dict(
                Norm=float(np.median(spec.flux.value)),
                piv_wv=1215.*(1+zqso),
                #piv_wv2=915.*(1+zqso),
                igm='True')
        if self.base_continuum is None:
            if zqso is not None:
                self.zqso = zqso
                # Read Telfer and apply IGM
                if template is not None:
                    tspec = lsi.readspec(template)
                    # assume wavelengths
                    tspec = XSpectrum1D.from_tuple(
                        (tspec.wavelength.value * (1 + zqso),
                        tspec.flux.value))
                else:
                    tspec = pycq.get_telfer_spec(zqso=zqso,
                              igm=(self.conti_dict['igm']=='True'))
                # Rebin
                self.continuum = tspec.rebin(spec.wavelength)
                # Reset pivot wave
                self.conti_dict['piv_wv'] = 915.*(1+zqso)
                #self.conti_dict['piv_wv'] = 1215.*(1+zqso)
                #self.conti_dict['piv_wv2'] = 915.*(1+zqso)
            else:
                self.zqso = None
                self.continuum = XSpectrum1D.from_tuple((
                    spec.wavelength,np.ones(len(spec.wavelength))))
            self.base_continuum = self.continuum.flux
        self.update_conti()

        self.spec_widg.continuum = self.continuum

        # Full Model (LLS+continuum)
        self.full_model = XSpectrum1D.from_tuple((
            spec.wavelength,np.ones(len(spec.wavelength))))
        if self.smooth is None:
            self.smooth = smooth

        # Initialize as needed
        if lls_fit_file is not None:
            self.update_boxes()
            self.update_model()

        # Outfil
        wbtn = QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.write_out)
        #self.out_box = QtGui.QLineEdit()
        #self.out_box.setText(self.outfil)
        #self.connect(self.out_box, QtCore.SIGNAL('editingFinished ()'), self.set_outfil)

        # Quit
        buttons = QWidget()
        wqbtn = QPushButton('Write\n Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.write_quit)
        qbtn = QPushButton('Quit', self)
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

    def on_list_change(self):
        self.update_boxes()

    def create_status_bar(self):
        self.status_text = QLabel("XFitLLS v0.5.0")
        self.statusBar().addWidget(self.status_text, 1)

    @pyqtSlot()
    def setbzN(self):
        '''Set the column density or redshift from the box
        '''
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
        if self.abssys_widg.all_abssys[idx].NHI > 40:  # Max
            self.abssys_widg.all_abssys[idx].NHI = 40.
            self.Nwidget.set_text(40.)
        # z
        try:
            self.abssys_widg.all_abssys[idx].zabs = (
                float(self.zwidget.box.text()))
        except:
            self.abssys_widg.all_abssys[idx].zabs = 2.
        # b-value
        try:
            self.abssys_widg.all_abssys[idx].bval = (float(self.bwidget.box.text()))*u.km/u.s
        except:
            self.abssys_widg.all_abssys[idx].bval = 10 * u.km/u.s
        if self.abssys_widg.all_abssys[idx].bval < (1*u.km/u.s):
            self.abssys_widg.all_abssys[idx].bval = (1*u.km/u.s)
        self.abssys_widg.all_abssys[idx].comment = (
            self.Cwidget.box.text())
        # Update the lines
        for iline in self.abssys_widg.all_abssys[idx].lls_lines:
            iline.setz(self.abssys_widg.all_abssys[idx].zabs)
            iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
            iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
        # Update the rest
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
        cflux = self.base_continuum * self.conti_dict['Norm']
        # Double tilt
        if 'piv_wv2' in self.conti_dict.keys():
            lowwv = self.continuum.wavelength.value < self.conti_dict['piv_wv2']
            self.continuum.flux[lowwv] = (cflux[lowwv] * (self.continuum.wavelength.value[lowwv]/
                                            self.conti_dict['piv_wv2'])**self.conti_dict['tilt2'])
            self.continuum.flux[~lowwv] = (cflux[~lowwv] * (self.continuum.wavelength.value[~lowwv]/
                                            self.conti_dict['piv_wv'])**self.conti_dict['tilt'])
        else:
            self.continuum.flux = (cflux * (self.continuum.wavelength.value/
                    self.conti_dict['piv_wv'])**self.conti_dict['tilt'])
        if self.lls_model is not None:
            self.full_model.flux = self.lls_model * self.continuum.flux
        # For plotting
        self.spec_widg.spec.co = self.continuum.flux.value

    def update_model(self):
        '''Update absorption model '''
        from linetools.analysis import voigt as lav

        if len(self.abssys_widg.all_abssys) == 0:
            self.lls_model = None
            self.spec_widg.model = None
            return
        # use finer wavelength array to resolve absorption features.
        wa = self.full_model.wavelength
        # Angstroms
        # should really make this a constant velocity width array instead.
        if not self.skip_wveval:
            wa1 = np.arange(wa[0].value, wa[-1].value, self.dw) * wa.unit
        else:
            wa1 = wa
        all_tau_model = igmlls.tau_multi_lls(wa1,
           self.abssys_widg.all_abssys, skip_wveval=self.skip_wveval)

        # Loop on forest lines
        for forest in self.all_forest:
            tau_Lyman = lav.voigt_from_abslines(wa1, forest.lines,
                ret='tau', skip_wveval=self.skip_wveval)
            all_tau_model += tau_Lyman

        # Flux and smooth
        flux = np.exp(-1. * all_tau_model)
        if self.smooth > 0:
            if not self.skip_wveval:
                mult = np.median(np.diff(wa.value)) / self.dw
                flux = lsc.convolve_psf(flux, self.smooth * mult)
            else:
                flux = lsc.convolve_psf(flux, self.smooth)
        if not self.skip_wveval:
            self.lls_model = np.interp(wa.value, wa1.value, flux)
        else:
            self.lls_model = flux

        # Finish
        self.full_model.flux = self.lls_model * self.continuum.flux
        # Over-absorbed
        self.spec_widg.bad_model = np.where( (self.lls_model < 0.7) &
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

    def on_click(self, event):
        """ Over-loads click events
        """
        if event.button == 3:  # Set redshift; equivalent to 'g'
            idx = self.get_sngl_sel_sys()
            if idx is None:
                return
            wrest = event.xdata/(1+self.abssys_widg.all_abssys[idx].zabs)
            awrest = np.array([iline.wrest.value for iline in self.abssys_widg.all_abssys[idx].lls_lines])
            imn = np.argmin(np.abs(wrest-awrest))
            newz = event.xdata/awrest[imn]-1.
            self.abssys_widg.all_abssys[idx].zabs = newz
            # Update the lines
            self.llist['z'] = self.abssys_widg.all_abssys[idx].zabs
            for iline in self.abssys_widg.all_abssys[idx].lls_lines:
                iline.setz(self.abssys_widg.all_abssys[idx].zabs)
                iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
                iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
            # Update the model
            self.update_model()
            # Draw by default
            self.update_boxes()
            self.draw()


    def on_key(self,event):
        if event.key in ['C','1','2','!','@']:  # Set continuum level
            if event.key == 'C':
                imin = np.argmin(np.abs(
                    self.continuum.wavelength.value-event.xdata))
                self.conti_dict['Norm'] = float(event.ydata /
                    (self.base_continuum[imin].value*(event.xdata/
                        self.conti_dict['piv_wv'])**self.conti_dict['tilt']))
            elif event.key == '1':
                self.conti_dict['tilt'] += 0.1
            elif event.key == '2':
                self.conti_dict['tilt'] -= 0.1
            elif event.key == '!':
                self.conti_dict['tilt2'] += 0.1
            elif event.key == '@':
                self.conti_dict['tilt2'] -= 0.1
            self.update_conti()
        elif event.key == 'A':  # New LLS
            # Generate
            z = event.xdata/911.7633 - 1.
            self.add_LLS(z, bval=20.*u.km/u.s, NHI=17.3)
        elif event.key == 'F': # New LLS automagically
            self.auto_plls(event.xdata, event.ydata)
        elif event.key in ['L','a','N','n','v','V','D','$','g']: # LLS-centric
            idx = self.get_sngl_sel_sys()
            if idx is None:
                return
            if event.key == 'L': #LLS
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/911.7633 - 1.
            elif event.key == 'a': #Lya
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/1215.6700-1.
            elif event.key == 'g': # Move nearest line to cursor
                wrest = event.xdata/(1+self.abssys_widg.all_abssys[idx].zabs)
                awrest = np.array([iline.wrest.value for iline in self.abssys_widg.all_abssys[idx].lls_lines])
                imn = np.argmin(np.abs(wrest-awrest))
                newz = event.xdata/awrest[imn]-1.
                self.abssys_widg.all_abssys[idx].zabs = newz
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
                for iline in self.abssys_widg.all_abssys[idx].lls_lines:
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
            self.script_quit = True
        else:
            self.spec_widg.on_key(event)

        # Draw by default
        self.update_boxes()
        self.draw()

    def draw(self):
        self.spec_widg.on_draw(no_draw=True)
        # Add text?
        for kk,lls in enumerate(self.abssys_widg.all_abssys):
            # Label
            ipos = self.abssys_widg.all_items[kk].rfind('_')
            ilbl = self.abssys_widg.all_items[kk][ipos+1:]
            # Add text
            for wv,lbl in self.plt_wv:
                idx = np.argmin(np.abs(self.continuum.wavelength-wv*(1+lls.zabs)))
                self.spec_widg.ax.text(wv.value*(1+lls.zabs),
                    self.continuum.flux[idx],
                    '{:s}_{:s}'.format(ilbl,lbl), ha='center',
                    color='blue', size='small', rotation=90.)
        # Ticks for selected LLS
        idxl = self.get_sngl_sel_sys()
        if idxl is not None:
            lls = self.abssys_widg.all_abssys[idxl]
            # Label
            ipos = self.abssys_widg.all_items[idxl].rfind('_')
            ilbl = self.abssys_widg.all_items[idxl][ipos+1:]
            for line in lls.lls_lines:
                if line.wrest < 915.*u.AA:
                    continue
                idx = np.argmin(np.abs(self.continuum.wavelength-
                    line.wrest*(1+lls.zabs)))
                self.spec_widg.ax.text(line.wrest.value*(1+lls.zabs),
                    self.continuum.flux[idx],
                    '-{:s}'.format(ilbl), ha='center',
                    color='red', size='small', rotation=90.)
        # Draw
        self.spec_widg.canvas.draw()

    def add_forest(self,inp,z):
        '''Add a Lya/Lyb forest line
        '''
        from linetools.isgm.abssystem import GenericAbsSystem
        forest = GenericAbsSystem((0.*u.deg,0.*u.deg), z, [-300.,300.]*u.km/u.s)
        # NHI
        NHI_dict = {'6':12.,'7':13.,'8':14.,'9':15.}
        forest.NHI=NHI_dict[inp]
        # Lines
        for name in ['HI 1215','HI 1025', 'HI 972']:
            aline = AbsLine(name,
                linelist=self.llist[self.llist['List']])
            # Attributes
            aline.attrib['N'] = 10**forest.NHI * u.cm**-2
            aline.attrib['b'] = 20.*u.km/u.s
            aline.setz(forest.zabs)
            # Append
            forest.lines.append(aline)
        # Append to forest lines
        self.all_forest.append(forest)

    def add_LLS(self,z, NHI=17.3,bval=20.*u.km/u.s, comment='None', model=True):
        """Generate a new LLS
        """
        #
        new_sys = LLSSystem((0*u.deg,0*u.deg),z,[-300.,300]*u.km/u.s,NHI=NHI)
        new_sys.bval = bval # This is not standard, but for convenience
        new_sys.comment = comment
        new_sys.fill_lls_lines(bval=bval, do_analysis=0)
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        # Name
        self.count_lls += 1
        new_sys.label = 'LLS_Sys_{:d}'.format(self.count_lls)
        # Add
        self.abssys_widg.add_fil(new_sys.label)
        self.abssys_widg.all_abssys.append(new_sys)
        self.abssys_widg.abslist_widget.item(
            len(self.abssys_widg.all_abssys)).setSelected(True)

        # Update
        self.llist['Plot'] = False # Turn off metal-lines
        if model:  # For dealing with initialization
            self.update_model()

    def auto_plls(self,x,y):
        """Automatically fit a pLLS

        Parameters:
        ----------
        x,y: floats
          x,y values in the GUI
        """
        spec = self.spec_widg.spec # For convenience
        if len(self.abssys_widg.all_abssys) > 0:
            conti= self.full_model
        else:
            conti= self.continuum
        # Generate toy LLS from click
        ximn = np.argmin(np.abs(spec.wavelength.value-x))
        NHI = 17.29 + np.log10(-1.*np.log(y/conti.flux.value[ximn]))
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        #print('NHI={:g}'.format(NHI))
        z = x/(911.7)-1
        plls = LLSSystem((0*u.deg,0*u.deg),z,[-300.,300]*u.km/u.s,NHI=NHI)
        plls.bval = 20*u.km/u.s
        plls.fill_lls_lines(bval=20*u.km/u.s, do_analysis=0)

        # wrest, Tau model, flux
        wrest = spec.wavelength/(1+plls.zabs)
        tau = igmlls.tau_multi_lls(spec.wavelength,[plls])
        emtau = np.exp(-1. * tau)
        lls_flux = lsc.convolve_psf(emtau, 3.)
#xdb.xplot(wrest, lls_flux)

        # zmin (next highest LLS or zem)
        if len(self.abssys_widg.all_abssys) != 0:
            zlls = [lls.zabs for lls in self.abssys_widg.all_abssys if lls.zabs > plls.zabs]
            if len(zlls) == 0:
                zmin = self.zqso+0.01
            else:
                zmin = np.min(np.array(zlls)) - 0.01
        else:
            zmin = self.zqso+0.01

        # Pixels for analysis and rolling
        # NEED TO CUT ON X-Shooter ARM
        apix = np.where( (wrest > 914*u.AA) & #(spec.wavelength<5600*u.AA) &
                        (spec.wavelength<(1+zmin)*1026.*u.AA))[0] # Might go to Lyb
        nroll = (np.argmin(np.abs(spec.wavelength-(911.7*u.AA*(1+zmin))))- # Extra 0.01 for bad z
                   np.argmin(np.abs(spec.wavelength-(911.7*u.AA*(1+plls.zabs)))))
        # Require nroll does not exceed length of spectrum
        if np.max(apix)+nroll > len(spec.wavelength):
            nroll = len(spec.wavelength) - np.max(apix) - 1
        gdpix = np.arange(np.min(apix)-nroll,np.max(apix)+nroll+1)
        roll_flux = np.concatenate([np.ones(nroll),lls_flux[apix], np.ones(nroll)])
        roll_msk = roll_flux < 0.7

        # Generate data arrays
        wave_pad = spec.wavelength[gdpix]
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
        flux_pad = spec.flux[gdpix]
        sig_pad = spec.sig[gdpix]
        if len(self.abssys_widg.all_abssys) > 0:
            conti_pad = conti.flux[gdpix]
        else:
            conti_pad = conti.flux[gdpix]

        # Generate matricies
        flux_matrix = np.zeros((len(roll_flux),nroll))
        sig_matrix = np.zeros((len(roll_flux),nroll))
        conti_matrix = np.zeros((len(roll_flux),nroll))

        roll_matrix = np.zeros((len(roll_flux),nroll))
        mask_matrix = np.zeros((len(roll_flux),nroll))
        for kk in range(nroll):
            roll_matrix[:,kk] = np.roll(roll_flux,kk)
            mask_matrix[:,kk] = np.roll(roll_msk,kk)
            flux_matrix[:,kk] = flux_pad
            conti_matrix[:,kk] = conti_pad
            sig_matrix[:,kk] = sig_pad

        # Model -- Multiply by continuum
        model = roll_matrix * conti_matrix

        # Condition
        idx = np.where( (model < (flux_matrix-sig_matrix*1.5)) & (mask_matrix==True))
        bad_matrix = np.zeros((len(roll_flux),nroll))
        bad_matrix[idx] = 1

        # Sum on offsets and get redshift
        bad = np.sum(bad_matrix,0)
        ibest = np.argmin(bad)
        zbest = spec.wavelength[ibest+ximn]/(911.7*u.AA)-1 # Quantity

        # Add pLLS?
        if bad[ibest] < 10:
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_LLS(zbest.value, bval=20.*u.km/u.s, NHI=NHI)
        else:
            print('No viable pLLS found with our criteria!')


    #def refine_abssys(self):
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

    # Read from a JSON file
    def init_LLS(self,fit_file,spec):
        import json
        # Read the JSON file
        with open(fit_file) as data_file:
            lls_dict = json.load(data_file)
        # Init continuum
        try:
            self.conti_dict = lls_dict['conti_model']
        except KeyError: # Historic
            self.conti_dict = lls_dict['conti']
        else:
            try:
                self.base_continuum = Quantity(lls_dict['conti'])
            except:
                print('Will generate a new base continuum')
                self.base_continuum = None
            else:
                self.continuum = XSpectrum1D.from_tuple((
                    spec.wavelength,np.ones(len(spec.wavelength))))
        #self.update_conti()
        # Check spectra names
        if spec.filename != lls_dict['spec_file']:
            warnings.warn('Spec file names do not match!')
        # LLS
        for key in lls_dict['LLS'].keys():
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_LLS(lls_dict['LLS'][key]['z'],
                NHI=lls_dict['LLS'][key]['NHI'],
                bval=lls_dict['LLS'][key]['bval']*u.km/u.s,
                comment=lls_dict['LLS'][key]['comment'], model=False)
        self.smooth = lls_dict['smooth']
        try:
            self.zqso = lls_dict['zqso']
        except KeyError:
            self.zqso = None
        # Updates
        #self.update_boxes()
        #self.update_model()
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

    # Write
    def write_out(self):
        import json, io
        # Create dict
        out_dict = dict(LLS={},conti_model=self.conti_dict, conti=list(self.base_continuum.value),
            spec_file=self.spec_widg.spec.filename,smooth=self.smooth)
        if self.zqso is not None:
            out_dict['zqso'] = self.zqso
        # Load
        for kk,lls in enumerate(self.abssys_widg.all_abssys):
            key = '{:d}'.format(kk+1)
            out_dict['LLS'][key] = {}
            out_dict['LLS'][key]['z'] = lls.zabs
            out_dict['LLS'][key]['NHI'] = lls.NHI
            out_dict['LLS'][key]['bval'] = lls.lls_lines[0].attrib['b'].value
            out_dict['LLS'][key]['comment'] = str(lls.comment).strip()
        # Write
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        clean_dict = ltu.jsonify(out_dict)
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(ustr(json.dumps(clean_dict, sort_keys=True, indent=4,
                separators=(',', ': '))))
        self.flag_write = True

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.quit()

    # Quit
    def quit(self):
        self.close()



