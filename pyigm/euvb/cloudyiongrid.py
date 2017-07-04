"""

Code to generate large grids of cloudy models under different assumptions

Written by Michele Fumagalli in Durham, Summer 2015
michele.fumagalli@durham.ac.uk

Copyright (C) 2015 Michele Fumagalli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


""" 

#general import 
import numpy as np 
from astropy import constants as const
from scipy import interpolate
# import mypython as mp ##comment out cbw
from scipy import integrate 
from scipy.interpolate import interp1d
import os
import re

from pkg_resources import resource_filename


class Cldyion():

    """
    Define a class to run ionization modeling to do cloudy grids of elements 
    By default, use the HM2005. Other UVB that can be included are: HM2012, FG09

    """

    def __init__(self,uvb='HM05'):

        #This will be constants 
        self.column=False
        self.metal=False
        self.redshift=False
        self.density=False
        self.jgal=False
        self.jqso=False
        self.temperature=False
        self.fstar=False
        self.dustscale=False
        self.root=''
        self.path=''
        self.uvbtype=uvb
        self.HM05write=''
        self.subdir=''

        #this is stuff in output of a model (dictionaries)
        self.ins={}  #mean ion fraction
        self.clm={}  #column densities 
        self.prof={} #radial profiles
        self.abund={} #abundances


        #A dictionary that stores values of frequency log [Hz] and mean intensity (log 4 pi J_nu [erg/s/Hz/cm^2])
        self.source={'lgnu':0,'lgfpiJ':0}

        #A dictionary that stores correction fraction compared to the
        #Asplund 2009 abundances (1 is default, <1 is depleted)
        #Runs over the elements one wants to deplete
        self.abundance=False


    def minimal(self,column,metal,redshift,density):

        """
        Specify a minimal cloudy model (see Fumagalli et al. 2015)

        column     -> the HI column density [log]
        metal      -> the log metallicity in solar unit [X/H] 
        redshift   -> the redshift
        density    -> the total hydrogen volume density [log]

        """

        #Start by defining defaults
        self.column=column
        self.metal=metal
        self.redshift=redshift
        self.density=density 
        self.root="n{0}_r{1}_c{2}_Z{3}".format(self.density,self.redshift,self.column,self.metal)

        #leave elements as they are in minimal model 
        #Initialize the HM05 UVB at the given redshift 
        self.inituvb()

    def cie(self,column,metal,redshift,density,temperature):


        """
        Specify a cloudy model with constant temperature [CIE] (see Fumagalli et al. 2015)

        column     -> the HI column density [log]
        metal      -> the log metallicity in solar unit [X/H] 
        redshift   -> the redshift
        density    -> the total hydrogen volume density [log]
        temperature -> the log of the temperature in K 

        """

        #first init the minimal model 
        self.minimal(column,metal,redshift,density)

        #update as needed
        self.temperature=temperature
        self.root="n{0}_r{1}_c{2}_Z{3}_T{4}".format(self.density,self.redshift,self.column,self.metal,self.temperature)

    def dust(self,column,metal,redshift,density,fstar):


        """
        Specify a cloudy model with dust depletion (see Fumagalli et al. 2015)

        column     -> the HI column density [log]
        metal      -> the log metallicity in solar unit [X/H] 
        redshift   -> the redshift
        density    -> the total hydrogen volume density [log]
        fstar      -> the strength of the dust deplition [linear]

        """

        #first init the minimal model 
        self.minimal(column,metal,redshift,density)

        #update as needed global variable
        self.fstar=fstar
        self.root="n{0}_r{1}_c{2}_Z{3}_f{4}".format(self.density,self.redshift,self.column,self.metal,self.fstar)

        #now compute the dust scaling
        if(self.fstar > 0):
            #If fstar > 0 ,add dust
            self.dustscale=self.fstar*(10**self.metal)
        else:
            #Else, do not include dust
            self.dustscale=False

        #now compute the abundances depletion following the method by Jenkins et al. 
        #Going to f = -1 brings the corrections close to 0
        ejelm=['Carbon','Nitrogen','Oxygen','Magnesium','Silicon','Phosphorus','Chlorine','Titanium',\
                        'Chromium','Manganese','Iron','Nickel','Copper','Zinc']
        ejAx=[-0.101,-0.000,-0.225,-0.997,-1.136,-0.945,-1.242,-2.048,-1.447,-0.857,-1.285,-1.490,-0.710,-0.610]
        ejBx=[-0.193,-0.109,-0.145,-0.800,-0.570,-0.166,-0.314,-1.957,-1.508,-1.354,-1.513,-1.829,-1.102,-0.279]
        ejzx=[0.803,0.550,0.598,0.531,0.305,0.488,0.609,0.430,0.470,0.520,0.437,0.599,0.711,0.555]

        #init depletion
        self.abundance={}

        #compute values
        for ee in range(len(ejelm)):
            depl=ejBx[ee]+ejAx[ee]*(self.fstar-ejzx[ee])
            if(depl < 0):
                self.abundance[ejelm[ee]]=10**depl
            else:
                self.abundance[ejelm[ee]]=1.0


    def modsource(self,column,metal,redshift,density,jgal,jqso):


        """
        Specify a cloudy model with varying SED properties 

        column     -> the HI column density [log]
        metal      -> the log metallicity in solar unit [X/H] 
        redshift   -> the redshift
        density    -> the total hydrogen volume density [log]
        jgal       -> weight to add in the galaxy contribution
        jqso       -> weight to add in the qso contribution

        """

        #import numpy as np
        #from astropy import constants as const
        #from scipy import interpolate
        #import mypython as mp

        #first init the minimal model 
        self.minimal(column,metal,redshift,density)

        #make copy of uvb
        uvb_nu=self.source['lgnu']
        uvb_f=self.source['lgfpiJ']

        #next update coefficient and name
        self.jgal=jgal
        self.jqso=jqso
        self.root="n{0}_r{1}_c{2}_Z{3}_g{4}_q{5}".format(self.density,self.redshift,self.column,self.metal,self.jgal,self.jqso)

        ##################
        ###  quasar    ###
        ##################

        #now load qso sed
        qso_nu=[]
        qso_f=[]

        fil=open('Scripts/richards_qsosed.txt')

        for ll in fil:
            fields=ll.split()
            try:
                qso_nu.append(float(fields[0])) #log Hz
                qso_f.append(float(fields[1]))  #log erg/s
            except:
                #skip header
                hhd=0

        fil.close()

        #extend to max frequency for UVB with alpha=-1
        qso_nu.append(np.max(self.source['lgnu']))
        qso_f.append(qso_f[-1])

        #go to np arrays
        qso_f=np.array(qso_f)
        qso_nu=np.array(qso_nu)

        #make qso sed comparable to the UVB in units 
        #Cuba is in log of Hz; log of erg/s/cm2/Hz

        #Assume inner radius of 100 kpc for flux conversion
        inrad=1e5 #pc
        qso_f=qso_f-qso_nu-np.log10(4*np.pi)-2*np.log10(inrad*const.pc.to('cm').value)

        #Now apply scaling 
        qso_f=qso_f+self.jqso

        #add to model the quasar contrinution 
        qsoint=interpolate.interp1d(qso_nu,qso_f,fill_value=-99,bounds_error=False)
        qso_f_int=qsoint(self.source['lgnu'])
        self.source['lgfpiJ']=np.log10(10**self.source['lgfpiJ']+10**qso_f_int)

        ##################
        ###  galaxy    ###
        ##################
        # CBW # gal=mp.seds.sb99.Sb99()
        gal.loadmod('sfr1_ge_z014',path='./Scripts/')

        #grab the last time 
        gal_lam=np.array(gal.spectrum[-1]['WAVE']) #A
        gal_f=np.array(gal.spectrum[-1]['TOTAL']) #log erg/s/A

        #convert in unit similar to qso
        gal_nu=np.log10(const.c.to('cm/s').value*1e8/gal_lam) #log Hz
        gal_f=np.log10(10**gal_f*gal_lam/10**gal_nu) #log erg/s/Hz

        #go to flux at given distance 
        gal_f=gal_f-np.log10(4*np.pi)-2*np.log10(inrad*const.pc.to('cm').value) #log erg/s/Hz/cm2

        #apply scaling 
        gal_f=gal_f+self.jgal

        #add to model 
        galint=interpolate.interp1d(gal_nu,gal_f,fill_value=-99,bounds_error=False)
        gal_f_int=galint(self.source['lgnu'])
        self.source['lgfpiJ']=np.log10(10**self.source['lgfpiJ']+10**gal_f_int)


        ##DEBUG sequence 
        #import matplotlib as mpl
        #import matplotlib.pyplot as pl
        #pl.plot(uvb_nu,uvb_f,linestyle=':', linewidth=2)
        #pl.plot(self.source['lgnu'],qso_f_int,linestyle=':', linewidth=2)
        #pl.plot(self.source['lgnu'],gal_f_int,linestyle=':', linewidth=2)
        #pl.plot(self.source['lgnu'],self.source['lgfpiJ'], linewidth=4)
        #pl.show()
        #exit()

        #pl.plot(self.source['lgnu'],10**(self.source['lgfpiJ']-uvb_f))
        #pl.show()

        #import pickle
        #src={'lgnu':self.source['lgnu'],'all':self.source['lgfpiJ'],'uvb':uvb_f,'gal':gal_f_int,'qso':qso_f_int}
        #fil=open('source_debug.pkl','w')
        #pickle.dump(src,fil)
        #fil.close()
        ## END of DEBUG sequence

    def writein(self,basedir):

        """
        Write the init file 

        """
        #from scipy import integrate 
        #from astropy import constants as const 
        #from scipy.interpolate import interp1d
        #import numpy as np

        ##CBW separate into sub-directories
        self.subdir="n{0}/r{1}/".format(self.density,self.redshift)
        self.path=basedir+"/"+self.subdir
        if not os.path.isdir(self.path+"/"):
            os.makedirs(self.path+"/")

        ##CBW don't write to stdout every time
        # print('Write {}'.format(self.root))

        fil=open(self.path+"/"+self.root+".in","w")

        #Write header
        fil.write("##Header\n")
        fil.write("# Metallicity {0}\n".format(self.metal))
        fil.write("# Density {0}\n".format(self.density))
        fil.write("# Redshift {0}\n".format(self.redshift))
        fil.write("# Column {0}\n".format(self.column))
        fil.write("# Jgal {0}\n".format(self.jgal))
        fil.write("# Jqso {0}\n".format(self.jqso))
        fil.write("# Temp {0}\n".format(self.temperature))
        fil.write("# Fstar {0}\n".format(self.fstar))

        #set the metallicity
        fil.write("##Abundances\n")
        #turn on default
        fil.write("abundances GASS10 no grains\n") 
        #scale all elements
        fil.write("metals {0} log \n".format(self.metal)) 
        #loop over elemnts for individual corrections
        if(self.abundance):
            for dep in self.abundance:
                fil.write("element scale factor {0} {1}\n".format(dep,self.abundance[dep]))

        #if needed specify dust
        if(self.dustscale):
            fil.write("grains ism {}\n".format(self.dustscale))
            fil.write("no grain qheat\n")

        #density/geometry
        fil.write("##Geometry/density\n")
        fil.write("hden {0}\n".format(self.density))

        #Include constant temperature for CIE
        if(self.temperature):
            fil.write("constant temperature {0} log\n".format(self.temperature))

        #radiation
        fil.write("##Radiation\n")
        fil.write("cmb, z={0}\n".format(self.redshift))

        #if using HM05, write the call to the spectrum here
        if(self.HM05write != ""):
            fil.write(self.HM05write)

        # #compute integral over range free of difficulties
        # #Define min-max in Hz
        # maxnu=(const.c.to('cm/s')/(1250*1e-8)).value
        # ionnu=(const.c.to('cm/s')/(2000*1e-8)).value

        # #make sure interval is in range of model, and convert to Ry
        # ionnu = np.max((ionnu,10**np.min(self.source['lgnu'])))
        # maxry=(4.1356675e-15*maxnu/13.605698066)
        # minry=(4.1356675e-15*ionnu/13.605698066)

        # #integrate
        # fint = interp1d(10**self.source['lgnu'],10**self.source['lgfpiJ'])
        # integral = np.log10(integrate.quad(fint,ionnu,maxnu))

        # #now write
        # self.writecontinuum(fil)
        # fil.write("intensity {0}, range {1} to {2} Ryd\n".format(integral[0],minry,maxry))

        #stop/save
        fil.write("##StopSave\n")
        fil.write("stop neutral column density {0}\n".format(self.column))
        #force no stop temperature
        fil.write("stop temperature off\n")
        #Add stop for CIE
        if(self.temperature):
            fil.write("stop column density 24\n")
        fil.write("iterate to convergence\n")
        fil.write('set save prefix "{0}"\n'.format(self.root))
        #fil.write('save last continuum units Angstrom ".con"\n')
        fil.write('save last column densities ".col"\n')
        fil.write('save last ionization means ".ion"\n')
        #fil.write('save last overview ".ovr"\n')
        fil.write('save last abundances ".abn"\n')
        #fil.write('save last radius ".rad"\n')

        fil.close()

    def inituvb(self):

        """
        This function is a dispatcher of the different UVB models
        currently supported by the code

        """

        #Initialize the HM12 UVB at the given redshift 
        if(self.uvbtype == 'HM12'):
            self.initcuba()

        #Initialize the HM05 UVB by printing the
        #  command to the Cloudy *.in file
        elif(self.uvbtype == 'HM05'):
            self.initHM05()

        #Initialize the HM05 UVB that creates the UVB at
        #  a given redshift, used with getionisation.py
        elif(self.uvbtype == 'HM05logU'):
            self.initHM05logU()

        #Initialize the FC09 UVB at the given redshift 
        elif(self.uvbtype == 'FG09'):
            self.initcafg()

        else:
            print('Unknown UVB!')
            exit()


    def initcuba(self):

        """
        Initialise the CUBA HM12 UVB at a given redshift

        """

        #import numpy as np
        #from astropy import constants as const

        #load cuba
        input_spectrum = resource_filename("pyigm", "/data/euvb/UVB_HM12.dat")
        uvb=open(input_spectrum)
        flag=0

        #store the input cuba
        cuba_red=[]
        cuba_wave=[]

        #pick the flux at three redshifts near to the z of interest
        cuba_one=[]
        cuba_two=[]
        cuba_three=[]

        #loop cuba file 
        for line in uvb:
            if not "#" in line:

                #parse cuba redshift first
                if(flag == 0):
                    line=line.strip()
                    line=line.split(" ")
                    for l in line:
                        cuba_red.append(float(l))
                    cuba_red=np.array(cuba_red)
                    flag = flag+1


                #now read the left and right spectrum for input redshift
                else:

                    #wave
                    line=line.strip()
                    line=line.split(" ")
                    cuba_wave.append(float(line[0]))

                    #check edge of the grid
                    index=np.argmin(abs(self.redshift - cuba_red))


                    #if edge of the grid - stay at that redshift
                    if((index < 1) | (index > cuba_red.size-2)):
                        cuba_one.append(float(line[index+1]))
                        cuba_two.append(float(line[index+1]))
                        cuba_three.append(float(line[index+1]))
                        cuba_int_red=[cuba_red[index],cuba_red[index],cuba_red[index]]

                    else: 
                        #prepare for interpolation by loading z +/- 1 setp
                        cuba_one.append(float(line[index]))
                        cuba_two.append(float(line[index+1]))
                        cuba_three.append(float(line[index+2]))
                        cuba_int_red=[cuba_red[index-1],cuba_red[index],cuba_red[index+1]]

        uvb.close()

        #turn to arrays - interpolate in log space
        cuba_wave=np.array(cuba_wave)
        cuba_one=np.log10(np.array(cuba_one))
        cuba_two=np.log10(np.array(cuba_two))
        cuba_three=np.log10(np.array(cuba_three))
        cuba_int_red=np.array(cuba_int_red)
        cuba_int=np.zeros(cuba_one.size)

        #interpolate to give redshift
        for ii in range(cuba_wave.size):
            cuba_int[ii]=np.interp([self.redshift],cuba_int_red,[cuba_one[ii],cuba_two[ii],cuba_three[ii]])
            #print(cuba_wave[ii], cuba_int[ii])

        #cuba has the same wave twice at discontinuities.
        #cloudy does not like that, so perturbe these values 
        for ii in range(cuba_wave.size-1):
            if(cuba_wave[ii]==cuba_wave[ii+1]):
                cuba_wave[ii]=cuba_wave[ii]-0.005

        #import matplotlib.pyplot as plt
        #cuba_wave=np.log10(cuba_wave)
        #plt.plot(cuba_wave,cuba_one,color='red')
        #plt.plot(cuba_wave,cuba_two,color='green')
        #plt.plot(cuba_wave,cuba_three,color='blue')
        #plt.plot(cuba_wave,cuba_int,color='black')
        #plt.show()

        #store the interpolated UVB
        lgnu=const.c.to('cm/s')/(cuba_wave*1e-8)
        lgnu=np.log10(lgnu.value) #Log of Hz
        cuba_int=cuba_int+np.log10(4*np.pi) #get rid of /str and go to Log of erg/s/cm2/Hz
        self.source['lgnu']=lgnu
        self.source['lgfpiJ']=cuba_int


    def initHM05(self):

        self.HM05write="table HM05 redshift {0}\n".format(self.redshift)
        self.HM05write="table HM05, z={0}\n".format(self.redshift)


    def initHM05logU(self):

        """
        Initialise the HM05 UVB at a given redshift

        """

        #import numpy as np
        #from astropy import constants as const

        #load hm05cbw
        input_spectrum = resource_filename("pyigm", "/data/euvb/UVB_HM05.dat")
        uvb=open(input_spectrum)
        flag=0
        
        #store the input hm05cbw
        hm05cbw_red=[]
        hm05cbw_wave=[]
        
        #pick the flux at three redshifts near to the z of interest
        hm05cbw_one=[]
        hm05cbw_two=[]
        hm05cbw_three=[]
        
        #loop hm05cbw file 
        for line in uvb:
            if not "#" in line:
                
                #parse hm05cbw redshift first
                if(flag == 0):
                    line=line.strip()
                    line=line.split()
                    for l in line:
                        hm05cbw_red.append(float(l))
                    hm05cbw_red=np.array(hm05cbw_red)
                    flag = flag+1
                    
                    
                #now read the left and right spectrum for input redshift
                else:
                    
                    #wave
                    line=line.strip()
                    line=line.split()
                    hm05cbw_wave.append(float(line[0]))
                    
                    #check edge of the grid
                    index=np.argmin(abs(self.redshift - hm05cbw_red))
                    
                    
                    #if edge of the grid - stay at that redshift
                    if((index < 1) | (index > hm05cbw_red.size-2)):
                        hm05cbw_one.append(float(line[index+1]))
                        hm05cbw_two.append(float(line[index+1]))
                        hm05cbw_three.append(float(line[index+1]))
                        hm05cbw_int_red=[hm05cbw_red[index],hm05cbw_red[index],hm05cbw_red[index]]
                        
                    else: 
                        #prepare for interpolation by loading z +/- 1 setp
                        hm05cbw_one.append(float(line[index]))
                        hm05cbw_two.append(float(line[index+1]))
                        hm05cbw_three.append(float(line[index+2]))
                        hm05cbw_int_red=[hm05cbw_red[index-1],hm05cbw_red[index],hm05cbw_red[index+1]]
                        
        uvb.close()
        
        #turn to arrays - interpolate in log space
        hm05cbw_wave=np.array(hm05cbw_wave)
        hm05cbw_one=np.log10(np.array(hm05cbw_one))
        hm05cbw_two=np.log10(np.array(hm05cbw_two))
        hm05cbw_three=np.log10(np.array(hm05cbw_three))
        hm05cbw_int_red=np.array(hm05cbw_int_red)
        hm05cbw_int=np.zeros(hm05cbw_one.size)
        
        #interpolate to give redshift
        for ii in range(hm05cbw_wave.size):
            hm05cbw_int[ii]=np.interp([self.redshift],hm05cbw_int_red,[hm05cbw_one[ii],hm05cbw_two[ii],hm05cbw_three[ii]])
            #print(hm05cbw_wave[ii], hm05cbw_int[ii])
            
        #hm05cbw has the same wave twice at discontinuities.
        #cloudy does not like that, so perturbe these values 
        for ii in range(hm05cbw_wave.size-1):
            if(hm05cbw_wave[ii]==hm05cbw_wave[ii+1]):
                hm05cbw_wave[ii]=hm05cbw_wave[ii]-0.005
                
        #import matplotlib.pyplot as plt
        #hm05cbw_wave=np.log10(hm05cbw_wave)
        #plt.plot(hm05cbw_wave,hm05cbw_one,color='red')
        #plt.plot(hm05cbw_wave,hm05cbw_two,color='green')
        #plt.plot(hm05cbw_wave,hm05cbw_three,color='blue')
        #plt.plot(hm05cbw_wave,hm05cbw_int,color='black')
        #plt.show()
        
        #store the interpolated UVB
        lgnu=const.c.to('cm/s')/(hm05cbw_wave*1e-8)
        lgnu=np.log10(lgnu.value) #Log of Hz
        hm05cbw_int=hm05cbw_int+np.log10(4*np.pi) #get rid of /str and go to Log of erg/s/cm2/Hz
        self.source['lgnu']=lgnu
        self.source['lgfpiJ']=hm05cbw_int


    def initcafg(self):

        """
        Initialise the FG09 UVB at a given redshift

        """

        #all available redshift for UVB 
        avail_z=np.arange(0,10,0.05)

        #pick closest in redshift
        zindex=np.argmin(abs(avail_z-self.redshift))
        nameuvb='fg_uvb_dec11/fg_uvb_dec11_z_{}.dat'.format(avail_z[zindex])
        try:
            fgtab=np.loadtxt(nameuvb,dtype={'names': ('nu','f'),'formats':('f10','f10')})
        except:
            print('Could not read UVB file {}'.format(nameuvb))
            exit()

        #from Ry to Hz
        self.source['lgnu']=np.log10(fgtab['nu'])+np.log10(3.2898419499e15)

        #from (10^-21 erg s^-1 cm^-2 Hz^-1 sr^-1) to Log of erg/s/cm2/Hz
        self.source['lgfpiJ']=np.log10(fgtab['f'])-21+np.log10(4*np.pi)

        #now reverse freq to have high freq first and low freq last (as in HM12)
        #This is nothing special, but it is needed for how the writing procedure
        #fixes the boundary of the SED
        self.source['lgnu'] = self.source['lgnu'][::-1]
        self.source['lgfpiJ'] = self.source['lgfpiJ'][::-1]

    def writecontinuum(self,fil):

        """ 
        Parse the continuum for an input cloudy file

        """

        #import numpy as np

        #make copies     
        spec = 10**np.array(self.source['lgfpiJ'])
        specclean = np.copy(spec)
        logfreq=np.copy(self.source['lgnu'])

        #first point
        specclean[spec == 0.0] = np.amin(spec[spec > 0])*1e-6
        logL_nu = np.log10(specclean)

        #put a floor at -40 when writing to avoid cloudy complaints 
        fil.write("interpolate")
        fil.write(" ({0:f} {1:f})".format(7.51, 
                                          max(np.amin(logL_nu)-6,-40)))
        fil.write(" ({0:f} {1:f})".format(logfreq[-1]-0.01, 
                                          max(np.amin(logL_nu)-6,-40)))
        #stuff in between 
        for i in range(len(logL_nu)):
            fil.write("\ncontinue")
            fil.write(" ({0:f} {1:f})".format(logfreq[-i-1],
                                              max(logL_nu[-i-1],-40)))
        #last point 
        fil.write("\ncontinue ({0:f} {1:f})".
                  format(logfreq[0]+0.01, max(np.amin(logL_nu)-6,-40)))

        fil.write(" ({0:f} {1:f})\n".format(22.4,max(np.amin(logL_nu)-6,-40)))

    def checkrun(self,stop=False):

        """
        Check if a model has run successfully
        Check stop check if the stop was due to lower column

        """

        #import os

        #check if file exists 
        status = False
        warning = True
        filname=self.path+'/'+self.root+'.out'
        filthere=os.path.isfile(filname)

        #check exist status
        linestop=''
        if filthere:
            status=True
            #print('Parse output for warnings')
            file1=open(filname,"r")
            #loop entire file
            for line in file1:
                #check exit condition
                if ("Cloudy exited OK" in line):
                    warning = False
                #grab stopping line
                if "Calculation stopped because" in line:
                    linestop=line
            file1.close()

            #if exit condition not clean, raise warning
            if warning:
                print("Warning: cloudy issues for {}".format(self.root))

            #check if stopped because of unwanted reasons:
            if stop:
                if ("H0 column dens reached" not in linestop) and ("H0-H0+H2 column dens reached" not in linestop):
                    print("Stop: condition not met in {0}".format(self.root))

        return status 

    def getion(self):

        """
        Get the ionization fractions from a model

        """
        #import re
        #import numpy as np

        #first check if everything is fine
        status=self.checkrun()

        #if so proceed
        if status:

            #open safe
            filname=self.path+'/'+self.root+'.ion'
            try:
                file1=open(filname,"r")
            except:
                print("Ion not found")
                return

            #now parse
            for line in file1:

                #check if a good data 
                if '-' in line:

                    #split it
                    fields=re.split(' |-|\n',line)
                    #remove unwanted and make it array numbers
                    fields=filter(None, fields)
                    values=np.array([-1*float(x) for x in fields[1:] if '.' in x])
                    self.ins[fields[0]]=values

            #close and return
            file1.close()

        else:
            print("Check model... Something is wrong!")
            return

    def getcolumn(self):

        """
        Get the column density from a model

        """

        #import re
        #import numpy as np

        #first check if everything is fine
        status=self.checkrun()

        #if so proceed
        if status:

            #open safe
            filname=self.path+'/'+self.root+'.col'
            try:
                file1=open(filname,"r")
            except:
                print("Column not found")
                return

            #now parse
            header=file1.readline()
            value=file1.readline()

            #extract
            clmd=np.array([float(x) for x in value.split()])
            header=np.array((header.split())[2:])

            #parity check 
            if (clmd.size != header.size):
                print("Something wrong in parser column")
                return

            for ii in range(header.size):
                if(clmd[ii] > 0):
                    self.clm[header[ii]]=np.log10(clmd[ii])

            #close and return
            file1.close()

        else:
            print("Check model... Something is wrong!")
            return


    def getprofile(self):

        """
        Get the radial overview from a model

        """

        #import re
        #import numpy as np

        #first check if everything is fine
        status=self.checkrun()

        #if so proceed
        if status:

            #open safe
            filname=self.path+'/'+self.root+'.ovr'
            try:
                file1=open(filname,"r")
            except:
                print("Overview not found")
                return

            #parse header
            header=file1.readline()
            header=header.replace("#","")
            header=header.split()

            #init dictionary
            for ff in header:
                self.prof[ff]=np.empty(0)

            #now parse
            for line in file1:
                line=line.split()

                ii=0
                for ff in header:
                    self.prof[ff]=np.append(self.prof[ff],float(line[ii]))
                    ii=ii+1

            #close and return
            file1.close()

        else:
            print("Check model... Something is wrong!")
            return


    def getabundance(self):

        """
        Get the gas phase abundance in each zone

        """

        #import re
        #import numpy as np

        #first check if everything is fine
        status=self.checkrun()

        #if so proceed
        if status:

            #open safe
            filname=self.path+'/'+self.root+'.abn'
            try:
                file1=open(filname,"r")
            except:
                print("Abundance not found")
                return

            #parse header
            header=file1.readline()
            header=header.replace("#","")
            header=header.split()
            header=header[1:]

            #init dictionary
            for ff in header:
                self.abund[ff]=np.empty(0)

            #now parse
            for line in file1:
                line=line.split()

                ii=0
                for ff in header:
                    self.abund[ff]=np.append(self.abund[ff],float(line[ii]))
                    ii=ii+1

            #close and return
            file1.close()

        else:
            print("Check model... Something is wrong!")
            return
