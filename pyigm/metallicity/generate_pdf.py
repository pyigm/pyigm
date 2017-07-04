
"""
Generate a PDF of a given quantity also accounting for boostrap errors

"""
from __future__ import print_function, absolute_import, division, unicode_literals

#some global staff first
import glob
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.rcParams.update({'font.size': 16})
import os.path
import h5py

from astropy.table import Table

from xastropy.xutils import xdebug as xdb

mblue='#1E90FF'
mred='#DC143C'

#now the work function
def work_generate_pdf(modelsdir, quantity, data_fil, outh5):

    print('Processing ', modelsdir, quantity)

    if(quantity == 'met'):
        minval=-5
        maxval=1.5
        binsize=0.15
        precomp=False
    elif(quantity == 'dens'):
        minval=-5
        maxval=1
        binsize=0.15
        precomp=False
    elif(quantity == 'jgal'):
        minval=-5.4
        maxval=2.4
        binsize=0.15
        precomp=False
    elif(quantity == 'jqso'):
        minval=-6.8
        maxval=1.0
        binsize=0.15
        precomp=False
    elif(quantity == 'fstar'):
        minval=-2
        maxval=2
        binsize=0.15
        precomp=False
    elif(quantity == 'temp'):
        minval=3.6
        maxval=6.4
        binsize=0.15
        precomp=False
    elif(quantity == 'upar'):
        precomp=True
        filein=modelsdir+'/'+modelsdir+'_ionisation.pkl'
        if(os.path.isfile(filein)):
            pdfcom=pickle.load(open(filein))
        else:
            print('Need to compute {} first'.format(filein))
            exit()
        #load precomputed histograms
        binsize=pdfcom['U_base'][1]-pdfcom['U_base'][0]
        minval=np.min(pdfcom['U_base'])
        maxval=np.max(pdfcom['U_base'])+binsize
        base_hist=pdfcom['U_base']
        pdf_hist=pdfcom['U_hist']
    elif(quantity == 'xhi'):
        precomp=True
        filein=modelsdir+'/'+modelsdir+'_ionisation.pkl'
        if(os.path.isfile(filein)):
            pdfcom=pickle.load(open(filein))
        else:
            print('Need to compute {} first'.format(filein))
            exit()
        #load precomputed histograms
        binsize=pdfcom['XHI_base'][1]-pdfcom['XHI_base'][0]
        minval=np.min(pdfcom['XHI_base'])
        maxval=np.max(pdfcom['XHI_base'])+binsize
        base_hist=pdfcom['XHI_base']
        pdf_hist=pdfcom['XHI_hist']
    elif(quantity == 'si2'):
        precomp=True
        filein=modelsdir+'/'+modelsdir+'_ionisation.pkl'
        if(os.path.isfile(filein)):
            pdfcom=pickle.load(open(filein))
        else:
            print('Need to compute {} first'.format(filein))
            exit()
        #load precomputed histograms
        binsize=pdfcom['XSiII_base'][1]-pdfcom['XSiII_base'][0]
        minval=np.min(pdfcom['XSiII_base'])
        maxval=np.max(pdfcom['XSiII_base'])+binsize
        base_hist=pdfcom['XSiII_base']
        pdf_hist=pdfcom['XSiII_hist']
    else:
        print('Option not found!')
        exit()
    
    #some util 
    nsample=1000
    figname=modelsdir

    #load all data
    info = Table.read(data_fil, format='ascii')
    #info=np.loadtxt(data_fil,dtype={'names':('name','zabs','ezabs','ion','logN','elogN','flag','sample'),
    #                           'formats':('S30', 'f8', 'f8','S4','f6','f6','i4','S5')})

    #if quantity is not precomputed, compute it
    if(precomp == False):
        #load the models
        listin=glob.glob(modelsdir+'/*emcee*hd5')
        nsys=len(listin)
        print('Found {} fits'.format(nsys))

        #This defines the bin edges, including the rightmost edge
        base_hist=np.arange(minval,maxval,binsize)
        pdf_hist=[]
        sample_tag=[]
        outh5[quantity]['left_edge_bins'] = base_hist[0:-1]
        outh5[quantity]['right_edge_bins'] = base_hist[1:]
        outh5[quantity]['left_edge_bins'].attrs['BINSIZE'] = binsize

        #first, loop over the sample to construct PDFs for individual chains
        print('Generating PDFs for individual systems')

        for fit in listin:

            # Load
            fh5=h5py.File(fit, 'r')
            #xdb.set_trace()
            #fl=open(fit)
            #result=pickle.load(fl)
            #fl.close()
    
            #find sample tag
            #match=np.where(info['name'] == result['info']['name'])
            match=np.where(info['name'] == fh5['inputs'].attrs['name'])
            sample_tag.append(info['sample'][match])
            
            #populate histo
            #qnt=result['tags'].index(quantity)
            qnt= np.where(fh5['outputs']['tags'].value == quantity)[0][0]
            #hist,edge=np.histogram(result['pdfs'][:,qnt],bins=base_hist)
            hist,edge=np.histogram(fh5['outputs']['pdfs'][:,qnt],bins=base_hist)

            #normalise 
            hist=hist/binsize/len(fh5['outputs']['pdfs'][:,qnt])
            
            #Append
            pdf_hist.append(hist)

            # Write to hdf5
            outh5[quantity][fh5['inputs'].attrs['name']] = hist
            # Inputs
            try:
                outh5['inputs'].create_group(fh5['inputs'].attrs['name'])
            except ValueError:  # No repeat
                pass
            else:
                for key in fh5['inputs']:
                    outh5['inputs'][fh5['inputs'].attrs['name']][key] = fh5['inputs'][key].value
    else:
        print("Work with precomputed quantity")
        nsys=len(pdfcom['name'])
        print('Found {} fits'.format(nsys))

        #Now need to find the sample tag
        sample_tag=[]

        for fit in pdfcom['name']:

            #find sample tag
            match=np.where(info['name'] == fit)
            sample_tag.append(info['sample'][match][0])

    #now loop over trials for entire sample
    pdf_trials=np.zeros((nsample,base_hist.size-1))

    print('Bootstrap all sample')

    for trial in range(nsample):

        #now generate random indexes with repetition
        thisindex=np.random.random_integers(0,high=nsys-1,size=nsys)

        #fill the pdf
        pdf_thisset=np.zeros(base_hist.size-1)

        for ii in thisindex:
            pdf_thisset=pdf_thisset+pdf_hist[ii]
        
        #normalise and stack
        pdf_thisset=pdf_thisset/nsys
        pdf_trials[trial,:]=pdf_thisset

    #now compute percentiles
    pdf_med=np.percentile(pdf_trials,50,axis=0)
    pdf_90=np.percentile(pdf_trials,90,axis=0)
    pdf_75=np.percentile(pdf_trials,75,axis=0)
    pdf_25=np.percentile(pdf_trials,25,axis=0)
    pdf_10=np.percentile(pdf_trials,10,axis=0)

    #plot
    plt.plot(base_hist[0:-1],pdf_med)
    plt.plot(base_hist[0:-1],pdf_90)
    plt.plot(base_hist[0:-1],pdf_10)
    plt.show()

    #now extract the statistical subset
    sample_tag=np.array(sample_tag)
    substat=np.where(sample_tag == 'P15')
    substat=substat[0]
    nstat=len(substat)
    pdf_hist_sub=[]

    for ii in substat:
        pdf_hist_sub.append(pdf_hist[ii])
    
    #now repeat sampling for the statistical subset
    pdf_trials=np.zeros((nsample,base_hist.size-1))

    return
    #####################################################################
    #####################################################################
    #####################################################################

    boot_stat = False
    if boot_stat:
        print('Bootstrap stat sample')

        for trial in range(nsample):

            #now generate random indexes with repetition
            try:
                thisindex=np.random.random_integers(0,high=nstat-1,size=nstat)
            except ValueError:
                xdb.set_trace()

            #fill the pdf
            pdf_thisset=np.zeros(base_hist.size-1)

            for ii in thisindex:
                pdf_thisset=pdf_thisset+pdf_hist_sub[ii]

            #normalise to nsys (not nstat to have fraction) and stack
            pdf_thisset=pdf_thisset/nsys
            pdf_trials[trial,:]=pdf_thisset

        #now compute percentiles
        pdf_med_stat=np.percentile(pdf_trials,50,axis=0)
        pdf_90_stat=np.percentile(pdf_trials,90,axis=0)
        pdf_75_stat=np.percentile(pdf_trials,75,axis=0)
        pdf_25_stat=np.percentile(pdf_trials,25,axis=0)
        pdf_10_stat=np.percentile(pdf_trials,10,axis=0)


        #write
        nameout="{0}/{1}_{2}_pdfdata.pkl".format(modelsdir,figname,quantity)
        save={'base':base_hist,'nsample':nsample,'minval':minval,'maxval':maxval,'binsize':binsize,'quantity':quantity,
              'pdf50':pdf_med,'pdf90':pdf_90,'pdf75':pdf_75,'pdf25':pdf_25,'pdf10':pdf_10,'nsys':nsys,'nstat':nstat,
              'pdf50stat':pdf_med_stat,'pdf90stat':pdf_90_stat,'pdf75stat':pdf_75_stat,'pdf25stat':pdf_25_stat,'pdf10stat':pdf_10_stat}
        pickle.dump(save,open(nameout,'w'))

        #plot
        plt.plot(base_hist[0:-1],pdf_med_stat)
        plt.plot(base_hist[0:-1],pdf_90_stat)
        plt.plot(base_hist[0:-1],pdf_10_stat)
        plt.show()

if __name__ == "__main__":
    
    #start with general quantities
    workers=[]
    #allmod=['minimalnhi','sourcenhi','cafgnhi','dustnhi']
    #allmod=['MCMC_TST']
    allmod=['./Lehner13-MCMC-CRC/MCMC_FULL/']
    alltype=['met', 'dens']
    mcmc_fil = './Lehner13-MCMC-CRC/MgIIsurvey_MTL.ascii'
    hdf5_outfil = './MgIIsurvey_MTL.hdf5'

    #allmod=['dustnhi']
    #alltype=['fstar']
   
    #allmod=['dustnhi']
    #alltype=['upar','xhi','si2']
    #alltype=['xhi']

    outh5 = h5py.File(hdf5_outfil, 'w')
    outh5.create_group('inputs')

    for mm in allmod:
        for tt in alltype:
            outh5.create_group(tt)
            work_generate_pdf(mm, tt, mcmc_fil, outh5)

    print("Writing {:s}".format(hdf5_outfil))
    outh5.close()
