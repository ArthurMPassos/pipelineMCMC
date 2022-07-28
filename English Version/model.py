__author__ = "Yuri Neto, Felipe Pinho, Beatriz Duque, Adriana Valio (Orientadora)"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "starsandexoplanets@gmail.com"
__status__ = "Production"

#%matplotlib nbagg
#%matplotlib inline

import scipy
from scipy import interpolate
from scipy import ndimage
import scipy.signal as sps
from scipy.optimize import minimize
from PyAstronomy.pyasl import binningx0dt
import collections
import glob
#from astropy.io import fits
import pandas as pd
import matplotlib.patches as mpatches
import lmfit
from lmfit import Model
import emcee
import corner

import requests
from logging import Handler, Formatter
import logging
import datetime
from pytz import timezone

import timeit

from IPython.display import display, Math
from IPython.display import display, Math

from estrela_nv1 import estrela #star and eclipse:: extensions of auxiliary programs that perform the calculation of the light curve.
from eclipse_nv1 import Eclipse
from verify import Validar, calSemiEixo, converte

import numpy
import matplotlib as mpl
from matplotlib import pyplot
import matplotlib.pyplot as plt

import lightkurve as lk
from lightkurve import search_lightcurvefile


class Modelo:

    def __init__(self,estrela, eclipse):
        '''
         star parameter :: star class
         eclipse parameter::eclipse class
         '''
         #collecting star objects
        self.u1 = estrela.getu1()
        self.u2 = estrela.getu2()
        self.n = estrela.getTamanhoMatriz()
        self.r = estrela.getRaioStar() #star radius in pixel
        self.r_Sun = estrela.getRaioSun() #star radius in Rsun
        self.mx = estrela.getIntensidadeMaxima()
        self.star_name = estrela.getStarName() #star name
        self.cadence = estrela.getCadence() #star cadence (short or long)

        #collecting Eclipse objects
        self.raioPlan = eclipse.getRaioPlan()
        self.R_jup = eclipse.getRplanJup()
        self.AU = eclipse.getSemiEixo() #axle shaft in AU
        self.semiEixoRaioStar = eclipse.getsemiEixoRaioStar() #semi-axis with respect to the radius of the star
        self.porb = eclipse.getPeriodo()
        self.inc = eclipse.getInc()
        self.ecc,self.anom = eclipse.getEccAnom()
        self.lua = eclipse.getLua()

        #variables that will be returned from the rd_data function
        self.time = 0.
        self.flux = 0.
        self.flux_err = 0.


    def rd_data(self,plot,save_data):
                
        ''' 
        Function created to access starlight curve data and extract time and flux.
        
        lightkurve packeage:: used to open the Kepler Light curves and
        also open TESS light curve/ general light curves
        Documentation:: https://docs.lightkurve.org/api/index.html documentation of using lightkurve
        
        expected outputs::time, flux
        
        plot = 1 (plots the light curve, to not plot just enter any value)
        save_data = 1 plots the light curve (to not plot, just enter any value)
        '''
    ##--------------------------------------------------------------------------------------------------------------------------------------------------##
    # PDCSAP_FLUX is used because the analysis will be carried out in transit.
        lc = search_lightcurvefile(self.star_name, cadence = self.cadence).download_all().SAP_FLUX
        time = [] # time = array with the time data
        flux = [] # flux = array with the flux data
        flux_err = [] # flux_err = array with the flow error data
        time_temp = [] #temporary variables 
        flux_temp = []
        flux_err_temp = []

        for i in range(0, len(lc)-1):
            flux_temp.append(lc[i].flux)
            flux_err_temp.append(lc[i].flux_err)
            time_temp.append(lc[i].time.to_value('bkjd','float'))

        for i in range(0, len(lc)-1):   
            
    # Limit indices to read time values
            flux_temp[i] = flux_temp[i][0:time_temp[i].size]
            flux_err_temp[i] = flux_err_temp[i][0:time_temp[i].size]
        
    # Deletes all NaN (Not a Number) values from the data
            time_temp[i] = time_temp[i][~numpy.isnan(flux_err_temp[i])]
            flux_temp[i] = flux_temp[i][~numpy.isnan(flux_err_temp[i])]
            flux_err_temp[i] = flux_err_temp[i][~numpy.isnan(flux_err_temp[i])]
        
    # Normalize each quarter
            flux_err_temp[i] = flux_err_temp[i]/ abs(numpy.median(flux_temp[i]))
            flux_temp[i] = flux_temp[i]/ abs(numpy.median(flux_temp[i]))

        for i in range(0, len(lc)-1):
            flux = numpy.append(flux, flux_temp[i])
            time = numpy.append(time, time_temp[i])
            flux_err = numpy.append(flux_err,flux_err_temp[i])
        
    #plot of the complete light curve
        if plot == (1):
            plt.rcParams['figure.figsize'] = 10,4
            graf1, ax = plt.subplots()

            ax.set_xlabel('Time (BJD - 2454833)')
            ax.set_ylabel('Normalized Flux')

            ax.set_title("Light Curve - " + self.star_name)
            ax.set_xlim(min(time), max(time))
            ax.set_ylim(min(flux), max(flux))
        
            #ax.plot(self.time, self.flux, "k.", ms=2)
            
            ax.errorbar(time, flux, yerr = flux_err, fmt = '.k', capsize = 0,alpha = 0.5)
            
    #save the data to a .dat file
        if save_data == 1:
            numpy.savetxt('%s_LC.dat'%self.star_name, numpy.c_[(time, flux, flux_err)])

        self.time = time 
        self.flux = flux
        self.flux_err = flux_err

        return self.time, self.flux, self.flux_err

    def det_x0(self, plot):
        
        '''
        Function to get the center of the first transit.
        porb is used to create a phased light curve a smooth is applied to the phased light curve.
        time parameter::
        flux parameter::
        porb parameter:: orbital period of the planet (per in days)
        plot parameter::
        
        returns
         x0 = value of the center of the first transit
         nt = number of possible transits
         plot = 1 (plots the light curve, to not plot just type any value)
        '''
        time = self.time 
        porb = self.porb
        flux = self.flux
        
        phase = (time % porb)/ porb
        jj = numpy.argsort(phase)
        ff = phase[jj]

        smoothed_LC = scipy.ndimage.filters.uniform_filter(flux[jj], size = 100) # equivalent to idl's smooth with edge_truncade
        smoothed_LC[0:200] = 1
        smoothed_LC[len(flux[jj])-200:len(flux[jj])] = 1
        x = phase[jj]
        y = 1 - smoothed_LC
        yh = 0.002

        kk = numpy.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # central value of phased transits

        self.x0 = (numpy.fix(time[0] / porb) + fa0) * porb # first transit center time
        
        self.nt = numpy.fix(max(time - time[0])/ porb) + 1 # number of possible transits

        # plot of the complete light curve
        if plot == 1:
        
            plt.rcParams['figure.figsize'] = 10,8
            graf2, ax = plt.subplots(2, 1, sharex = False, gridspec_kw = {'height_ratios': [1, 1]})
            graf2.subplots_adjust(hspace = 0.5)
        
            ax[0].set_ylim(0.9, 1.1)
            ax[0].set_title("Phased LC")
        
            ax[0].set_xlabel('Time')
            ax[0].set_ylabel('Normalized Flux')

            ax[0].plot(phase, flux, "k.", ms = 2)
            ax[0].plot(phase[jj], smoothed_LC, "r.", ms = 2)        
        
            ax[1].set_ylim(0.9, 1.1)
            ax[1].set_title("x0 = first transit center time")
        
            ax[1].set_xlabel('Time (BJD - 2454833)')
            ax[1].set_ylabel('Normalized Flux')
        
            ax[1].set_xlim(min(time) - 2, min(time) + 2)
            ax[1].set_ylim(0.95, 1.01)
        
            ax[1].plot(time, flux, "k.", ms=2)

            plt.axvline(x = self.x0,linewidth = 2, color='r')
        
        return self.x0, self.nt

    #--------------------------------------------------#
    def limb(self, plot):
        '''
        Function that generates a synthesized star with
        limb obscuration given by 1-u1*(1-cos(mu))-u2*(1-cos(mu))^2), where mu=heliocentric angle
        
        -- limbo darkening coefficient --
        parameter u1 :: linear coefficient
        parameter u2 :: coefficient of the quadratic term
        
        returns
         parameter wl:: array with star intensity
         if plot = 1, plot the star profile (to not plot, just enter any value)
        '''
        #PERGUNTAR

        if self.u1 == 999:
            self.u1 = 0.59
        if self.u2 == 999:
            self.u2 = 0.0
        
        self.wl = numpy.zeros((self.n, self.n))
        x = numpy.arange(self.n/2-self.n, self.n - self.n/2, 1)
        
        for j in range (0, self.n):     
            z = numpy.sqrt(x**2 + (j - self.n/2.)**2)
            kk = numpy.where(z <= self.r)
            if kk[0].size > 0:
                    m = numpy.cos(numpy.arcsin(z[kk]/self.r))
                    self.wl[kk, j] = self.mx*(1 - self.u1*(1 - m) - self.u2*(1 - m)**2)
        
        data = self.wl[:,int(self.n/2)]
        smoothed = scipy.ndimage.filters.uniform_filter(data, size = 5)

        if plot == 1:
            plt.rcParams['figure.figsize']= 12, 5
            graf1,ax = plt.subplots()
            ax.plot(x,smoothed)
            ax.plot(x, self.wl[:, 428],'k')
            ax.set_title('limb')


        return self.wl, self.u1, self.u2 
    #--------------------------------------------------#

    def eclipse_model(self):
        '''
        Call of auxiliary programs for the creation of the light curve model, which may contain:
         - One planet or more
         - One stain or more
         - A moon or more
        
         parameter u1 :: limbo darkening coefficient 1
         parameter u2 :: limbo darkening coefficient 1
         parameter per:: transit period in days
         parameter a :: semi-axis in UA
         parameter inc:: tilt angle in degrees
         parameter rp :: radius of the planet in relation to the radius of the star
        
         returns
         lc_model parameter:: light curve
         ts_model parameter:: transit time in Hours
        '''
        estrela_1 = estrela(self.r,self.r_Sun, self.mx , self.u1, self.u2, self.n)  # create the star object
        Nx1 = estrela_1.getNx() # collect parameters from the star array
        Ny1 = estrela_1.getNy()
        raioEstrelaPixel1 = estrela_1.getRaioStar() # collect star ray in pixel
        estrelaManchada1 = estrela_1.getEstrela() # collect spotted star

        eclipse1 = Eclipse(Nx1, Ny1, raioEstrelaPixel1, estrelaManchada1)  #create the eclipse object

        eclipse1.setTempoHoras(1.)

        eclipse1.criarEclipse(self.semiEixoRaioStar, self.AU, self.raioPlan,self.R_jup, self.porb,self.inc,self.lua,self.ecc,self.anom, False)

        self.lc_model = numpy.array(eclipse1.getCurvaLuz())
        self.ts_model = numpy.array(eclipse1.getTempoHoras())
        
        return self.lc_model, self.ts_model

    #--------------------------------------------------#
    def retornaParametros(self):
        return self.u1,self.u2,self.porb,self.time,self.flux,self.flux_err,self.raioPlan,self.AU,self.inc, self.x0, self.nt,self.ts_model
    
    def setTime(self,time):
        self.time = time

    def setFlux(self,flux):
        self.flux = flux
    
    def setFluxErr(self,flux_err):
        self.flux_err = flux_err


class Tratamento :

    def __init__(self,modelo):

        '''
        Function to extract the traffics individually from the light curve
        
        time parameter :: total light curve time
        flux parameter :: total light curve flux
        flux_err parameter:: flux error
        parameter u1 :: limbo darkening coefficient 1
        parameter u2 :: limbo 2 darkening coefficient
        parameter porb :: orbit period in days
        parameter AU :: orbital semi-axis in UA
        parameter radiusPlan:: radius of the planet in relation to the radius of the star
        parameter inc:: tilt angle in degrees
        parameter x0 ::
        parameter nt ::
        '''
        self.modelo = modelo
        self.u1,self.u2,self.porb,self.time,self.flux,self.flux_err,self.raioPlan,self.AU,self.inc,self.x0,self.nt,self.ts_model = modelo.retornaParametros()

    def cut_transit_single(self):
        
        '''
        returns
        
        parameter dur :: duration of transit in hours
        parameter t_split :: time in hours (same for all transits)
        parameter n_f_split ::normalized traffic light curve
        parameter n_f_err_split :: normalized flow error
        '''

        if self.u1 == 999:
            self.u1 = 0.59
        if self.u2 == 999:
            self.u2 = 0.0
        
        self.wl, self.u1, self.u2 = self.modelo.limb(0)    
        lc0, ts0 = self.modelo.eclipse_model()
        
        # duration of transit in hours
        
        x = ts0
        y = 1 - lc0
        #yh = 0.002
        yh=max(y)/2.
        kk = numpy.where(y >= yh)
        x1 = min(x[kk])
        x2 = max(x[kk])
        meio = (x1 + x2)/ 2 # center value of phased transits
        self.dur = 2 * numpy.abs(ts0[min(numpy.where(lc0 < 1))[0]]) # total transit time
        # in hours
        
        # given the central value of each transit, I determine points (+-) at a distance of 45%
        # of the period. 
        mm = []
        mp = []
        ttt = []

        for i in range(0, int(self.nt)):
            mm.append(self.x0 + self.porb*i - (self.porb*0.45))
            mp.append(self.x0 + self.porb*i + (self.porb*0.45))
            
        # index what each value represents in the time array
        for i in range(0, int(self.nt)):
            ttt.append(numpy.where((self.time >= mm[i]) & (self.time <= mp[i])))
        
        # separate the time and flow arrays with these values.
        self.t_split = []
        f_split = []
        f_err_split = []
        for i in range(0, int(self.nt)):
            self.t_split.append((self.time[ttt[i]]-self.x0-(self.porb*i))*24)   # time in hours with the transit center = 0
            f_split.append(self.flux[ttt[i]])         
            f_err_split.append(self.flux_err[ttt[i]])  
            
        # slope elimination
        
        size = [] # number of points present in up to 3* the transit duration
        for u in range(int(self.nt)):
            size.append(len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]))  
            
        # considered data that present a flow greater than 0
        # and time data with amount of points with 0.9*numpy.mean(size)
        self.n_f_split = []
        n_f_err_split = []
        for i in range(0, int(self.nt)):
            if len((f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                ss = numpy.polyfit(self.t_split[i], f_split[i], 1)
                zz = numpy.polyval(ss, self.t_split[i])
                self.n_f_split.append(numpy.array(f_split[i] - zz + 1))
                n_f_err_split.append(f_err_split[i])
            else:
                self.n_f_split.append(f_split[i])
                n_f_err_split.append(f_err_split[i])
                
        # renormalization
        w_flux = []
        for i in range(0, int(self.nt)):
            if len((self.n_f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                w_flux =  numpy.append(w_flux, self.n_f_split[i])
        m0 = numpy.median(w_flux) 
        for i in range(0, int(self.nt)):
            if len((self.n_f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                self.n_f_split[i] = self.n_f_split[i]/m0
                n_f_err_split[i] = n_f_err_split[i]/m0

        return self.dur, self.t_split, self.n_f_split, n_f_err_split
        #t_split = tim
        #n_f_split = lcurve
    #--------------------------------------------------#

    def transit_smooth(self,ntransit, selection):
    
        '''
        Function for a smooth curve with n transits
        
        parameter ntransit :: number of transits to use in the smoothed curve
        parameter selection :: 0, use a random choice of all transits
        parameter deepest_transit :: 1, use the deepest transits
        
        returns
         parameter time_phased[bb] :: time
         smoothed_LC[bb] parameter:: Smoothed light curve
        '''
        
        if selection == 0:
            tran_selec = numpy.random.randint(int(self.nt), size=(1, ntransit))[0]
        
        else:
            deepest_transit = []
            for i in range(0, int(self.nt)):
                if len(self.n_f_split[i]) > 0:
                    deepest_transit.append(numpy.mean(self.n_f_split[i]))
                else:
                    deepest_transit.append(900)
            tran_selec = numpy.argsort(deepest_transit)[0:ntransit]    
        
        lc = []
        t = []
        
        for i in tran_selec:
            lc = numpy.append(lc, self.n_f_split[i])
            t = numpy.append(t, (self.t_split[i]+self.porb*24*i)/24+self.x0)
        
        phase = (t % self.porb)/ self.porb
        jj = numpy.argsort(phase)
        ff = phase[jj]

        self.smoothed_LC = scipy.ndimage.filters.uniform_filter(lc[jj], size = 100) # equivalent to idl's smooth with edge_truncade

        x = phase[jj]
        y = 1 - self.smoothed_LC
        yh = 0.002

        kk = numpy.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # central value of phased transits

        self.time_phased = (ff - fa0)*self.porb*24

        bb = numpy.where((self.time_phased >= min(self.ts_model)) & (self.time_phased <= max(self.ts_model)))
        
        return self.time_phased[bb], self.smoothed_LC[bb]

        def gettime_phased(self):
            return self.time_phased


class Ajuste:

    def __init__(self,tratamento, time, flux, nwalkers, niter, burnin):

        self.u1_p0 = 0.5
        self.u2_p0 = 0.1
        self.a_p0 = 0.05
        self.inc_p0 = 88.
        self.rp_p0 = 1

        self.time = time
        self.flux = flux

        self.flux_err = numpy.var(self.flux)
        self.data = (self.time, self.flux, self.flux_err)

        self.nwalkers = nwalkers
        self.niter = niter
        self.burnin = burnin


        self.initial = numpy.array([self.u1_p0, self.u2_p0, self.a_p0, self.inc_p0, self.rp_p0])
        self.ndim = len(self.initial)

        self.p0 = [numpy.array(self.initial) + 1e-4 * numpy.random.randn(self.ndim) for i in range(self.nwalkers)]

        self.tratamento = tratamento

    #--------------------------------------------------#
    #----------------------MCMC------------------------#
    #--------------------------------------------------#
    def eclipse_mcmc(self, time, theta):
        rsun = 1.
        u1, u2, semiEixoUA, anguloInclinacao, raioPlanJup = theta
        periodo = 1.
        raioStar, raioPlanetaRstar, semiEixoRaioStar = converte(rsun,raioPlanJup,semiEixoUA)
        
        estrela_ = estrela(373, raioStar, 240., u1, u2, 856)
        Nx = estrela_.getNx()
        Ny = estrela_.getNy()
        raioEstrelaPixel = estrela_.getRaioStar()
        estrelaManchada = estrela_.getEstrela()
        
        eclipse = Eclipse(Nx,Ny,raioEstrelaPixel,estrelaManchada)
        eclipse.setTempoHoras(1.)
        eclipse.criarEclipse(semiEixoRaioStar, semiEixoUA, raioPlanetaRstar, raioPlanJup, periodo, anguloInclinacao, 0,0 , 0, False, False)
        lc0 = numpy.array(eclipse.getCurvaLuz()) 
        ts0 = numpy.array(eclipse.getTempoHoras()) 
        return interpolate.interp1d(ts0,lc0,fill_value="extrapolate")(time)
    #--------------------------------------------------#
    def lnlike(self, theta, time, flux, flux_err):
        return -0.5 * numpy.sum(((flux - self.eclipse_mcmc(time, theta))/flux_err) ** 2)
    #--------------------------------------------------#
    def lnprior(self, theta):
        u1, u2, semiEixoUA, anguloInclinacao, rp = theta
        if 0.0 < u1 < 1.0 and 0.0 < u2 < 1.0 and 0.001 < semiEixoUA < 1 and 80. < anguloInclinacao < 90 and 0.01 < rp < 5:
            return 0.0
        return -numpy.inf
    #--------------------------------------------------#
    def lnprob(self, theta, time, flux, flux_err):
        lp = self.lnprior(theta)
        if not numpy.isfinite(lp):
            return -numpy.inf
        return lp + self.lnlike(theta, time, flux, flux_err)
    #--------------------------------------------------#
    def main(self):
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=self.data)

        print("Running burn-in...")
        self.p0, _, _ = self.sampler.run_mcmc(self.p0, self.burnin, progress=True)
        self.sampler.reset()

        print("Running production...")    
        self.pos, self.prob, self.state = self.sampler.run_mcmc(self.p0, self.niter, progress=True)

        return self.sampler, self.pos, self.prob, self.state
    #--------------------------------------------------#
