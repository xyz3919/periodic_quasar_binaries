import numpy as np
from numpy import mean, median
from astropy.stats import sigma_clip
import scipy.signal as sig
#from lomb_scargle_red_fix import lomb
from astroML.time_series import lomb_scargle
import os, subprocess
from scipy.optimize import curve_fit
import emcee
import scipy.optimize as op


class quasar_drw:

    """ 
    Originally written by Wei-Ting Liao 2018
    Modified by Yu-Ching (Tony) Chen Oct. 2018
    This version directly fits state space 
    """
    
    def __init__(self, time, signal, error, redshift, preprocess=True):
        self.time     = np.array(time,   dtype=np.float64)
        self.signal   = np.array(signal, dtype=np.float64)
        self.error    = np.array(error,  dtype=np.float64)
        self.redshift = float(redshift)
        self._preprocessed = False
        
        if ( len(time) != len(signal) ) or ( len(time)!= len(error) ):
            print("[quasar_drw] Error in input data: time, signal, error must have the same length.")
        
        self._sort_data()
        if preprocess == True:
            self._preprocess()
        
        self.__initiate()
        

    
    def __initiate(self):
        ## parameters for periodogram
        if (len(self.time) >= 2 and float( np.max(self.time) - np.min(self.time) ) > 0.0 ):
            self.__Tspan    = float( np.max(self.time) - np.min(self.time) )
            self.__Ndata    = len(self.signal)
            self.__psd_freq = \
                np.linspace(1.0/self.__Tspan, self.__Ndata/(2.0*self.__Tspan), 10*self.__Ndata)
               # np.linspace(1.0/self.__Tspan, self.__Ndata/(2.0*self.__Tspan), self.__Ndata) 
            self.__dt = self.__Tspan / float(self.__Ndata)
            self.__df = self.__psd_freq[1] - self.__psd_freq[0]
        else:
            pass
        
    
    def get_Tspan(self):
        if (len(self.time) >= 2):
            return self.__Tspan
        else:
            pass
            
    def get_Ndata(self):
        if (len(self.time) >= 2):
            return self.__Ndata
        else:
            pass
            
    def get_psd_freq(self):
        return self.__psd_freq
        
    def get_psd_time(self):
        return 1.0/self.__psd_freq
        
    def get_psd_time_err(self):
        period = self.get_psd_time()
        return (self.__df * period**2.)/2.
            
    
    def _preprocess(self):
        #self._sort_data()
        self._no_outlier()
        self._bin_data()
        self._preprocessed = True
        self.__initiate()
        
    
    
    def get_lc(self):
        """ output: time, signal, error """
        return (self.time, self.signal, self.error)
        
        
    def get_redshift(self):
        return self.redshift
        
    
    def add_lc(self, time, signal, error, preprocess=True):
        self.time   = np.array(np.append(self.time,   time),   dtype=np.float64) 
        self.signal = np.array(np.append(self.signal, signal), dtype=np.float64)
        self.error  = np.array(np.append(self.error,  error),  dtype=np.float64)
        
        self._sort_data()
        
        self._preprocessed = False
        
        if (preprocess == True):
            self._preprocess()
            
        self.__initiate()
        
         
        
    def ls_astroML(self):
        """
        calculate periodogram using generalized Lomb-Scargle periodogram from AstroML
        function description: http://www.astroml.org/modules/generated/astroML.time_series.lomb_scargle.html
        example: http://www.astroml.org/book_figures/chapter10/fig_LS_example.html
        """
        LS_lc = lomb_scargle(self.time, self.signal, self.error, self.__psd_freq*(2.0*np.pi), generalized=True)
        
        return 1.0/self.__psd_freq, LS_lc

    
    ### ********************************* ###
    ###  helper functions for preprocess  ###
    ### ********************************* ###
        
    def _sort_data(self):
        
        # take away points w/o data
        idx = self.error > 0.
        time   = self.time[idx]
        signal = self.signal[idx]
        error  = self.error[idx]
        
        idx = self.time > 0.
        time   = self.time[idx]
        signal = self.signal[idx]
        error  = self.error[idx]
        
        
        # sort
        idx = np.argsort(time)
        time   = time[idx]
        signal = signal[idx]
        error  = error[idx]
        
        # restore data
        self.time   = time
        self.signal = signal
        self.error  = error
    
    
    def _no_outlier(self, sigma=5, iters=100):
        idx = ((np.abs(self.signal) < 100.) & (self.signal > 0.))
        self.time   = self.time[idx]
        self.signal = self.signal[idx]
        self.error  = self.error[idx]
        
        after_clip = sigma_clip(self.signal, sigma=sigma, iters=iters, cenfunc=median, copy=True)
        
        idx = ~(after_clip.mask)
        self.time   = self.time[idx]
        self.signal = self.signal[idx]
        self.error  = self.error[idx]

        
    def _bin_data(self):
        time2   = []
        signal2 = []
        error2  = []
        count   = 0
    
        while(count < len(self.time)):
            idx = ( np.floor(self.time) == np.floor(self.time[count]) )
            signal_temp = self.signal[idx]
            error_temp  = self.error[idx]
            nn          = len(signal_temp)
        
            signal_temp, error_temp = self.__mag2flux(signal_temp, error_temp)
            signal_temp, error_temp = self.__weighted_mean(signal_temp, error_temp)
            signal_temp, error_temp = self.__flux2mag(signal_temp, error_temp)
        
            time2.append( np.floor(self.time[count]) ) 
            signal2.append( signal_temp )
            error2.append( error_temp )
        
            count += nn
        
        self.time   = np.asarray(time2)
        self.signal = np.asarray(signal2)
        self.error  = np.asarray(error2)
    
    
    ## bin input signal_temp to just one data point
    def _bin(self, signal_temp, error_temp):
        if (len(signal_temp) > 1):
            signal_temp, error_temp = self.__mag2flux(signal_temp, error_temp)
            signal_temp, error_temp = self.__weighted_mean(signal_temp, error_temp)
            signal_temp, error_temp = self.__flux2mag(signal_temp, error_temp)
        
        return signal_temp, error_temp
        

    def __mag2flux(self, signal, error):
        flux = 10.**(-1.*signal/2.5)
        return 10.**(-1.*signal/2.5), np.abs( -flux*error*np.log(10.)/2.5 )
    
    
    def __flux2mag(self, signal, error):
        return -2.5*np.log10(signal), np.abs( -2.5* error/signal/np.log(10.))
        
    
    def __weighted_mean(self, signal, error):
        signal_mean = np.sum(signal/error**2.) / np.sum(1./error**2.) 
        error_mean  = np.sqrt( np.sum(error**2.) ) / np.sqrt( np.float(len(signal)) )
        return signal_mean, error_mean
    
    ### *********************************** ###
    ###  END of helper func for preprocess  ###
    ### *********************************** ###
    
    

    
    
    ##### ------------------------------- #####
    ##### --- END of quasar_drw class --- #####
    ##### ------------------------------- #####




