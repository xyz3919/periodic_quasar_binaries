import numpy as np
import pandas as pd

from quasar_drw import quasar_drw as qso_drw
from plot import plot
import query as query

from scipy.optimize import leastsq

def read_DES(name,band):

    file_DES = "data/DES/"+name+"/"+band+".csv"
    data_DES =  pd.read_csv(file_DES,comment="#")
    data_DES_return = data_DES[["MJD_OBS","FLUX_PSF","FLUX_ERR_PSF"]].values
    return data_DES_return[:,0],22.5-2.5*np.log10(data_DES_return[:,1]),data_DES_return[:,2]*1.09/data_DES_return[:,1]

def read_DES_Y5(name,band):

    file_DES = "data/DES_Y5/"+name+"/"+band+".csv"
    data_DES =  pd.read_csv(file_DES,comment="#")
    data_DES_return = data_DES[["MJD_OBS","FLUX_PSF","FLUX_ERR_PSF"]].values
    return data_DES_return[:,0],22.5-2.5*np.log10(data_DES_return[:,1]),data_DES_return[:,2]*1.09/data_DES_return[:,1]

def read_SDSS(name,band):

    file_SDSS = "data/SDSS/"+name+".csv"
    data_SDSS =  pd.read_csv(file_SDSS,comment="#")
    data_SDSS_return = data_SDSS[["mjd_"+band,"psfmag_"+band,"psfMagErr_"+band]].values 
    return data_SDSS_return[:,0],data_SDSS_return[:,1],data_SDSS_return[:,2]

def sin_func(x,_freq_max,amplitude,ref_day,median):

    return amplitude*np.sin(2*np.pi*x/_freq_max+ref_day)+median

def residual(p, x, y,_freq_max):
    return y - sin_func(x, _freq_max, *p)

def main(ra,dec,name):

    query_sdss = query.Stripe82()
    query_sdss.q(ra,dec,name,dist=5.0)
    bands = ["g","r","i","z"]
    z = 1.3

    periodogram = plot(2,2)
    lightcurve = plot(4,1,figsize=(8,8),sharex=True)

    for band in bands:

        mjd_DES,mag_DES,magerr_DES = read_DES(name,band)
        lc = qso_drw(mjd_DES,mag_DES,magerr_DES, z, preprocess=True)

        mjd_SDSS,mag_SDSS,magerr_SDSS = read_SDSS(name,band) 
        lc2 = qso_drw(mjd_SDSS,mag_SDSS,magerr_SDSS, z, preprocess=True) 
        time, signal, error = lc2.get_lc()
        lc.add_lc(time, signal, error, preprocess=False)

        mjd_DES_Y5,mag_DES_Y5,magerr_DES_Y5 = read_DES_Y5(name,band)
        lc3 = qso_drw(mjd_DES_Y5,mag_DES_Y5,magerr_DES_Y5, z, preprocess=True)
        time, signal, error = lc3.get_lc()
        lc.add_lc(time, signal, error, preprocess=False)
        
        time, signal, error = lc.get_lc()
        
        lightcurve.plot_light_curve(time,signal,error,band)
        ## periodogram
        _freq, psd = lc.ls_astroML()
        periodogram.plot_periodogram(_freq, psd,band)
 
        _freq_max = _freq[np.where(psd==np.max(psd))]
        p0 = [0.5,20,51000]
        popt,pcov = leastsq(residual, p0,\
                                  args=(time,signal,np.array([_freq_max]*len(time))))
        yn = sin_func(time,[_freq_max]*len(time),*popt)
        plot_fit_curve(time,yn,band)

    lightcurve.savefig("output/","lightcurve.png",name)
    periodogram.savefig("output/","periodogram.png",name)

if __name__ == "__main__":

#    main(41.1321,-0.428554,"J024431.71-002542.79")
    main(43.0611,-0.470451,"J025214.66-002813.62")
