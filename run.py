import numpy as np
import pandas as pd
import os
import sys
from quasar_drw import quasar_drw as qso_drw
from plot import plot
import query as query

from scipy import stats
from scipy.optimize import curve_fit

def read_quasar_catalog(filename):

    df = pd.read_csv(filename)
    return df

def read_DES(name,band):

    file_DES = "data/DES/"+name+"/"+band+".csv"
    data_DES =  pd.read_csv(file_DES,comment="#")
    data_DES_return = data_DES[["MJD_OBS","FLUX_PSF","FLUX_ERR_PSF"]].values
    return data_DES_return[:,0],22.5-2.5*np.log10(data_DES_return[:,1]),data_DES_return[:,2]*1.09/data_DES_return[:,1]
#    return data_DES_return[:,0],data_DES_return[:,1],data_DES_return[:,2]

def read_DES_Y5(name,band):

    file_DES = "data/DES_Y5/"+name+"/"+band+".csv"
    data_DES =  pd.read_csv(file_DES,comment="#")
    data_DES_return = data_DES[["MJD_OBS","FLUX_PSF","FLUX_ERR_PSF"]].values
    return data_DES_return[:,0],22.5-2.5*np.log10(data_DES_return[:,1]),data_DES_return[:,2]*1.09/data_DES_return[:,1]
#    return data_DES_return[:,0],data_DES_return[:,1],data_DES_return[:,2]

def check_SDSS(name,band):

    file_SDSS = "data/SDSS/"+name+".csv"
    if os.path.exists(file_SDSS): return True
    else: return False

def check_DES(name,band):

    file_DES = "data/DES/"+name+"/"+band+".csv"
    if os.path.exists(file_DES): return True
    else: return False

def check_DES_Y5(name,band):

    file_DES = "data/DES_Y5/"+name+"/"+band+".csv"
    if os.path.exists(file_DES): return True
    else: return False

def read_SDSS(name,band):

    file_SDSS = "data/SDSS/"+name+".csv"
    data_SDSS =  pd.read_csv(file_SDSS,comment="#")
    data_SDSS_return = data_SDSS[["mjd_"+band,"psfmag_"+band,"psfMagErr_"+band]].values 
    return data_SDSS_return[:,0],data_SDSS_return[:,1],data_SDSS_return[:,2]
#    return data_SDSS_return[:,0],10**((data_SDSS_return[:,1]-22.5)/(-2.5)),data_SDSS_return[:,2]*10**((data_SDSS_return[:,1]-22.5)/(-2.5))/1.09

def create_dir(directory):

    if not os.path.exists(directory):
        os.makedirs(directory)

def fitting(time,signal,error,period):

    def sin_func(x,amplitude,ref_day,median):
        return amplitude*np.sin(2*np.pi*x/period+ref_day)+median
    p0 = [1,50000,20]
    popt, pcov = curve_fit(sin_func, time, signal,p0=p0)
    
    xn = np.linspace(np.min(time)-100,np.max(time)+100,10000)
    yn = sin_func(xn,*popt)
    return xn,yn
  
def check_period_max_amp(_freq, psd):
    # only search for peak amp between 1.5 year and ? year.
    try:
        period = _freq[_freq<5*365][_freq>1.5*365]
        amp = psd[_freq<5*365][_freq>1.5*365]
        period_max = period[np.where(amp==np.max(amp))]

        return period_max

    except: return None
 
def save_freq_amp(_freq, psd,filename):

     save_data = zip(_freq,psd)
     np.savetxt(filename,save_data,delimiter=",",header="period,amplitude")

def save_freq_confidence_level(_freq,boundary_all,filename):

     period = np.array([_freq]).T
     psd_array = np.array([boundary_all]).T
     print period,psd_array
     save_data = np.concatenate((period, psd_array), axis=1)
     np.savetxt(filename,save_data,delimiter=",",header="period,confidence_level")

def clean_parameters_list(parameters_list):

    num_parameters = len(parameters_list[0])
    for i in range(0,num_parameters-1):
        array = parameters_list[:,i]
        sigma = np.std(array)
        mean = np.mean(array)
        parameters_list = parameters_list[(array<mean+2*sigma) & (array>mean-2*sigma)] 

    return parameters_list


def tailored_simulation(lc,time,signal,band,z,name,output_dir,periodogram,lightcurve,random_state):

    psd_mock_all = []
    parameters_list =  lc.fit_drw_emcee(nwalkers=200, burnin=200, Nstep=400,random_state=random_state)
#    parameters_list =  lc.fit_drw_emcee(nwalkers=10, burnin=10, Nstep=20,random_state=random_state)
    parameters_list_good = clean_parameters_list(parameters_list)
    for i in range(10000):
#    for parameters in parameters_list:
        tau,c,b = np.exp(parameters_list_good[random_state.randint(len(parameters_list_good))])
#        tau,b,c = np.exp(parameters)
        mock_time,mock_signal = lc.generate_mock_lightcurve(tau,b,c,time,z,random_state=random_state)
        #print np.mean(signal)
        mock_signal_correct = mock_signal#+np.mean(signal)
        #lightcurve.plot_mock_curve(mock_time,mock_signal_correct,band)
        #print mock_signal_correct
        _freq_mock, psd_mock = lc.periodogram(mock_time,mock_signal)
        psd_mock_all.append(psd_mock)
        periodogram.plot_mock_periodogram(_freq_mock, psd_mock,band)

    psd_at_each__freq = zip(*psd_mock_all)
    percentiles = [68.27,95.45,99.0,99.74]
    boundary_all = []
    confidence_levels = []
    _freq, psd_true = lc.ls_astroML()
    for percentile in percentiles:
        bounday_psd_at_each__freq = [np.percentile(psd,50+percentile/2.) for psd in psd_at_each__freq]
        boundary_all.append(bounday_psd_at_each__freq)
    for i in range(len(_freq)):
        confidence_level_at_each__freq = float(stats.percentileofscore(psd_at_each__freq[i],psd_true[i]))
        confidence_levels.append(100.-confidence_level_at_each__freq)
    save_freq_confidence_level(_freq,confidence_levels,output_dir+name+"/confidence_"+band+".csv")
    periodogram.plot_confidence_level(_freq_mock, psd_mock_all, band)



#####################
##  main functions ##
#####################

def main(ra,dec,name,z):

 
    output_dir = "output/"
    random_state = np.random.RandomState(0)
#    query_sdss = query.Stripe82()
#    query_sdss.q(ra,dec,name,dist=5.0)
    bands = ["g","r","i","z"]

    periodogram = plot(2,2)
    lightcurve = plot(4,1,figsize=(8,8),sharex=True)

    create_dir(output_dir+name)
    for band in bands:

        if check_DES(name,band):
            mjd_DES,mag_DES,magerr_DES = read_DES(name,band)
            lc = qso_drw(mjd_DES,mag_DES,magerr_DES, z, preprocess=True)
            if check_SDSS(name,band):
                mjd_SDSS,mag_SDSS,magerr_SDSS = read_SDSS(name,band) 
                lc2 = qso_drw(mjd_SDSS,mag_SDSS,magerr_SDSS, z, preprocess=True) 
                time, signal, error = lc2.get_lc()
                lc.add_lc(time, signal, error, preprocess=False)
            else: print "SDSS not found !"
        else:
            print  "DES not found !"
            if  check_SDSS(name,band):
                mjd_SDSS,mag_SDSS,magerr_SDSS = read_SDSS(name,band)
                lc = qso_drw(mjd_SDSS,mag_SDSS,magerr_SDSS, z, preprocess=True)
            else: print "SDSS not found !"
        if check_DES_Y5(name,band):
            mjd_DES_Y5,mag_DES_Y5,magerr_DES_Y5 = read_DES_Y5(name,band)
            lc3 = qso_drw(mjd_DES_Y5,mag_DES_Y5,magerr_DES_Y5, z, preprocess=True)
            time, signal, error = lc3.get_lc()
            lc.add_lc(time, signal, error, preprocess=False)
        else: "DES_Y5 not found !"
        
        time, signal, error = lc.get_lc()
        lightcurve.plot_light_curve(time,signal,error,band)

        ## periodogram
        if len(time) > 3:
            _freq, psd = lc.ls_astroML()
            periodogram.plot_periodogram(_freq, psd,band)
  
            _freq_max = check_period_max_amp(_freq, psd)
            if _freq_max is not None:
                xn,yn = fitting(time,signal,error,_freq_max)
                lightcurve.plot_fit_curve(xn,yn,band)
            tailored_simulation(lc,time,signal,band,z,name,output_dir,periodogram,lightcurve,random_state)
            save_freq_amp(_freq, psd,output_dir+name+"/periodogram_"+band+".csv")
    lightcurve.savefig(output_dir+name,"/lightcurve.png",name)
    periodogram.savefig(output_dir+name,"/periodogram.png",name)
    
def is_candidate_good(name):

    upper_period = 7*365
    lower_period = 1.5*365
    amp_cut = 0.8
    bands = ["g","r","i","z"]
    search_dir = "test/"
    pass_times = 0
    for band in bands:
        if os.path.exists(search_dir+name+"/periodogram_"+band+".csv"):
            periodogram = pd.read_csv(search_dir+name+"/periodogram_"+band+".csv",names=["period","amplitude"],skiprows=1)
            periodogram1 = periodogram[periodogram["period"]<upper_period]
            periodogram2= periodogram1[periodogram1["period"]>lower_period]
            N = len(np.where(periodogram2["amplitude"]>amp_cut)[0])
            # check baseline is long enough
            if N>0:
                pass_times += 1
        else: print "Can not read csv."
    if pass_times == 4:
        print name
        return True
    return False

def is_candidate_strong(name):

    upper_period = 8*365
    lower_period = 500
    bands = ["g","r","i","z"]
    search_dir = "output/"
    for band in bands:
        if os.path.exists(search_dir+name+"/confidence_"+band+".csv"):
            periodogram = pd.read_csv(search_dir+name+"/periodogram_"+band+".csv",names=["period","confidence"],skiprows=1)
            periodogram1 = periodogram[periodogram["period"]<upper_period]
            periodogram2= periodogram1[periodogram1["period"]>lower_period]
            peak_confidence = max(periodogram2["confidence"])
            N = len(np.where(periodogram2["amplitude"]>amp_cut)[0])
            # check baseline is long enough
            if N>0:
                pass_times += 1
        else: print "Can not read csv."
    if pass_times == 4:
        print name
        return True
    return False
        
def find_good_cadidates():

    """ Making simple power cut in periodogram """
    can = open("candidates.txt","a")
    print "Finding strong candidates :"
    df_quasar_list = read_quasar_catalog("strip82_catalog.csv")
    for index, row in df_quasar_list.iterrows():
        name = row["name"]
        if is_candidate_good(name):
            can.write(name+"\n")
def record_confidence_peak():

    upper_period = 8*365
    lower_period = 500
    bands = ["g","r","i","z"]
    search_dir = "output/"
    f_confidence = open("statistics/candidates_confidence.txt","w")
    f_confidence.write("name,peak_g,peak_r,peak_i,peak_z\n")
    df_quasar_list = read_quasar_catalog("strip82_catalog.csv")
    for index, row in df_quasar_list.iterrows():
        name = row["name"]
        f_confidence.write(name)
        for band in bands:
            if os.path.exists(search_dir+name+"/confidence_"+band+".csv"):
                periodogram = pd.read_csv(search_dir+name+"/confidence_"+band+".csv",names=["period","confidence"],skiprows=1)
                upper_period =max(periodogram["period"])/3
                periodogram1 = periodogram[periodogram["period"]<upper_period]
                periodogram2= periodogram1[periodogram1["period"]>lower_period]
                if len(periodogram2) != 0 :
                    peak_confidence = min(periodogram2["confidence"])
                    f_confidence.write(","+str(peak_confidence))
                else: f_confidence.write(",999")
            else: f_confidence.write(",999")
        f_confidence.write("\n")
    f_confidence.close()

def show_number_significance():

    data_significance = pd.read_csv("statistics/candidates_confidence.txt")
    bands = ["g","r","i","z"]
    pvalues = np.logspace(np.log10(100),np.log10(0.0001),1000)
    number_sig_plot = plot(1,1)
    for band in bands:
        N_list = []
        for pvalue in pvalues:
            N = len(np.where(data_significance["peak_"+band]<pvalue)[0])
            N_list.append(N)
        number_sig_plot.plot(pvalues,N_list,band,log=True)
    number_sig_plot.savefig("statistics/","number_significance.png","Number vs. P-value")



def find_strong_candidates():

    """ Compare the power with mock light curves """
#    cand_strong = open("statistics/strong_candidates.txt","w")
    cand = pd.read_csv("statistics/candidates_confidence.txt")
    bands = ["g","r","i","z"]
    for band in bands:
        cand = cand[cand["peak_"+band]<2.25]
    cand.to_csv("statistics/strong_candidates.csv")

def record_success(name,success):

    """ record the susseceful run """

    f = open("finished.list","a")
    f.write(str(name)+","+str(success)+"\n")
    f.close()

 
if __name__ == "__main__":

#    df_quasar_list = read_quasar_catalog("strip82_catalog.csv")
#    for index, row in df_quasar_list.iterrows():
#        ra,dec,name = row["ra"],row["dec"],row["name"]
#        main(ra,dec,name)
        
#    find_strong_cadidates()
#    main(43.0611,-0.470451,"J025214.66-002813.62",1.3)
#    record_confidence_peak()
#    show_number_significance()
#    find_strong_candidates()
    print sys.argv[1].split(",")
    name,ra,dec,z = sys.argv[1].split(",")
    ra = float(ra)
    dec = float(dec)
    z = float(z)
    try:
        main(ra,dec,name,z)
        record_success(name,True)
    except :
        record_success(name,False)
