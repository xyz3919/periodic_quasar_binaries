import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

class plot:

    """
    All the useful plot functions.
    
    by Yu-Ching (Tony) Chen
    ycchen@illinois.edu
    """

    def __init__(self, cols, rows,**args):

        self.f,self.axes  = plt.subplots(cols, rows,**args)

    def plot(self,x,y,band,log=False):

        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        self.axes.plot(x,y,label=band,c=color_list[band])
        if log: 
            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
        self.axes.legend()



    def plot_periodogram(self,_freq, psd,band):

        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.plot(_freq/365,psd,label=band,c=color_list[band])
        ax.set_xlim(0.8,10)
        ax.set_ylim(0,1)
        ax.set_xscale("log")
        ax.set_xticks([1,2,4,8,10])
        ax.set_xticklabels([1,2,4,8,10])
        ax.fill_betweenx([0.0, 1.05], 0.8,  500./365., color='grey', alpha='0.5')
        ax.fill_betweenx([0.0, 1.05], max(_freq)/365/3,  max(_freq)/365, color='grey', alpha='0.5')
        if band == "i" or  band == "z":
            ax.set_xlabel("Period(yr)")
        if band == "g" or  band == "i":
            ax.set_ylabel("Power")
        ax.annotate(band, xy=(-12, -12), xycoords='axes points',
                    size=12, ha='right', va='top', color=color_list[band],
                    bbox=dict(boxstyle='round', fc='w'))
#        ax.legend()


    def plot_mock_periodogram(self,_freq, psd,band):

        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.plot(_freq/365,psd,label=band,c="grey",linewidth=0.01)

    def plot_confidence_level(self,_freq, psd_total,band):
        ax_list = {"g":self.axes[0,0],"r":self.axes[0,1],\
                   "i":self.axes[1,0],"z":self.axes[1,1]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        psd_at_each__freq = zip(*psd_total)         
        percentiles = [68.27,95.45,99.0,99.74,99.99]
        for percentile in percentiles:
           bounday_psd_at_each__freq = [np.percentile(psd,50+percentile/2.) for psd in psd_at_each__freq]
           ax.plot(_freq/365,bounday_psd_at_each__freq,"--",c="black",linewidth=0.3)

    def plot_light_curve(self,time,signal,error,band):

        ax_list = {"g":self.axes[0],"r":self.axes[1],\
                   "i":self.axes[2],"z":self.axes[3]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.errorbar(time,signal,yerr=error,fmt='o',\
                    label=band,c=color_list[band])
        ax.set_ylim(np.max(signal)+0.1,np.min(signal)-0.1)
        if band == "z":
            ax.set_xlabel("MJD")
        ax.set_ylabel("Magnitude") 
        ax.annotate(band, xy=(-12, -12), xycoords='axes points',
                    size=12, ha='right', va='top', color=color_list[band],
                    bbox=dict(boxstyle='round', fc='w'))

    def plot_fit_curve(self,time,signal,band):

        ax_list = {"g":self.axes[0],"r":self.axes[1],\
                   "i":self.axes[2],"z":self.axes[3]}
        color_list = {"g":"g","r":"orange",\
                      "i":"brown","z":"purple"}
        ax = ax_list[band]
        ax.plot(time,signal,\
                label=band,c=color_list[band])

    def plot_mock_curve(self,time,signal,band):

        ax_list = {"g":self.axes[0],"r":self.axes[1],\
                   "i":self.axes[2],"z":self.axes[3]}
        ax = ax_list[band]
        ax.plot(time,signal,label=band,c="grey",linewidth=0.001)

    def savefig(self,dir_output,name,title):

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        self.f.suptitle(title)
        self.f.savefig(dir_output+name)
        plt.close()
        

