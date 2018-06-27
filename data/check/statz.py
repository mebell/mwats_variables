import json
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare
import os
from scipy import stats
import math
import datetime 

mod = []
avg_s = []
chi = []
slopes = []
avgm = []
fit_errors = []

def days_to_hmsm(days):
    hours = days * 24.
    hours, hour = math.modf(hours)
    mins = hours * 60.
    mins, min = math.modf(mins)
    secs = mins * 60.
    secs, sec = math.modf(secs)
    micro = round(secs * 1.e6)
    return int(hour), int(min), int(sec), int(micro)

def jd_to_date(jd):
    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I        
    C = B + 1524    
    D = math.trunc((C - 122.1) / 365.25)    
    E = math.trunc(365.25 * D)    
    G = math.trunc((C - E) / 30.6001)    
    day = C - E + F - math.trunc(30.6001 * G)    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715        
    return year, month, day

def jd_to_datetime(jd):
    year, month, day = jd_to_date(jd)
    frac_days,day = math.modf(day)
    day = int(day)
    hour,min,sec,micro = days_to_hmsm(frac_days)
    return datetime.datetime(year,month,day,hour,min,sec,micro)

def chi2(flux, error):
    w_mean = np.average(flux)
    add_this = []
    for s in range(len(flux)):
        d = (flux[s] - w_mean)**2 / (error[s])**2
        add_this.append(d)
    chi = sum(add_this)
    return chi / (len(fluxes) -1) 

def re_plot(s_time, s_flux, s_error, name, m, x):
        #### Plot linear ######
        plt.figure(figsize=(3.5, 4))
	ax = plt.subplot(211)
	plt.errorbar(range(0, len(s_flux)), s_flux, s_error, marker='.', color='k', linestyle='None', capsize=0)
	for f in range(len(s_time)-1):
	    if (s_time[f+1] - s_time[f]) > 10:
		plt.plot([f,f],[0, max(s_flux)*1.2], marker=None, color='0.6', linestyle=':', linewidth=2.0)
	plt.xlabel('Epoch Number')
	plt.ylabel('Flux Density (Jy)')
        #plt.text(5, 1, "$M=$"+str(m)+'% $\chi^{2}_{r}$='+str(x))
        plt.text(5, 1, "$M=$"+str(m)+'%')
	plt.ylim(0, max(s_flux)*1.2)
        locator = plt.MaxNLocator(nbins=5)
        ax.yaxis.set_major_locator(locator)
        locator = plt.MaxNLocator(nbins=5)
        ax.xaxis.set_major_locator(locator)
	plt.title(name.replace('_',' '))
        ######## Plot in time  ########
	ax = plt.subplot(212)
        dtime = []
        for t in s_time:
            test = jd_to_datetime(t+50000)
            #year = float(test.strftime('%Y')) + float(test.strftime('%j')) / 366
            dtime.append(test)
	plt.errorbar(dtime, s_flux, s_error, marker='.', color='k', linestyle='None', capsize=0)
	plt.xlabel('Date')
	plt.ylabel('Flux Density (Jy)')
        #plt.gcf().autofmt_xdate()
	plt.ylim(0, max(s_flux)*1.2)
        ##### Plot the average #####
        div = [0]
        for f in range(len(s_time)-1):
            if (s_time[f+1] - s_time[f]) > 10:
               div.append(f+1)
        div.append(len(s_flux))
        avg_s  = []
        avg_jd = []
        avg_jd_dt = []
        for p in range(len(div)-1):
            avs = np.mean(s_flux[div[p]:div[p+1]])
            avjd = np.mean(s_time[div[p]:div[p+1]])
            avjd_dt = jd_to_datetime(avjd+50000)
            avg_s.append(avs)
            avg_jd_dt.append(avjd_dt)
            avg_jd.append(avjd)
        #plt.plot(avg_jd_dt, avg_s,'k.-')
	plt.subplots_adjust(hspace=0.5)
        avg_m = (np.std(avg_s) / np.mean(avg_s)) * 100
        #### Do the fit to the average data
        slope, intercept, r_value, p_value, std_err = stats.linregress(avg_jd, avg_s)
        ys = []
        for t in range(len(s_time)):
            y = slope * s_time[t] + intercept # remove y minus offset
            ys.append(y)
        plt.plot(dtime, ys, 'r-')
        sigma = slope / std_err
        plt.text((min(dtime)), 0.8, "$\overline{M}=$"+str(round(avg_m,1))+'% $\\nabla_{S}=$'+str(round(sigma,1)))
        plt.xticks(rotation=30)
        locator = plt.MaxNLocator(nbins=7)
        ax.xaxis.set_major_locator(locator)
        locator = plt.MaxNLocator(nbins=5)
        ax.yaxis.set_major_locator(locator)
        xmin, xmax = plt.gca().get_xlim()
        xrange = xmax - xmin
        xmargin = xrange * 0.05
        xmin = xmin - xmargin
        xmax = xmax + xmargin
        plt.gca().set_xlim(xmin, xmax)        
        plt.tight_layout()
        plt.savefig(name+'_avg.png')
        plt.close()
        ###############################
        # Print all the details to HTML
        print '<tr>'
        print '   <td class="lalign">'+name+'</td><td>'+str(round(np.abs(sigma),1))+'</td><td>'+str(round(np.abs(slope*10000),3))+'</td><td>'+str(round(avg_m,1))+'</td>'
        print '   <td><img src="data/'+name+'_avg.png" alt="" border=1 height=300 width=300></img></td>'
        print '</tr>'
        ###############################
        return avg_m, slope, std_err

for file in glob.glob("GLEAM*.txt"):
    source = json.load(open(file))
    fluxes = source['peak_flux']
    errors = source['s_error']
    s_time = source['jd_time']
    #### Mod index and Chi Squared ####
    avg = np.mean(fluxes)
    avg_s.append(avg)
    std = np.std(fluxes)
    m = (std/avg)*100
    mod.append((std/avg)*100) 
    x = chisquare(fluxes, ddof=len(fluxes))[0]
    chi.append(x)
    ###### Re_plot #######
    name=file.split('.txt')[0]
    am, slope, fit_error = re_plot(s_time, fluxes, errors, name, round(m, 1), round(x, 1))
    avgm.append(am)
    slopes.append(slope)
    fit_errors.append(fit_error)
    sig = slope / fit_error

