#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import warnings
warnings.filterwarnings("ignore")
from obspy.core import read
import operator
import argparse
from sys import exit
import os.path
plt.style.use('bmh')

#parser = argparse.ArgumentParser(description=" Calculation of depth and Vp/Vs ratio following Zhu & Kanamori (2000)")
#parser.add_argument("-f1",help="f1 frequency for the band-pass filter",dest="f1",default=0.5, type=float)
#parser.add_argument("-f2",help="f2 frequency for the band-pass filter",dest="f2",default=1.5, type=float)
#parser.add_argument("infile",help="Input file",type=str)
#parser.add_argument("-fil",help="Apply filter if exists",dest="fil",action="store_true",default=False)
#parser.add_argument("-phf",help="Plots trace with final phases",dest="phf",action="store_true",default=False)
#param = parser.parse_args()

def t2sp(t1,dt,t0=0):
    return int(round(t1/dt))-int(t0/dt)

Hc = np.arange(13,25,1) # Conrad depth domain
Hm = np.arange(30,40,1) # Moho depth domain
K = np.arange(1.70,1.90,0.1) # Kappa domain
w1=0.50; w2=0.10; w3=0.30; w4=.10 # Weights for the stack calculation
Vps=[7.0]
Vps = np.arange(6.0,7.0,.2) # Vp domain

st = read("data/*pws*.sac")
def get_doubletimes(sacfile):
    dt = sacfile.stats.delta
    sta =sacfile.stats.station
    lat = sacfile.stats.sac.stla
    lon = sacfile.stats.sac.stlo
    t0 = sacfile.stats.sac.b
    t1 = sacfile.stats.sac.e
    t = np.arange(t0,t1+dt,dt)
    apil = []
    k_f = []
    hc_f = []
    hm_f = []
    tpmp = []
    tpb = []
    tsb = []
    for Vp in Vps:
        for k in K:
            for hc in Hc:
                for hm in Hm:
                    t_pmp = round((2. * hm) / Vp,2)
                    t_pb  = round((2. * hc) / Vp,2)
                    t_sms = round((2. * k * hm) / Vp,2)
                    t_sb  = round((2. * k * hc) / Vp,2)
                    a_pmp = sacfile.data[t2sp(t_pmp,dt,t0=t0)]  # Amp of PmP
                    a_pb  = sacfile.data[t2sp(t_pb,dt,t0=t0)]   # Amp of Pb
                    a_sms = sacfile.data[t2sp(t_sms,dt,t0=t0)]  # Amp of SmS
                    a_sb  = sacfile.data[t2sp(t_sb,dt,t0=t0)]   # Amp of Sb
                    apil.append(a_pmp * w1 + a_pb * w2 + a_sms * w3 + a_sb * w4)
                    k_f.append(k)
                    hc_f.append(hc)
                    hm_f.append(hm)
                    tpmp.append(t_pmp)
                    tpb.append(t_pb)
    
    idx = np.where(apil == min(apil))[0][0]
    hm = hm_f[idx]
    hc = hc_f[idx]
    k = k_f[idx]
    t_pmp = tpmp[idx]
    t_pb = tpb[idx]
    return sta,t_pb, t_pmp, hc, hm, k, lat, lon

def plot_times(sacfile,t_b,t_pmp,hc,hm,k):
    fig1 = plt.figure(1,figsize=(10,5))
    dt = sacfile.stats.delta
    sta =sacfile.stats.station
    t0 = sacfile.stats.sac.b
    t1 = sacfile.stats.sac.e
    t = np.arange(t0,t1+dt,dt)
    a_pmp = sacfile.data[t2sp(t_pmp,dt,t0=t0)]  # Amp of PmP
    a_pb =  sacfile.data[t2sp(t_pb,dt,t0=t0)]   # Amp of Pb
    title="%s | Conrad: %4.2f km   |  Moho: %4.2f km   | Vp/Vs= %3.2f" % (sta,hc,hm,k)
    plt.title(title)
    plt.plot(t, sacfile.data, 'k-')
    plt.plot(t_pmp,a_pmp, 'r*-',label="PmP", ms=10)
    plt.plot(t_pb,a_pb, 'b*-',label="Pb",ms=10)
    plt.xlim([t0,t1])
    plt.legend(loc=1)
    of = "fig_" + sta + ".png"
    plt.savefig(of,dpi=300)
    plt.close(fig1)
#    plt.show()

sta1 = ""
d = {}
tmpdict = {}
for tr in st:
    (sta,t_pb, t_pmp, hc, hm, k, lat, lon) = get_doubletimes(tr)
    tmpdict["tpb"] = t_pb
    tmpdict["tpmp"] = t_pmp
    tmpdict["hc"] = hc
    tmpdict["hm"] = hm
    tmpdict["k"] = k
    tmpdict["lat"] = lat
    tmpdict["lon"] = lon
    sta1=sta
    if not sta1 in d.keys():
        d[sta] = {"tpb":[t_pb], "tpmp":[t_pmp], "hc":[hc], "hm":[hm], "k":[k], "lat":[lat],"lon":[lon]}
    else:
         for k in d[sta].keys():
             d[sta][k].append(tmpdict[k])
for stat in d.keys():
    means = []
    stdvs = []
    for k,v in d[stat].items():
        means.append(np.mean(v))
        stdvs.append(np.std(v))
    if os.path.isfile("data/" + stat + ".PWS.sac"):
        sacfile = read("data/" + stat + ".PWS.sac")
        plot_times(sacfile[0],means[0],means[1],means[2],means[3],means[4])
        # Info: # means[0]- t_pmp, means[1] - t_pb, means[2] - Hc, means[3] -Hm, means[4]-K
    print("%s %7.4f %7.4f %4.2f %4.2f" % (str(stat),d[stat]['lon'][0], d[stat]['lat'][0], means[0], stdvs[0]))
#    print("%s %4.2f %4.2f" % (str(stat),means[0], means[1]))
