#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys

n_bins=75
if len(sys.argv)<1: exit();

if sys.argv[len(sys.argv)-1]=='-c':
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
else:
    fig, axs = plt.subplots(1, len(sys.argv)-1, sharey=True, tight_layout=True)
K=0
enl=[]
mV,MV=-400,400
xmin,xmax=mV-5,MV+5
xa=np.linspace(xmin,xmax,1000)
plt.ylabel("Fraction of molecules")
#plt.title("Energy distribution.")
for f in sys.argv[1:]:
    if f=='-c': continue
    if not enl: enl=list(np.loadtxt(f,delimiter=' ',usecols=(1,)))
    elif sys.argv[len(sys.argv)-1]=='-c': enl+=list(np.loadtxt(f,delimiter=' ',usecols=(1,)))
    else: enl=list(np.loadtxt(f,delimiter=' ',usecols=(1,)))
    enl=[x for x in enl if x<400]
    if sys.argv[len(sys.argv)-1]=='-c': pass
    else:
        mean=sum(enl)/len(enl); stdev=np.std(enl);
        ya=norm.pdf(xa,mean,stdev) #*len(enl)
        if len(sys.argv)<3:
            myplt=axs.hist(enl,bins=n_bins,range=(mV-5,MV+5),density=True,alpha=0.5,color='y')
            axs.plot(xa,ya,'b',linewidth=1,linestyle='-')
            axs.set_xlabel("Energy (kJ/mol)")
            axs.set_title("Energy distribution. Total molecules: "+str(len(enl)))
            axs.axvline(x=mean,linestyle=':',color='k',label="mean")
            axs.axvline(x=mean-3*stdev,linestyle=':',color='r',label="cutoff")
        else:
            myplt=axs[K].hist(enl,bins=n_bins,range=(mV-5,MV+5),density=True,alpha=0.5,color='y')
            axs[K].plot(xa,ya,'b',linewidth=1,linestyle='-')
            axs[K].set_xlabel("Energy (kJ/mol)")
            axs[K].set_title("Energy distribution. Total molecules: "+str(len(enl)))
            axs[K].axvline(x=mean,linestyle=':',color='k',label="mean")
            axs[K].axvline(x=mean-3*stdev,linestyle=':',color='r',label="cutoff")
        print(f,"\tAverage energy: ",mean,"\tStd Dev: ",stdev)
        K+=1

if sys.argv[len(sys.argv)-1]=='-c':
        axs.hist(enl,bins=n_bins,range=(mV-5,MV+5),density=True,alpha=0.5,color='y')
        mean=sum(enl)/len(enl); stdev=np.std(enl);
        ya=norm.pdf(xa,mean,stdev) #*len(enl)
        axs.plot(xa,ya,'b',linewidth=1,linestyle='-')
        axs.set_xlabel("Energy (kJ/mol)")
        axs.set_title("Energy distribution. Total molecules: "+str(len(enl)))
        axs.axvline(x=mean,linestyle=':',color='k',label="mean")
        axs.axvline(x=mean-3*stdev,linestyle=':',color='r',label="cutoff")
        print("Totalling: Average energy: ",mean,"\tStd Dev: ",stdev)
plt.legend()
plt.show()
