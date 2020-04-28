import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

mydata = np.loadtxt('rhospectra.dat',skiprows=0,unpack=True)
x=mydata[0]
y=mydata[1]
mydatap = np.loadtxt('rhospectra_polemass.dat',skiprows=0,unpack=True)
xp=mydatap[0]
yp=mydatap[1]
yscale('log')
plt.plot(x,y,linestyle='-',linewidth=3,color='r')
plt.plot(xp,yp,linestyle='-',linewidth=3,color='b')

#sangwookdata = np.loadtxt('sangwook/pTspec_Msample21_2214.dat',skiprows=1,unpack=True)
sangwookdata = np.loadtxt('sangwook/pTspec_Msample21_113.dat',skiprows=1,unpack=True)
sx=sangwookdata[0]
sy=sangwookdata[1]/125000.0
sangwookdata = np.loadtxt('sangwook/pTspec_Msample00_113.dat',skiprows=1,unpack=True)
sxp=sangwookdata[0]
syp=sangwookdata[1]/125000.0
yscale('log')
plt.plot(sx,sy,linestyle='--',linewidth=2,color='r')
plt.plot(sxp,syp,linestyle='--',linewidth=2,color='b')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,2.1,0.2), minor=False)
ax.set_xticklabels(np.arange(0,2.1,0.1), minor=False, family='serif')
ax.set_xticks(np.arange(0,2.1,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.5)

#ax.set_yticks(np.arange(0.00001), minor=False)
#ax.set_yticklabels(np.arange(0,200000,40000), minor=False, family='serif')
#ax.set_yticks(np.arange(0,200000,10000), minor=True)
plt.ylim(0.000001,0.01)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$p$ (GeV/$c$)', fontsize=18, weight='normal')
plt.ylabel('$(1/V)dN/d^3p$ [GeV$/c$]$^{-3}$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('spectra_rho.pdf',format='pdf')
os.system('open -a Preview spectra_rho.pdf')
#plt.show()
quit()
