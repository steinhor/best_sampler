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

mydata = np.loadtxt('rhomassdist.dat',skiprows=0,unpack=True)
x=mydata[0]
y=mydata[1]
mydatap = np.loadtxt('rhomassdist_polemass.dat',skiprows=0,unpack=True)
xp=mydatap[0]
yp=mydatap[1]
yscale('log')
plt.plot(x,y,linestyle='-',linewidth=2,color='r')
plt.plot(xp,yp,linestyle='-',linewidth=2,color='b')

mydata_d = np.loadtxt('deltamassdist.dat',skiprows=0,unpack=True)
x=mydata_d[0]
y=mydata_d[1]
mydatap_d = np.loadtxt('deltamassdist_polemass.dat',skiprows=0,unpack=True)
xp=mydatap_d[0]
yp=mydatap_d[1]
yscale('log')
plt.plot(x,y,linestyle='-',linewidth=2,color='r')
plt.plot(xp,yp,linestyle='-',linewidth=2,color='b')

#sangwookdata = np.loadtxt('sangwook/pTspec_Msample21_2214.dat',skiprows=1,unpack=True)
#sangwookdata = np.loadtxt('sangwook/pTspec_Msample21_113.dat',skiprows=1,unpack=True)
#sx=sangwookdata[0]
#sy=sangwookdata[1]
#sangwookdatap = np.loadtxt('sangwook/pTspec_Msample00_113.dat',skiprows=1,unpack=True)
#sxp=sangwookdatap[0]
#syp=sangwookdatap[1]
#yscale('log')
#plt.plot(sx,sy,linestyle='--',linewidth=2,color='r')
#plt.plot(sx,sy,linestyle='--',linewidth=2,color='b')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,2.1,0.25), minor=False)
ax.set_xticklabels(np.arange(0,2.1,0.25), minor=False, family='serif')
ax.set_xticks(np.arange(0,2.1,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.25,2)

#ax.set_yticks(np.arange(0.00001), minor=False)
#ax.set_yticklabels(np.arange(0,200000,40000), minor=False, family='serif')
#ax.set_yticks(np.arange(0,200000,10000), minor=True)
plt.ylim(0.001,10)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$M$ (GeV)', fontsize=18, weight='normal')
plt.ylabel('$dN/dM$ (normalized to unity)',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('massdist.pdf',format='pdf')
os.system('open -a Preview massdist.pdf')
#plt.show()
quit()
