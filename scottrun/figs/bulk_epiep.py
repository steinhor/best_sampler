import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
#sformatter=ScalarFormatter(useOffset=True,useMathText=False)
#sformatter.set_scientific(True)
#sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,4.5))
fig = plt.figure(1)
ax = fig.add_axes([0.14,0.14,0.84,0.84])

mydata = np.loadtxt('bulk_rho0.dat',skiprows=1,unpack=True)
PIdata=-mydata[0]
PIoverPtarget=-mydata[1]
PIcheck=mydata[2]
echeck=mydata[3]
Tbulk=mydata[4]
pbulkscale=1.0/(1.0+mydata[5]/3.0)
Epion=mydata[6]-0.138
Enucleon=mydata[7]-0.938

mydata2 = np.loadtxt('bulk_rho0.1.dat',skiprows=1,unpack=True)
PIdata2=-mydata2[0]
PIoverPtarget2=-mydata2[1]
PIcheck2=mydata2[2]
echeck2=mydata2[3]
Tbulk2=mydata2[4]
pbulkscale2=1.0/(1.0+mydata2[5]/3.0)
Epion2=mydata2[6]-0.138
Enucleon2=mydata2[7]-0.938

dPde=0.152181
Pressure=0.0404889
epsilon=0.24014
#PIcheck=PIcheck-dPde*(echeck-1.0)*epsilon/(Pressure*PIoverPtarget);

plt.plot(100.0*PIoverPtarget,1000*Epion,marker='o',markersize=8,linestyle='--',linewidth=1,color='g')
plt.plot(100.0*PIoverPtarget,1000*Enucleon,marker='o',markersize=8,linestyle='--',linewidth=1,color='g')
plt.plot(100.0*PIoverPtarget2,1000*Epion2,marker='s',markersize=8,linestyle='--',linewidth=1,color='r')
plt.plot(100.0*PIoverPtarget2,1000*Enucleon2,marker='s',markersize=8,linestyle='--',linewidth=1,color='r')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticklabels(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.25)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(0,500,50), minor=False)
ax.set_yticklabels(np.arange(0,500,50), minor=False)
ax.set_yticks(np.arange(0,500,10), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.ylim(200,400)

plt.xlabel('$-\Pi_{\\rm hydro}/P$ (percent)', fontsize=18, weight='normal')
plt.ylabel('$\langle E\\rangle-M$  (MeV)',fontsize=18)

plt.text(1.2,260,"protons",color='k',fontsize='22',ha='right')
plt.text(1.2,362,"pions",color='k',fontsize='22',ha='right')


plt.savefig('bulk_epiep.pdf',format='pdf')
os.system('open -a Preview bulk_epiep.pdf')
#plt.show()
quit()
