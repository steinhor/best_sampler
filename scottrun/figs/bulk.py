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
plt.figure(figsize=(6,7))
fig = plt.figure(1)
ax = fig.add_axes([0.16,0.12,0.8,0.4])
dPde=0.152181
Pressure=0.0404889
epsilon=0.24014

mydata = np.loadtxt('bulk.dat',skiprows=1,unpack=True)
PIdata=-mydata[0]
PIoverPtarget=-mydata[1]
PIcheck=mydata[2]
echeck=mydata[3]
Tbulk=mydata[4]
pbulkscale=1.0/(1.0+mydata[5]/3.0)
#PIcheck=PIcheck+dPde*(echeck-1.0)*epsilon/(Pressure*PIoverPtarget);

mydata2 = np.loadtxt('bulk_rho0.1.dat',skiprows=1,unpack=True)
PIdata2=-mydata2[0]
PIoverPtarget2=-mydata2[1]
PIcheck2=mydata2[2]
echeck2=mydata2[3]
Tbulk2=mydata2[4]
pbulkscale2=1.0/(1.0+mydata2[5]/3.0)
#PIcheck2=PIcheck2+dPde*(echeck2-1.0)*epsilon/(Pressure*PIoverPtarget2);

plt.plot(100.0*PIoverPtarget,PIcheck,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
#plt.plot(100.0*PIoverPtarget,echeck,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
plt.plot(100.0*PIoverPtarget2,PIcheck2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')
#plt.plot(100.0*PIoverPtarget2,echeck2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticklabels(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.25)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(0.8,1.5,0.04), minor=False)
ax.set_yticklabels(np.arange(0.8,1.5,0.04), minor=False)
ax.set_yticks(np.arange(0.8,1.5,0.01), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.ylim(0.95,1.07999)

text(1.2,0.99,'$\\rho_B=0$',color='green',fontsize='22',ha='right')
text(1.2,1.05,'$\\rho_B=0.1$',color='red',fontsize='22',ha='right')

plt.xlabel('$|\Pi_{\\rm hydro}/P|$ (percent)', fontsize=18, weight='normal')
plt.ylabel('$\Pi/\Pi_{\\rm hydro}$',fontsize=18)

ax = fig.add_axes([0.16,0.52,0.8,0.22])
plt.plot(100.0*PIoverPtarget,pbulkscale,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
plt.plot(100.0*PIoverPtarget2,pbulkscale2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.set_xticklabels([], minor=False)
plt.xlim(0.0,1.25)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(0.0,1.5,0.1), minor=False)
ax.set_yticklabels(np.arange(0.0,1.5,0.1), minor=False)
ax.set_yticks(np.arange(0.0,1.5,0.05), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.ylim(0.75,1.05)
plt.ylabel('$R_p$',fontsize=18)

ax = fig.add_axes([0.16,0.74,0.8,0.22])
plt.plot(100.0*PIoverPtarget,1000*Tbulk,marker='o',markersize='8',linestyle='--',linewidth=1,color='g')
plt.plot(100.0*PIoverPtarget2,1000*Tbulk2,marker='s',markersize='8',linestyle='--',linewidth=1,color='r')

ax.set_xticks(np.arange(0,1.5,0.2), minor=False)
ax.set_xticks(np.arange(0,1.5,0.05), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.set_xticklabels([], minor=False)
plt.xlim(0.0,1.25)

ax.tick_params(direction='in', which='both')
ax.yaxis.set_ticks_position('both')
ax.set_yticks(np.arange(100,220,20), minor=False)
ax.set_yticklabels(np.arange(100,220,20), minor=False)
ax.set_yticks(np.arange(100,220,10), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.ylim(145,210)
plt.ylabel('$T+\delta T$',fontsize=18)

plt.savefig('bulk.pdf',format='pdf')
os.system('open -a Preview bulk.pdf')
#plt.show()
quit()
