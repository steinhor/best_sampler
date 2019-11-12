import matplotlib.pyplot as plt
import numpy as np
import os

PI=3.14159265358979323

x2=[]
y2=[]

x3=[]
y3=[]

x2 , y2 = np.loadtxt('../build/corrections.txt', delimiter=' ', unpack=True)
x3 , y3 = np.loadtxt('../build/no_corrections.txt', delimiter=' ', unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

x = np.linspace(0, 1.0, 100)
T=.133294
m=.138
ax1.plot(x, np.exp(-np.sqrt(x**2+m*m)/T)/(1-np.exp(-np.sqrt(x**2+m*m)/T)))
ax1.plot(x, np.exp(-np.sqrt(x**2+m*m)/T))

ax1.scatter(x2,y2,c='r')
ax1.scatter(x3,y3,c='g')
plt.ylim(ymin=0)
plt.xlim(left=0,right=1.0)
plt.xlabel('p')
plt.ylabel('f')
plt.show()

os.remove("../build/corrections.txt")
os.remove("../build/no_corrections.txt")
