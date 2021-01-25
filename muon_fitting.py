import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import labellines
from lmfit import Model

n=1.33
c = 0.3 #m/ns
#qx=-50
#qy=-4
#qz=0
#t0=0
#theta=np.radians(80)
#phi=0
#ux=math.cos(theta)*math.cos(phi)
#uy=math.cos(theta)*math.sin(phi)
#uz=math.sin(theta)

def tc(qx,qy,qz,ux,uy,uz,t0):
    print("tc")
    tc=t0+(1/c)*(((qz*uz)-(ux*qx+uy*qy+uz*qz))/(1-uz**2))
    return tc

def zc(qx,qy,qz,ux,uy,uz,t0):
    print("zc")
    zc=((qz-uz*(ux*qx+uy*qy+uz*qz))/(1-uz**2))
    return zc

def x(t,qx,ux,t0):
    print("x")
    xval = qx+ux*c*(t-t0)
    return(xval)


def y(t,qy,uy,t0):
    print("y")
    yval = qy+uy*c*(t-t0)
    return(yval)

def d(t,qx,qy,ux,uy,t0):
    print("d")
    return(np.sqrt(x(t,qx,ux,t0)**2+y(t,qy,uy,t0)**2))

def d_gamma(z,qx,qy,qz,ux,uy,uz,t0):
    print("d_gamma")
    gamma_distance = ((n/np.sqrt(n**2-1))*np.sqrt((1-uz**2)*((z-zc(qx,qy,qz,ux,uy,uz,t0))/uz)**2+d(tc(qx,qy,qz,ux,uy,uz,t0),qx,qy,ux,uy,t0)**2))
    return gamma_distance

def t_gamma(z,qx,qy,qz,ux,uy,uz,t0):
    timing_of_gamma = (tc(qx,qy,qz,ux,uy,uz,t0)-t0)+(1/c)*((z-zc(qx,qy,qz,ux,uy,uz,t0))/uz+((uz*n**2-1)/(n*uz))*d_gamma(z,qx,qy,qz,ux,uy,uz,t0))
    return timing_of_gamma

heightA=[0,65,101,137,173,209,245,281,317,353,389,425,461,497,533,569,605,641,677]
height=np.array([101,137,173,209,281,317,353,389])

dom_time=np.array([55467339.,55467449.,55467580.,55467704.,55467945.,55468065.,55468191.,55468321.,])
area=np.array([1,1,2,3,10,5,2,1])*50

hmodel=Model(t_gamma)
params=hmodel.make_params(qx=-40,qy=-2,qz=5,ux=0.1,uy=0.05,uz=0.95,t0=dom_time[0])
#params['qx'].min=-50.1
#params['qx'].max=-49.9
#params['qy'].min=-3.99
#params['qy'].max=-4.01
params['qz'].min=-0.01
params['qz'].max=0.01

params['ux'].min=0
params['uy'].min=-1
params['uz'].min=0
params['ux'].max=1
params['uy'].max=1
params['uz'].max=1


#weights 1/error on my data
result = hmodel.fit(dom_time, params, z=height,weights=1)
print(result.fit_report())
print("Theta:",np.degrees(np.arcsin(result.params['uz'].value)))

name=['base','dom 1','dom 2','dom 3','dom 4','dom 5','dom 6','dom 7','dom 8','dom 9','dom 10','dom 11','dom 12','dom 13','dom 14','dom 15','dom 16','dom 17','dom 18']
f=plt.figure()
ax=f.add_subplot(111)
xvals=[]
hits=ax.scatter(dom_time,height,s=area, alpha=0.5,label='triggered hits')
for i in range(len(heightA)):
    ax.axhline(y=heightA[i], color='k', linestyle='--',alpha=0.1,label=name[i],linewidth=1)
    xvals.append(dom_time[-1]+70)
    
labellines.labelLines(ax.get_lines(),xvals=xvals)

ax.set_xlabel('time (ns)')
ax.set_ylabel('height z (m)')
#ax.axhline(y=zc(result.params['qx'].value,result.params['qy'].value,result.params['qz'].value,result.params['ux'].value,result.params['uy'].value,result.params['uz'].value,result.params['t0'].value), color='r', linestyle='-',alpha=0.1,linewidth=1)
#plt.errorbar(x, y, xerr=0.2, yerr=0.4)
fit_result=ax.plot(result.best_fit,height,label='fit result')
#ax.legend((hits),('triggered hits'))
#ax.set_yscale('log')
plt.xlim([dom_time[0]-100,dom_time[-1]+150])
plt.ylim([-10,700])
plt.show()
