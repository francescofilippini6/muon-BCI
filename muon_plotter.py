import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import labellines

n=1.33
c = 0.3 #m/ns
qx=-50
qy=-4
qz=0
#t0=0
theta=np.radians(80)
phi=0
ux=math.cos(theta)*math.cos(phi)
uy=math.cos(theta)*math.sin(phi)
uz=math.sin(theta)

#def tc(qx,qy,qz,ux,uy,uz):
def tc(ux,uy,uz):
    print("tc")
    tc=(1/c)*(((qz*uz)-(ux*qx+uy*qy+uz*qz))/(1-uz**2))
    return tc

#def zc(qx,qy,qz,ux,uy,uz):
def zc(ux,uy,uz):
    print("zc")
    zc=((qz-uz*(ux*qx+uy*qy+uz*qz))/(1-uz**2))
    return zc

#def x(t,qx,ux):
def x(t,ux):
    print("x")
    xval = qx+ux*c*t
    return(xval)

def y(t,uy):
#def y(t,qy,uy):
    print("y")
    yval = qy+uy*c*t
    return(yval)

#def d(t,qx,qy,ux,uy):
def d(t,ux,uy):
    print("d")
    return(np.sqrt(x(t,ux)**2+y(t,uy)**2))
#return(np.sqrt(x(t,qx,ux)**2+y(t,qy,uy)**2))

#def d_gamma(z,qx,qy,qz,ux,uy,uz):
def d_gamma(z,ux,uy,uz):
    print("d_gamma")
    gamma_distance = ((n/np.sqrt(n**2-1))*np.sqrt((1-uz**2)*((z-zc(ux,uy,uz))/uz)**2+d(tc(ux,uy,uz),ux,uy)**2))
    return gamma_distance

#def t_gamma(z,qx,qy,qz,ux,uy,uz):
def t_gamma(z,ux,uy,uz,t0):
    print("t_gamma")
    aaa=[]
    for za in z:
        timing_of_gamma = t0+(tc(ux,uy,uz))+(1/c)*((za-zc(ux,uy,uz))/uz+((uz*n**2-1)/(n*uz))*d_gamma(za,ux,uy,uz))
        aaa.append(timing_of_gamma)
    return aaa

#def inverse(t,qx,qy,qz,ux,uy,uz,t0):
#    print("inverse")
#    return(uz*(c*(t_gamma(t,qx,qy,qz,ux,uy,uz,t0)-tc(qx,qy,qz,t0,ux,uy,uz))-d_gamma(ta,qx,qy,qz,ux,uy,uz,t0)*(uz*n**2-1/uz*n)))+zc(qx,qy,qz,ux,uy,uz)

#def test(x,a,b):
#    return np.power(x,2)/a**2

heightA=[0,65,101,137,173,209,245,281,317,353,389,425,461,497,533,569,605,641,677]
height=[101,137,173,209,281,317,353,389]

dom_time=[55467339.,55467449.,55467580.,55467704.,55467945.,55468065.,55468191.,55468321.,]
area=np.array([1,1,2,3,10,5,2,1])*50
#time=np.linspace(55467300,55470000,100)

#params, _ =curve_fit(t_gamma,dom_time,height)
#print(params)


name=['base','dom 1','dom 2','dom 3','dom 4','dom 5','dom 6','dom 7','dom 8','dom 9','dom 10','dom 11','dom 12','dom 13','dom 14','dom 15','dom 16','dom 17','dom 18']
f=plt.figure()
ax=f.add_subplot(111)
xvals=[]
hits=ax.scatter(dom_time,height,s=area, alpha=0.5)
for i in range(len(heightA)):
    ax.axhline(y=heightA[i], color='k', linestyle='--',alpha=0.1,label=name[i],linewidth=1)
    #print(name[i])
    xvals.append(dom_time[-1]+70)
labellines.labelLines(ax.get_lines(),xvals=xvals)
#ax.text(0.02,1+i*38, name[i],fontsize=5)
ax.set_xlabel('time (ns)')
ax.set_ylabel('height z (m)')
ax.plot(t_gamma(height,ux,uy,uz,dom_time[0]-415),height)
#ax.legend((hits), ('hits'))
#plt.plot(t_gamma(height,time,qx,qy,qz,ux,uy,uz,t0),height)
plt.xlim([dom_time[0]-100,dom_time[-1]+150])
plt.ylim([-10,700])
plt.show()
