# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 15:01:32 2016

@author: Sharvari
"""
import matplotlib.pyplot as plt
from scipy import array
import scipy.optimize
Mdm=0.04175 #molar density gmol/Lit of methane ref 3
Hm=714.2857 #henry's coef of methane atm*lit/gmol ref2
Mdc=0.04186 #molar density gmol/Lit of co2 ref 3
Hc=29.76 #henry's coef of CO2 atm*lit/gmol ref 3
pwsat=0.03125 #atm ref3
pT=1 #atm
kla=0.008 #sec-1 ref1
Lc=0#initial guess for co2 in liq
Lm=0#initial guess for methane in liq
Lw=100#initial guess for water in liq
Gc=50#known co2 flow rate
Gm=50#known methane flow rate
Gw=0#known initial water in gas
L=10#meter
R=0.0821#atm/molK
T=298#K

Lt=[Lc,Lm,Lw]
Gt=[Gc,Gm,Gw]


kga=7*10^-7#atm-1*sec-1ref 4
Lc=0#initial guess for co2 in liq;
Lm=0#initial guess for methane in liq
Lw=100#initial guess for water in liq
Gc=50#known co2 flow rate
Gm=50#known methane flow rate
Gw=0#known initial water in gas
Gc1k3=[]
Gm1k3=[]
Gw1k3=[]
Lc1k3=[]
Lm1k3=[]
Lw1k3=[]
Lexp=[0,0,100]#known values of liquid flow rates at z=0
def absorption(Gc,Gm,Gw,Lc,Lm,Lw):
    for i in range (10):
        h=L/(10-1)
        yc=Gc/(Gc+Gm+Gw)
        ym=Gm/(Gc+Gm+Gw)
        yw=Gw/(Gc+Gm+Gw)
        xc=Lc/(Lc+Lm+Lw)
        xm=Lm/(Lc+Lm+Lw)
        #xw=Lw/(Lc+Lm+Lw)
        Gc1=Gc-h*kla*(pT*yc/Hc-xc*Mdc)
        Gm1=Gm-h*kla*(pT*ym/Hm-xm*Mdm)
        Gw1=Gw+h*kga*(pT*yw-pwsat)
        Lc1=Lc-h*kla*(pT*yc/Hc-xc*Mdc)
        Lm1=Lm-h*kla*(pT*ym/Hm-xm*Mdm)
        Lw1=Lw+h*kga*(pT*yw-pwsat)
        Gc=Gc1
        Gm=Gm1
        Gw=Gw1
        Lc=Lc1
        Lm=Lm1
        Lw=Lw1
       
        
    a=Lc;b=Lm;c=Lw;
    return array([a,b,c,Gc,Gm,Gw])
def error(x,a,b,c):
    
    a11=absorption(50,50,0,a,b,c)
    err1=a11[0]-Lexp[0]
    err2=a11[1]-Lexp[1]
    err3=a11[2]-Lexp[2]
    err=[err1,err2,err3]
    return err
Lguess=[10,10,80]
xdata=array([10])
ydata=array([0,0,0])
Lans=scipy.optimize.curve_fit(error,xdata,ydata,Lguess)
La=Lans[0] 
print La

error=(200-(La[0]+La[1]+La[2]+Gc+Gm+Gw))/200
print ("error is=")
print error
Gc=50;Gm=50;Gw=0;Lc=1.17431492e-03;Lm=4.89271499e-05;Lw=9.67744678e+01;
for i in range (10):
        h=L/(10-1)
        yc=Gc/(Gc+Gm+Gw)
        ym=Gm/(Gc+Gm+Gw)
        yw=Gw/(Gc+Gm+Gw)
        xc=Lc/(Lc+Lm+Lw)
        xm=Lm/(Lc+Lm+Lw)
        #xw=Lw/(Lc+Lm+Lw)
        Gc1=Gc-h*kla*(pT*yc/Hc-xc*Mdc)
        Gm1=Gm-h*kla*(pT*ym/Hm-xm*Mdm)
        Gw1=Gw+h*kga*(pT*yw-pwsat)
        Lc1=Lc-h*kla*(pT*yc/Hc-xc*Mdc)
        Lm1=Lm-h*kla*(pT*ym/Hm-xm*Mdm)
        Lw1=Lw+h*kga*(pT*yw-pwsat)
        Gc=Gc1
        Gm=Gm1
        Gw=Gw1
        Lc=Lc1
        Lm=Lm1
        Lw=Lw1
        Gc1k3.append(Gc)
        Gm1k3.append(Gm)
        Gw1k3.append(Gw)
        Lc1k3.append(Lc)
        Lm1k3.append(Lm)
        Lw1k3.append(Lw)

n=range(1,11)
plt.plot(n,Lw1k3)
plt.xlabel('length')
plt.ylabel('water flow rate')
plt.show()
'''Refrences:
ref1 :Advanced liq gas contactor technology for multi pollutant 
control systems.NSG.
ref2:Compilation of Henry’s Law Constants for Inorganic and
Organic Species of Potential Importance in Environmental
Chemistry
ref3:Physical Properties APPENDIX B of perry's handbook
ref4: Perry's handbook page 2113
'''



