import scipy as sc
import scipy.optimize as optimization
import os
import pandas
from scipy.integrate import quad
import scipy as sc
from scipy.integrate import odeint
from scipy import array

import perrysdata
import eosClass
Ethanol=perrysdata.Compound('Ethanol')
h20=perrysdata.Compound('Water')
print Ethanol.CpIG(300)
def derivative(y,t):
    dy1=-(y[0]-y[1])#/Ethanol.CpIG(y[0])
    dy2=-(y[0]-y[1])#/h20.CpIG(y[1])
    
    return([dy1,dy2])
t=sc.linspace(0,9,10)
ydata=300
def err(a,x):
  yinitial=([400,a])
  y=odeint(derivative,yinitial,t)
  
  return y[9,1]
xdata=10
([ans,err])=optimization.curve_fit(err,xdata,ydata,354)
#ans=fsolve(err,350)
print ans[0]
