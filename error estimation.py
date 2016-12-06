import scipy
import numpy as np
import scipy,scipy.optimize
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import array
'''
adsorption
'''
'''
x=array([0.0232,0.0318,0.0472,0.0510,0.0530,0.0616,0.0626,0.0664,0.0741]) #ceq of nacl
y=array([0.00154,0.00170,0.00176,0.00196,0.00235,0.00256,0.00225,0.00252,0.00259])#x/m
Nnaoh=0.048# N
Vnaoh=array([8,7.1,5.5,5.1,4.9,4,3.9,3.5,2.7])#in ml
dVnaoh=1#ml
err=Nnaoh*(dVnaoh/Vnaoh)/1000 # error in measurement
yerr=err*y
'''

import pandas as pd
location='C:\Users\BCS\Trial_import (1).xlsx'
df=pd.read_excel(location,0)
Xdata=df['X']
x=np.array([Xdata])[0]
print 'x values'
print x
YData=df['Y']
y=np.array([YData])[0]
print 'y values'
print y
Vdata=df['V']
V=np.array([Vdata])[0]
print 'volm of NaOH'
print V
Nnaoh=0.048# N
dVnaoh=0.1#ml
print 'vol of NaCl in ml'
Mdata=df['M']
M=np.array([Mdata])[0]
print M

delx=(M*10**(-3)*Nnaoh*dVnaoh/(2*1000)) # error in measurement of strength
delm=0.01 # error in mass
yerr=scipy.sqrt(delx**2+delm**2)# total err in y


'''
Fruendlich isotherm
x/m=k*c^(1/n)
'''
print 'Fruendlich isotherm'
def func(x, a, b):
     
     return a*x**(1/b)
ydata = func(x,2,0.1)
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
pguess=[2,0.1]
res=scipy.optimize.curve_fit(func,x,y,[2,0.1],yerr)
P=res[0]
print("k,n=")
print P
Pcov=res[1]
#print Pcov

ycalc=func(x, P[0],P[1])
r2=get_r2(x,y,ycalc)
print 'r2'
print r2
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(x,y,'ro');
ax.plot(x,ycalc,'b')
plt.xlabel('x/m (gm equivalent/gm)')
plt.ylabel('c (N)')
ax.title.set_text('Fruendlich isothem r2=%f'%(r2))
fig.canvas.draw()
plt.show()

'''
Langmuir isotherm
x/m=k1*c/(k2*c+1)
'''
print 'Langmuir isotherm'
def curve2(x,C,D):
    
    y=C*x/(1+D*x)
    return y

pguess=[2,0.1]
res=scipy.optimize.curve_fit(curve2,x,y,[2,0.1],yerr)
P=res[0]
print("k1,k2=")
print P
Pcov=res[1]
#print Pcov
ycalc=curve2(x, P[0],P[1])
r2=get_r2(x,y,ycalc)
print 'r2'
print r2
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(x,y,'ro');
ax.plot(x,ycalc,'b')
plt.xlabel('x/m (gm equivalent/gm)')
plt.ylabel('c (N)')
ax.title.set_text('Langmuir isothem r2=%f'%(r2))
fig.canvas.draw()
plt.show()

