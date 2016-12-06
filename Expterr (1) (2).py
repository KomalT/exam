# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 22:47:15 2016

@author: user
"""
import pandas as pd
import scipy,scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
Le=2.52  # m empty pipe
Lm=2.69  # m with mixing element
De=.0209 # m empty
Dm=.0198 # m with mixing element
Row=1000 # Kg/m3
RoHg=13600 #Kg/m3
RoC=1487 # Kg/m3
Meu=0.001 # Pa.s
DelRo=1   # Kg/m3
DelH=0.001# m
DelD=2*10**(-5)#m
DelL=0.002# m
DelQ=0.002# 
#location='C:\Users\user\Documents\fre1.xlsx'
#df=pd.read_excel(location,0)#C:\Users\user\Documents
'''
xlsx=pd.ExcelFile('C:\Users\BCS\kt(3).xlsx')
sheet1=xlsx.parse(0)
Qdata=sheet1.icol(3)
Q=np.array([Qdata])[0]
print Q
Hdata=sheet1.icol(2)
H=np.array([Hdata])[0]
print H
xdata=sheet1.icol(0)
x=np.array([xdata])[0]
print x
ydata=sheet1.icol(1)
y=np.array([ydata])[0]
print y
#Qdata=df['Q']
#print Qdata
#Q=np.array([Qdata])[0]
#print Q[0,0]
#print Q
#Hdata=df['H']
#H=np.array([Hdata])[0]
#Redata=df['RE']
#x=np.array([Redata])[0]
#fdata=df['F']
#y=np.array([fdata])[0]
'''
Q=array([19.104*10**-6,27.3*10**-6,32*10**-6,57.85*10**-6,70.87*10**-6,87.88*10**-6])
H=array([0.5*10**-2,1*10**-2,1.1*10**-2,1.8*10**-2,2.5*10**-2,3.5*10**-2])
x=array([1164.42,1663.97,1950.45,3526.04,4319.63,5356.41])
y=array([.03188,.03122,.025,.0125,.01158,.01055])
#constant=(2*(DelRo/(RoC-Row))**2)+((DelRo/Row)**2)+((DelL/Le)**2)+(5*DelD/De)**2
#ErrbyY=(9.81*3.14**2/32)*np.sqrt(constant+(DelH**2)/np.square(H)+((2*DelQ)**2)/np.square(Q))
e1=(9.81*3.14**2*De**5/(32*Row*Le))*np.divide(H,np.square(Q))*DelRo
#print e1
e2=(9.81*3.14**2*De**5/(32*Row*Le))*np.divide(H,np.square(Q))*RoHg*(-1/Row**2)*DelRo
e3=((RoHg-Row)*9.81*3.14**2*De**5/(32*Row*Le))*np.divide(1,np.square(Q))*DelH
e4=((RoHg-Row)*9.81*3.14**2*De**5/(32*Row))*np.divide(H,np.square(Q))*(-1/Le**2)*DelL
e5=((RoHg-Row)*9.81*3.14**2/(32*Row*Le))*np.divide(H,np.square(Q))*5*De**4*DelD
e6=((RoHg-Row)*9.81*3.14**2*De**5/(32*Row*Le))*np.divide(H,np.power(Q,3))*(-2*DelQ)
yerr=np.sqrt(np.square(e1)+np.square(e2)+np.square(e3)+np.square(e4)+np.square(e5)+np.square(e6))
#yerr=ErrbyY
def curve(x,A,B):
    #[A,B]=p
    
    return A*x**B#def err1(p,x,y,yerr):
    #err=(y-curve(x,p))/yerr
    #return err
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
pguess=[16,-1]
res=scipy.optimize.curve_fit(curve,x,y,[16,-1],yerr)
P=res[0]
print P
pcov=res[1]
print pcov
print pcov.diagonal()**.5
[A,B]=P
resid = (y-curve(x,A,B))/yerr
chisq=scipy.sum(resid**2)
chisqred = chisq/(len(resid)-len(P))
print chisqred
[A,B]=P
ycalc=curve(x,P[0],P[1])
r2=get_r2(x,y,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(x,y,'*')
ax.plot(x,ycalc,'b')    
ax.title.set_text('R2=%f'%(r2))
fig.canvas.draw()
plt.show()  

    
    
    
    
    
    