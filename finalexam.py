# -*- coding: utf-8 -*-
"""
Created on Sun May 01 13:03:15 2016

@author: Komal
"""

from scipy.integrate import odeint
from scipy import linspace
from scipy import array
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import scipy,scipy.optimize
import win32com.client
xl=win32com.client.gencache.EnsureDispatch("Excel.Application")
wb=xl.Workbooks('ExamProblemData2.1.xlsx')
sheet=wb.Sheets('ExamProblemData')
def getdat(sheet,Range):
    data=sheet.Range(Range).Value
    data=scipy.array(data)
    data=data.reshape((1,len(data)))[0]
    return data
t=getdat(sheet,"A2:A11")   #****** Re
y1=getdat(sheet,"B2:B11")   #******* f
y2=getdat(sheet,"C2:C11") 
print t
#t=array([0,9,18,26,35,44,53,62,71,80])
#y1=array([49.44,41.17,33.54,27.54,20.96,14.93,10.02,7.71,4.8,1.88])
#y2=array([0,1.14,2.88,3.44,4.13,4.48,4.75,5.57,5.75,5.8])

y=array([y1,y2])
#y=np.reshape(10,2)
y=np.transpose(y)
yinitial=array([49.44,0])
#3print y
Z=array([y1,y2,t])
#print Z
y1,y2,t=Z
print Z
def derivative(y,t,k1,k2,m,n):
    #k1,k2,m,n=A
    dcA=-k1*y[0]**n+k2*y[1]**m
    dcB=-k2*y[1]**m+k1*y[0]**n
    return array([dcA,dcB])
A0=array([0.5,0.5,1,1])
k1,k2,m,n=A0
yk=odeint(derivative,yinitial,t,args=(k1,k2,m,n))
#print yk
def error(A,Z):
    k1=A[0]
    k2=A[1]
    m=A[2]
    n=A[3]
    #y1,y2,t=Z
    y1=Z[0]
    y2=Z[1]
    t=Z[2]
    y=array([y1,y2])
    #t=array([Z[:,2]])
    #A=array([k1,k2,m,n])
    y3=odeint(derivative,yinitial,t,args=(k1,k2,m,n))
    y3=np.transpose(y3)
    #y30=y3[:,0]
    #y31=y3[:,1]
    #ya0=y1-y30
    #yb0=y2-y31
    #3print y
    r=y-y3
    return r[0]
([ans,err])=scipy.optimize.leastsq(error,A0,Z)
#print ans

k,l,m,n=ans
yk=odeint(derivative,yinitial,t,args=(k1,k2,m,n))
print yk[0]
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(t,y1,'*')
ax.plot(t,yk[:,0],'b')    
ax.title.set_text('P Vs x3')
fig.canvas.draw()
plt.show()

