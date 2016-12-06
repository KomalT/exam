import scipy,scipy.optimize,scipy.stats
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
'''
adsorption
'''
x=array([0.0232,0.0318,0.0472,0.0510,0.0530,0.0616,0.0626,0.0664,0.0741]) #ceq of nacl
y=array([0.00154,0.00170,0.00176,0.00196,0.00235,0.00256,0.00225,0.00252,0.00259])#x/m
Vnaoh=array([8,7.1,5.5,5.1,4.9,4,3.9,3.5,2.7])#in ml
Nnaoh=0.048# N
dVnaoh=1#ml
err=Nnaoh*(dVnaoh/Vnaoh)/1000 # error in measurement
yerr=err*y
'''
Fruendlich isotherm
x/m=k*c^(1/n)
'''
def curve1(x,p):
    [A,B]=p
    y=A*x**(1/B)
    return y
def err1(p,x,y,yerr):
    err=(y-curve1(x,p))/yerr
    return err
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2

pguess=[2,0.1]
res=scipy.optimize.leastsq(err1,pguess,args=(x,y,yerr))
print("k,n=")
print res
P=res[0]
ycalc=curve1(x,P)
r2=get_r2(x,y,ycalc)
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
def curve2(x,k):
    [C,D]=k
    y=C*x/(1+D*x)
    return y
def error(k,x,y,yerr):
    ycalc=curve2(x,k)
    err=(y-ycalc)/yerr
    return err
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
kguess=[0.5,8]
klsq=scipy.optimize.leastsq(error,kguess,args=(x,y,yerr))
k=klsq[0]
print("k1,k2=")
print k
ycalc=curve2(x,k)
r2=get_r2(x,y,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(x,y,'ro');
ax.plot(x,ycalc,'b')
plt.xlabel('x/m (gm equivalent/gm)')
plt.ylabel('c (N)')
ax.title.set_text('Langmuir isotherm r2=%f '%(r2))
fig.canvas.draw()
plt.show()

   
