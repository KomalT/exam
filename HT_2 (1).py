# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 10:55:45 2016

@author: sharvari
"""
import scipy
# Tube side Fluid- OIL (Hot Fluid)******
cpo=3000         # J/kg-K
doil=850         # kg/m^3
muo=0.00136      # Pa-s
ko=.45           # W/m-K
Tino=120.0       # deg C
Touto=70.0       # deg C
Fo=1180          # LPM
Fo=Fo/(1000.0*60)# m^3/sec
mo=Fo*doil       # kg/sec


# Shell side Fluid- Water (Cold Fluid)****
cpw=4000         # J/kg-K
dw=1000          # kg/m^3
muw=0.001        # Pa-s
kw=.6            # W/m-K
Tinw=32.0        # deg C
Toutw=37.0       # deg C

# ******Heat duty******
Q=mo*cpo*(Tino-Touto)
print 'Heat duty (kW)-', Q/1000

#***** Mass flowrate of Water******
mw=Q/(cpw*(Toutw-Tinw))

#*****delT LM******
delT=((Tino-Toutw)-(Touto-Tinw))/scipy.log((Tino-Toutw)/(Touto-Tinw))
print 'delT LM',delT


#LMTD Correction factor For Multipasses
R=(Tinw-Toutw)/(Touto-Tino)
S=(Touto-Tino)/(Tinw-Tino)
print 'R,S-',R,   S
num=scipy.sqrt((R**2+1))*scipy.log(((1-S)/(1-R*S)))
den=(R-1)*scipy.log(((2-S*((R+1)-scipy.sqrt(R**2+1)))/(2-S*((R+1)+scipy.sqrt(R**2+1)))))
Ft=num/den
print 'LMTD Correction factor For Multipasses-',Ft

# assumption for U
Uass=550      # J/(m^2*s*K)


#Optimizing Velocity
v=10
F=1
#while v<1 or v>2:

L=5*.3048 # m (16 foot-Selected depending on available space)
di=22*10**-3 #m (inner Diameter)
do=24*10**-3 #m (outer Diameter)


n=1
def compute_velocity(F,di,do,L,n):
    A=Q/(Uass*delT*F) #Heat transfer Area
    print'Heat transfer Area', A
    global N
    # Number of Tubes
    N1=A/(scipy.pi*do*L)
    if N1-int(N1)>=.5:
        N=int(N1+1)
        print'Number of Tubes-', N
    else:
        N=int(N1)
        print'Number of Tubes-', int(N1)
    #print N
    # Number of passes
    if n<>1:
        N=N/n
        #print N
    # Cross flow Area
    Acs=N*scipy.pi*.25*(di)**2
    print 'Cross flow Area-',Acs
    v=Fo/(Acs)
    print 'Velocity-', v
    return v
#n=2
vel=compute_velocity(F,di,do,L,n)
n=0
while vel<1 :
    F=Ft
    n=n+2
    vel=compute_velocity(F,di,do,L,n)
print 'Number of passes-',n    
print 'Velocity-', vel
PT=1.5*do #Pitch
Clr=0.5*do #Clearance

""" Square Layout :Water is a dirty fluid """
""" For n=4 and Square Pitch """
K1=.158
n1=2.263
N=N*n #Total Number of Tubes
#print N
Db=do*(N/K1)**(1/n1) # Tube Bundle Diameter
print 'Tube Bundle Diameter-',Db
c=0.04
 # Clearance=30-40 mm
Ds=Db+c 
print 'Shell Diameter',Ds

"""Baffle -Segmental Baffle (25% Baffle Cut)"""
Dbaffle=Ds-.004 #Clearance of 4-6 mm
print 'Baffle Diameter',Dbaffle
Baffle_spacing=Ds #Baffle_spacing: .2 Ds to Ds
print'Baffle spacing=',Baffle_spacing
tcs=(Ds*Clr*Baffle_spacing)/PT
print 'Cross flow area(Shell side)', tcs
vs=mw/(dw*tcs)# Shell side velocity
print 'Shell side velocity',vs

# Tube side hi calc
Ret=di*vel*doil/(muo)
print 'Re (Tube Side)',Ret
Prt=cpo*muo/ko
print 'Pr (Tube Side)',Prt
hi=.023*((Ret)**0.8)*((Prt)**.3)*ko/do
print 'Heat Transfer Coefficient(tube side)',hi
  
#Shell Side ho calc
De=4*(PT**2-scipy.pi*.25*do**2)/(scipy.pi*do)
print 'Equivalent Diameter',De
Res=De*vs*dw/muw
print 'Re (Shell Side)',Res
Prs=cpw*muw/kw
print 'Pr (Shell side)',Prs
ho=.023*((Res)**0.8)*((Prs)**.4)*kw/De
print 'Heat Transfer Coefficient(Shell side)',ho
U=1/((do**2/(di**2*hi))+(1/ho))
print U


    
    







