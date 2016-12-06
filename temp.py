# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 07 15:52:32 2016

@author: Komal
"""

# Code for flash column
from scipy.optimize import fsolve 
import os
import pandas
from scipy.integrate import quad
import scipy as sc

#Here we are importing excel file containing Data from Perry's hand book
filename = 'Data.xlsx' 
thisdir = os.getcwd()
if os.path.isdir("Data"): 
    os.chdir("Data")
xl_file = pandas.ExcelFile(filename) 
os.chdir(thisdir)
df = {}
df["Critical"] = xl_file.parse("Critical") 
df["LiqDens"] = xl_file.parse("LiqDens") 
df["Pvap"] = xl_file.parse("Pvap") 
df["Hvap"] = xl_file.parse("Hvap") 
df["CpIG"] = xl_file.parse("CpIG") 
df["HandGf"] = xl_file.parse("HandGf")


class Compound:
    def __init__(self, name): 
        self.Name = name 
        self.getCriticalConstants() 
        self.getHandGofFormation()
        self.constCp = self.getConstants('CpIG') 
        self.constPvap = self.getConstants('Pvap')
        self.constHvap = self.getConstants('Hvap') 
        self.constLiqDens = self.getConstants('LiqDens')
    def getCriticalConstants(self): 
        name = self.Name
        db = df['Critical']
        comp = db[db.Name == name] 
        self.MW = comp.MolWt.values[0] 
        self.Tc = comp.Tc.values[0] 
        self.Pc = comp.Pc.values[0] 
        self.Vc = comp.Vc.values[0] 
        self.Zc = comp.Zc.values[0] 
        self.acc = comp.Acc.values[0]
    def getHandGofFormation(self): 
        name = self.Name
        db = df['HandGf']
        comp = db[db.Name == name] 
        self.Hf = comp.Hf.values[0] 
        self.Gf = comp.Gf.values[0] 
        self.Tf = 298.15 #K
        self.Pf = 1.0e5 #Pa
        self.Sf = (self.Hf - self.Gf)/self.Tf 
        self.Sfabinit = comp.Sf.values[0] 
        self.Hcomb = comp.Hcomb.values[0]
    def getConstants(self, sheet):
        name = self.Name
        cnsts = ['C1','C2','C3','C4','C5'] 
        constants = []
        db = df[sheet]
        comp = db[db.Name == name]
        for cnst in cnsts:
            try:
                constants.append(comp[cnst].values[0]) 
            except KeyError:
                            pass 
        return constants
    def Cp(self, T):
        [C1, C2, C3, C4, C5] = self.constCp
        cp = C1 + C2*((C3/T)/sc.sinh(C3/T))**2 + C4*((C5/T)/sc.cosh(C5/T))**2 
        return cp
    def Pvap(self, T):
        [C1, C2, C3, C4, C5] = self.constPvap 
        p = C1 + C2/T + C3*sc.log(T) + C4*T**C5 
        return sc.exp(p)*1.1
    def Tsat(self,P):
        def f(T):
            f=P-self.Pvap(T) 
            return f
        z = fsolve(f,500)[0] 
        return z
    def Hvap(self, T):
        [C1, C2, C3, C4] = self.constHvap
        Tr = T/self.Tc
        exponent = C2 + C3*Tr + C4*Tr**2 
        hvap = C1*(1-Tr)**exponent 
        return hvap
    def getallnames():
        db = df["Critical"]
        listNames = db.Name.to_dict().values() 
        return listNames
    #if __name__ == "__main__": 
        #met = Compound('Methane')
#Getting constants for the required compunds 
Ethanol = Compound('Ethanol')
Acetaldehyde = Compound('Acetaldehyde') 
Hydrogen = Compound('Hydrogen') 
Methane=Compound('Methane')
Carbonmonoxide=Compound('Carbon monoxide')
Carbondioxide=Compound('Carbon dioxide')
Ac=Compound('Acetic acid')
EtOEt=Compound('Diethyl ether')
EtOAc=Compound('Ethyl acetate')
h20=Compound('Water')


#Defining Feed conditions 
z1=0.3217 #Ethanol 
z2=0.2852 #water 
z3=0.2993 #ethylene 
z4=0.002511 #butylene 
z5=0.001256
z6=0.001256
z7=0.001722
z8=0.002052
z9=0.00243
z10=0.08249

F=195

P=1.49*1.01325*10**5
Ts=25+273.16
Tf=312.8+273.16

#Feed mole fractions 
#x1=Fr1/Frt 
#x2=Fr2/Frt 
#x3=Fr3/Frt 
#x4=Fr4/Frt
#print 'feed mol frac', x1,x2,x3,x4

#Calculation of bubble point and dew point 
def Tbub(T):
    #K1=Eto.Pvap(T)/P
    #K2=h2o.Pvap(T)/P
    #K3=(Eth.Pvap(T))/P 
    #K4= but.Pvap(T)/P
    K1=Ethanol.Pvap(T)/P
    K2=Acetaldehyde.Pvap(T)/P
    K3=(Hydrogen.Pvap(T))/P 
    K4= (Methane.Pvap(T)*1)/P  
    K5= (Carbonmonoxide.Pvap(T)*1)/P  
    K6= (Carbondioxide.Pvap(T)*1)/P  
    K7= (Ac.Pvap(T)*1)/P  
    K8= (EtOEt.Pvap(T)*1)/P  
    K9= (EtOAc.Pvap(T)*1)/P  
    K10= (h20.Pvap(T)*1)/P  
    
    a=z2*K2+z3*K3+z4*K4+z1*K1+z5*K5+z6*K6+z7*K7+z8*K8+z9*K9+z10*K10
    return a-1
Tbu=float(fsolve(Tbub,300)) 
print 'Tbub',Tbu-273.15

def Tdew(T): 
    #K1=Eto.Pvap(T)/P 
    #K2=h2o.Pvap(T)/P 
    #K3=(Eth.Pvap(T))/P 
    #K4= but.Pvap(T)/P
    K1=Ethanol.Pvap(T)/P
    K2=Acetaldehyde.Pvap(T)/P
    K3=(Hydrogen.Pvap(T))/P 
    K4= (Methane.Pvap(T)*1)/P  
    K5= (Carbonmonoxide.Pvap(T)*1)/P  
    K6= (Carbondioxide.Pvap(T)*1)/P  
    K7= (Ac.Pvap(T)*1)/P  
    K8= (EtOEt.Pvap(T)*1)/P  
    K9= (EtOAc.Pvap(T)*1)/P  
    K10= (h20.Pvap(T)*1)/P  
    
    a=z2/K2+z3/K3+z4/K4+z1/K1+z5/K5+z6/K6+z7/K7+z8/K8+z9/K9+z10/K10
    return a-1
Tde=float(fsolve(Tdew,300.0)) 
print 'Tdew',Tde-273.15

#Solving the equation for flash column 
def funt(var):
    (f,T)=var
    K1=Ethanol.Pvap(T)/P
    K2=Acetaldehyde.Pvap(T)/P
    K3=(Hydrogen.Pvap(T))/P 
    K4= (Methane.Pvap(T)*1)/P  
    K5= (Carbonmonoxide.Pvap(T)*1)/P  
    K6= (Carbondioxide.Pvap(T)*1)/P  
    K7= (Ac.Pvap(T)*1)/P  
    K8= (EtOEt.Pvap(T)*1)/P  
    K9= (EtOAc.Pvap(T)*1)/P  
    K10= (h20.Pvap(T)*1)/P  
    
    rhs=(z1/(f*(K1-1)+1))+(z2/(f*(K2-1)+1))+(z3/(f*(K3-1)+1))+(z4/(f*(K4-1)+1))+(z5/(f*(K5-1)+1))+(z6/(f*(K6-1)+1))+(z7/(f*(K7-1)+1))+(z8/(f*(K8-1)+1))+(z9/(f*(K9-1)+1))+(z10/(f*(K10-1)+1)) -1
    #K1=Eto.Pvap(T)/P 
    xb1=(z1/(f*(K1-1)+1)) 
    yt1=xb1*K1
    #K2=h2o.Pvap(T)/P 
    xb2=(z2/(f*(K2-1)+1)) 
    yt2=xb2*K2
    #K3=(Eth.Pvap(T))/P 
    xb3=(z3/(f*(K3-1)+1)) 
    yt3=xb3*K3
    #K4=(but.Pvap(T)*10)/P 
    xb4=(z4/(f*(K4-1)+1)) 
    yt4=xb4*K4
    
    xb5=(z5/(f*(K5-1)+1)) 
    yt5=xb5*K5
    
    xb6=(z6/(f*(K6-1)+1)) 
    yt6=xb6*K6
    
    xb7=(z7/(f*(K7-1)+1)) 
    yt7=xb7*K7
    
    xb8=(z8/(f*(K8-1)+1)) 
    yt8=xb8*K8
    
    xb9=(z9/(f*(K9-1)+1)) 
    yt9=xb9*K9
    xb10=(z10/(f*(K10-1)+1)) 
    yt10=xb10*K10
    
    Hf1=((quad(Ethanol.Cp,Ts,Tf)[0]+Ethanol.Hvap(Ts))*z1)/1000
    Hf2=((quad(Acetaldehyde.Cp,Ts,Tf)[0]+Acetaldehyde.Hvap(Ts))*z2)/1000
    Hf7=((quad(Ac.Cp,Ts,Tf)[0]+Ac.Hvap(Ts))*z7)/1000
    Hf8=((quad(EtOEt.Cp,Ts,Tf)[0]+EtOEt.Hvap(Ts))*z8)/1000
    Hf9=((quad(EtOAc.Cp,Ts,Tf)[0]+EtOAc.Hvap(Ts))*z9)/1000
    Hf10=((quad(h20.Cp,Ts,Tf)[0]+h20.Hvap(Ts))*z10)/1000
    #Hf13=((quad(Bt.Cp,Ts,Tf)[0]+Bt.Hvap(Ts)-Bt.Hvap(Tf))*z13)/1000
    Hf3=((quad(Hydrogen.Cp,Ts,Tf)[0])*z3)/1000
    Hf4=((quad(Methane.Cp,Ts,Tf)[0])*z4)/1000
    Hf5=((quad(Carbonmonoxide.Cp,Ts,Tf)[0])*z5)/1000
    Hf6=((quad(Carbondioxide.Cp,Ts,Tf)[0])*z6)/1000
    #Hf11=((quad(Eth.Cp,Ts,Tf)[0])*z11)/1000
    #Hf12=((quad(Hex.Cp,Ts,Tf)[0])*z12)/1000
    #Hf15=((quad(BD.Cp,Ts,Tf)[0])*z15)/1000
    #Hf14=((quad(But.Cp,Ts,Tf)[0])*z14)/1000
    r1=Hf1+Hf2+Hf3+Hf4+Hf5+Hf6+Hf7+Hf8+Hf9+Hf10 #+Hf11+Hf12+Hf13+Hf14+Hf15
    
    HL1=((quad(Ethanol.Cp,Ts,T)[0]+Ethanol.Hvap(Ts)-Ethanol.Hvap(T)))/1000 
    HL2=((quad(Acetaldehyde.Cp,Ts,T)[0]+Acetaldehyde.Hvap(Ts)-Acetaldehyde.Hvap(T)))/1000 
    HL7=((quad(Ac.Cp,Ts,T)[0]+Ac.Hvap(Ts)-Ac.Hvap(T)))/1000 
    HL8=((quad(EtOEt.Cp,Ts,T)[0]+EtOEt.Hvap(Ts)-EtOEt.Hvap(T)))/1000 
    HL9=((quad(EtOAc.Cp,Ts,T)[0]+EtOAc.Hvap(Ts)-EtOAc.Hvap(T)))/1000 
    HL10=((quad(h20.Cp,Ts,T)[0]+h20.Hvap(Ts)-h20.Hvap(T)))/1000 
    #HL13=((quad(Bt.Cp,Ts,T)[0]+Bt.Hvap(Ts)-Bt.Hvap(T)))/1000 
    HL3=((quad(Hydrogen.Cp,Ts,T)[0])-Hydrogen.Hvap(T))/1000 
    HL4=((quad(Methane.Cp,Ts,T)[0])-Methane.Hvap(T))/1000 
    HL5=((quad(Carbonmonoxide.Cp,Ts,T)[0])-Carbonmonoxide.Hvap(T))/1000 
    HL6=((quad(Carbondioxide.Cp,Ts,T)[0])-Carbondioxide.Hvap(T))/1000 
    #HL11=((quad(Eth.Cp,Ts,T)[0])-Eth.Hvap(T))/1000 
    #HL12=((quad(Hex.Cp,Ts,T)[0])-Hex.Hvap(T))/1000 
    #HL15=((quad(BD.Cp,Ts,T)[0])-BD.Hvap(T))/1000 
    #HL14=((quad(But.Cp,Ts,T)[0])-But.Hvap(T))/1000 

    #Heth2=((quad(Eth.Cp,Ts,T)[0])-h2o.Hvap(T))/1000 
    r2=HL1*xb1+HL2*xb2+HL3*xb3+HL4*xb4+HL5*xb5+HL6*xb6+HL7*xb7+HL8*xb8+HL9*xb9+HL10*xb10  #+HL11*xb11+HL12*xb12+HL13*xb13+HL14*xb14+HL15*xb15
    
    HV1=((quad(Ethanol.Cp,Ts,T)[0]+Ethanol.Hvap(Ts)))/1000
    HV2=((quad(Acetaldehyde.Cp,Ts,T)[0]+Acetaldehyde.Hvap(Ts)))/1000
    HV7=((quad(Ac.Cp,Ts,T)[0]+Ac.Hvap(Ts)))/1000
    HV8=((quad(EtOEt.Cp,Ts,T)[0]+EtOEt.Hvap(Ts)))/1000
    HV9=((quad(EtOAc.Cp,Ts,T)[0]+EtOAc.Hvap(Ts)))/1000
    HV10=((quad(h20.Cp,Ts,T)[0]+h20.Hvap(Ts)))/1000
    #HV13=((quad(Bt.Cp,Ts,T)[0]+Bt.Hvap(Ts)))/1000
    HV3=((quad(Hydrogen.Cp,Ts,T)[0]))/1000
    HV4=((quad(Methane.Cp,Ts,T)[0]))/1000
    HV5=((quad(Carbonmonoxide.Cp,Ts,T)[0]))/1000
    HV6=((quad(Carbondioxide.Cp,Ts,T)[0]))/1000
    #HV11=((quad(Eth.Cp,Ts,T)[0]))/1000
    #HV12=((quad(Hex.Cp,Ts,T)[0]))/1000
    #HV15=((quad(BD.Cp,Ts,T)[0]))/1000
    #HV14=((quad(But.Cp,Ts,T)[0]))/1000
    r3=HV1*yt1+HV2*yt2+HV3*yt3+HV4*yt4+HV5*yt5+HV6*yt6+HV7*yt7+HV8*yt8+HV9*yt9+HV10*yt10 #+HV11*yt11+HV12*yt12+HV13*yt13+HV14*yt14+HV15*yt15
    b=(r1)-r3*f*F-r2*(1-f)*F-3632749 
    return [rhs,b]
sol=(fsolve(funt,(0.9,320.0))) 
print sol
f=sol[0]
T=sol[1]
print 'fraction of feed vaporized',f
print 'Temp',T-273.15
Top=f*F
print 'T',Top 
Bottom=F-Top 
print 'B',Bottom

K1=Ethanol.Pvap(T)/P
K2=Acetaldehyde.Pvap(T)/P
K3=(Hydrogen.Pvap(T))/P 
K4= (Methane.Pvap(T)*1)/P  
K5= (Carbonmonoxide.Pvap(T)*1)/P  
K6= (Carbondioxide.Pvap(T)*1)/P  
K7= (Ac.Pvap(T)*1)/P  
K8= (EtOEt.Pvap(T)*1)/P  
K9= (EtOAc.Pvap(T)*1)/P  
K10= (h20.Pvap(T)*1)/P  
 
xb1=(z1/(f*(K1-1)+1)) 
yt1=xb1*K1
#K2=h2o.Pvap(T)/P 
xb2=(z2/(f*(K2-1)+1)) 
yt2=xb2*K2
#K3=(Eth.Pvap(T))/P 
xb3=(z3/(f*(K3-1)+1)) 
yt3=xb3*K3
 #K4=(but.Pvap(T)*10)/P 
xb4=(z4/(f*(K4-1)+1)) 
yt4=xb4*K4
    
xb5=(z5/(f*(K5-1)+1)) 
yt5=xb5*K5
    
xb6=(z6/(f*(K6-1)+1)) 
yt6=xb6*K6
    
xb7=(z7/(f*(K7-1)+1)) 
yt7=xb7*K7
    
xb8=(z8/(f*(K8-1)+1)) 
yt8=xb8*K8
    
xb9=(z9/(f*(K9-1)+1)) 
yt9=xb9*K9
xb10=(z10/(f*(K10-1)+1)) 
yt10=xb10*K10

print 'Liquid compositions',xb1,xb2,xb3,xb4,xb5,xb6,xb7,xb8,xb9,xb10
print 'Vapor compositions',K1*xb1,K2*xb2,K3*xb3,K4*xb4,K5*xb5,K6*xb6,K7*xb7,K8*xb8,K9*xb9,K10*xb10

    
    

    
    
    
    
    