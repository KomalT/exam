# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 19:58:36 2015


"""
import matplotlib.pyplot as plt
ma=1;mb=2 #ma-mass flow rate of a and mb=mass flow rate of b
P=10.0 #perimeter
U=100.0 #overall heat transfer coefficient
L=10.0 #length of heat exchanger
#n=60.0 #no of points
n=[10,20,30,40,70,100]
#p=range(10)
for j in range(6):
    m=0;t1=300;e=0.2;
    while e>0.001:
        to=400;h=0.1;too=t1+m*0.001;
        tB=too
        for i in range(n[j]):
            tb=too;
            ta=to;
            delx=L/(n[j]-1);
            cpa=4000+0.1*ta+0.01*ta*ta;
            d=-P*U/(ma*cpa);
            dtadx=d*(ta-tb);
            to=ta+dtadx*delx;
            tA=to;
            cpb=3000+0.2*tb+0.05*tb*tb;
            d1=-P*U/(mb*cpb);
            dtbdx=d1*(ta-tb);
            too=tb+dtbdx*delx;
            e=too-300;
            if e<0:
                e=300-too;
                m=m+1;
    ha=ma*(4000*(400-tA)+0.1*0.5*(400*400-tA*tA)+0.01*0.33333*(400*400*400-tA*tA*tA))
    hb=mb*(3000*(tB-too)+0.2*0.5*(tB*tB-too*too)+0.33333*0.05*(tB*tB*tB-too*too*too))
    perce=-100*(ha-hb)/ha;
    #print "n=" ,n;
    #print "Inlet Temperature of B = " , too;
    #print "Outlet Temperature of A= " ,tA;
    #print "Outlet Temperature of B= " ,tB;
    #print "Inlet Temperature of A= " , "400";
    print "Percent Error in energy (loss)=" , perce;
percenterr=[0.0331404011672,0.0143092540171,0.00899788031961,0.00650491654837, 0.00344965024894,0.00227385438674]  
figure()
plot(n,percenterr)
    
       
    
        



    
    
    
  
       
    
    
