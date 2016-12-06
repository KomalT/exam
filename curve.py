# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 20:07:59 2016

@author: S
"""

import matplotlib.pyplot as plt
from scipy import array
from scipy import linspace


R=0.5
yr=R/(R+1)*x+0.9/(R+1)
print yr
def curve2(x,C):
    
    y=C*x/(1+(C-1)*x)
    return y
x=linspace(0,1,10)
y=linspace(0,1,10)
ycalc=curve2(x,2)
fig=plt.figure();
ax=fig.add_subplot(110)
ax.plot(x,ycalc,'b');
ax.plot(x,y,'b');
ax.plot([0.5,0.6],[0.5,0.7],'b');
fig.canvas.draw()
plt.show()
for x in range(0.5,1):
    yf
    xf=yf/((2-1)*yf+1)
    
