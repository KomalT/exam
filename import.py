# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 19:29:15 2016

@author: BCS
"""

import scipy
from scipy import array
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
import pandas as pd

location='C:\python68\sharya.xlsx'
df=pd.read_excel(location,0)
pdata=df['P']
P1=np.array([pdata])[0]
print P1[0]
P=P1[0]