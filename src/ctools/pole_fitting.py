#! /usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit

def pole_func(x,*args):
   y = args[0]
   for i in range(1,len(args),3):
      pole_re=args[i]
      pole_im=args[i+1]
      res=args[i+2]
      y = y + res*pole_im/((x-pole_re)**2 +pole_im**2)
   return y

def pole_fitting(xdata,ydata,npol,p0=None,method='dogbox'):

   if (p0==None):
      p0=[1.0 for x in range(3*npol+1)]
   #print('p0 ',p0)
   poles,pcov = curve_fit(pole_func,xdata,ydata,p0=p0,method=method)
   res=list(poles)
   res.append(pcov.trace())
   #print(res)
   return res

# main
#print('Entering Python exec')
params=pole_fitting(xdata,ydata,npol,p0)
#print('Closing Python exec')
#print();
   
