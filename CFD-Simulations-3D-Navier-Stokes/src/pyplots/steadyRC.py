# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 19:58:55 2015

@author: Mo
"""
#import matplotlib as plt
#import matplotlib as plt
#from numpy import *
#from scitools.std import * 
#plt.5
#import matplotlib as plt
#C = [1e-2, 7.5e-3, 5e-3, 2.25e-3, 1e-3, 7.5e-4, 5e-4, 2.25e-4, 1e-4]
#C = linspace(0.0005,0.0100, 9) 
#C = [0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005]
#R2 = [1.9, 2.2, 2.7, 4.4, 8.9, 11.2, 15.2, 26.3, 40]
#R84= [2.1, 2.4, 3.0, 4.7, 7.3, 8.5, 10.1, 13.2, 15.5]
#C = seq(1e-4, 1e-2, 0.0025)

xfmt = ScalarFormatter()
C = [1e-4, 2.25e-4, 5e-4, 7.5e-4, 1e-3, 2.25e-3, 5e-3, 7.5e-3, 1e-2] 
R2 = [40, 26.3, 15.2, 11.2, 8.9, 4.4, 2.7, 2.2, 1.9]
R84= [15.5, 13.2, 10.1,8.5 ,7.3 , 4.7, 3.0, 2.4,2.1 ]
p1 = Rectangle((0, 0), 1, 1, fc='r')
p2 = Rectangle((0, 0), 1, 1, fc='b')
#plot(C, R2,  '-ro')
#plot(C, R84, '-bo')

legend([p1, p2], ['Re=2.2', 'Re=84'])
xlabel('dC = 0.25')
ylabel('R')
axhspan(2.1, 15.5, facecolor='0.7', alpha=0.2)
hold('on')
title('The impact of Re on R=R3/R2 at steady-state')
grid(True)



Ra = pearsonr(R2,R84)
Rb = pearsonr(R2[1:],R84[1:])
Rc = pearsonr(R2[2:],R84[2:])
Rd = pearsonr(R2[3:],R84[3:])
Re = pearsonr(R2[4:],R84[4:])
Rf = pearsonr(R2[5:],R84[5:])
Rg = pearsonr(R2[6:],R84[6:])
Rh = pearsonr(R2[7:],R84[7:])
Ri = pearsonr(R2[8:],R84[8:])

Rcor = [Ra[0],Rb[0],Rc[0],Rd[0],Re[0],Rf[0],Rg[0],Rh[0],Ri[0]]#0.99999999999999999
Rstd = [Ra[1],Rb[1],Rc[1],Rd[1],Re[1],Rf[1],Rg[1],Rh[1],Ri[1]]
plot(C, Rcor,  '-ro')
#plot(C, Rstd, '-bo')
plt.axis([1e-4,1e-2,0.95,1.1])
"""
import numpy as np
import scipy.stats as stats
import pylab as pl
h = sorted([186, 176, 158, 180, 186, 168, 168, 164, 178, 170, 189, 195, 172,
     187, 180, 186, 185, 168, 179, 178, 183, 179, 170, 175, 186, 159,
     161, 178, 175, 185, 175, 162, 173, 172, 177, 175, 172, 177, 180])  #sorted
fit = stats.norm.pdf(h, np.mean(h), np.std(h))  #this is a fitting indeed

pl.plot(h,fit,'-o')

pl.hist(h,normed=True)      #use this to draw histogram of your data

pl.show() 

R2a = [40, 26.3, 15.2, 11.2, 8.9, 4.4, 2.7, 2.2, 1.9]
fit = stats.norm.pdf(R2, np.mean(R2), np.std(R2))  #this is a fitting indeed

pl.plot(R2a,fit,'-o')

pl.hist(R2a,normed=True)      #use this to draw histogram of your data

pl.show() 

"""
                  #use may also need add this 
#a = pearsonr(a1,a2)
"""
x = linspace(0,2*pi,20);
y = sin(x); yp = None
xi = linspace(x[0],x[-1],100);
yi = stineman_interp(xi,x,y,yp);

fig, ax = plt.subplots()
ax.plot(x,y,'ro',xi,yi,'-b.')
plt.show()

#plt.axhline(y=.5, xmin=0.25, xmax=0.75)
#plt.axhspan(5.45, 5.55, facecolor='0.7', alpha=0.7)

#axhline(y=.5, xmin=0.25, xmax=0.75)
#plt.axhspan(5., 5.5, facecolor='0.7', alpha=0.7)
#p = plt.axvspan(1.25, 1.55, facecolor='g', alpha=0.5)
#
#plt.axis([0,1000,4,6])
#plt.axis(0, 100)
"""