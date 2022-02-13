# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 15:24:30 2015

@author: Mo
"""

#f1 = open("12e5Re220dt300.txt").readlines()
#f2 = open("1e5Re220dt300dpdx18.txt").readlines()
#f1 = open("1e5Re220T1000dpdx18.txt").readlines()
f1 = open("test.pckl").readlines()
#f2 = open("1e5Re220T1000dpdx18.txt").readlines()
#N=len(f1)
t = range(len(f1))
R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; R7=[]; R8=[]; R9=[]; R10=[]
for l in f1:
    #l = float(l)
    #print l1
    R1.append(l)


plt.plot(t, R1, 'magenta')

plt.xlabel('t')
plt.ylabel('R=R3/R2')
p1 = plt.Rectangle((0, 0), 1, 1, fc='magenta')
#p2 = plt.Rectangle((0, 0), 1, 1, fc='g')
plt.legend([p1], ['C=1.2e-5'])
plt.title('R')

plt.hold('on')