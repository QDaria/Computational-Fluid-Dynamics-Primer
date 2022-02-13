# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 18:10:31 2015

@author: Mo
"""

f1 = open("fbi842D1e2dt100ts500.txt").readlines()
f2 = open('fbi842D1e3dt10ts5000.txt').readlines()
f3 = open('fbi842D1e4dt100ts500.txt').readlines()
#f1 = open("Rf2e4an3DRe2dt10ts500europa5020Oct29.txt").readlines()
#f2 = open('Rf15e4an3DRe2dt10ts500europa5020Oct29.txt').readlines()
#f3 = open('Rf1e4an3DRe2dt10ts500europa5020Oct29.txt').readlines()
#f2 = open("1e5Re220T1000dpdx18.txt").readlines()

#f1 = open("Rf2e4an3DRe2dt100ts300europa5020Oct29.txt").readlines()
#f2 = open('Rf15e4an3DRe2dt100ts300europa5020Oct29.txt').readlines()
#f3 = open('Rf1e4an3DRe2dt100ts300europa5020Oct29.txt').readlines()
#N=len(f1)
t = range(len(f1))
R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; R7=[]; R8=[]; R9=[]; R10=[]
for l in f1:
    l = float(l)
    #print l1
    R1.append(l)

for l in f2:
    l = float(l)
    #print l1
    R2.append(l)
    
for l in f3:
    l = float(l)
    #print l1
    R3.append(l)

plt.plot(t[0:50], R1[0:50], 'magenta')
plt.plot(t[0:50], R2[0:50], '-g')
plt.plot(t[0:50], R3[0:50], '-b')
#plt.plot(t[0:500], R1[0:500], 'magenta')
#plt.plot(t[0:500], R2[0:500], '-g')
##plt.plot(t[0:500], R3[0:500], '-b')
#plt.plot(t[0:250], R1[0:250], 'magenta')
#plt.plot(t[0:250], R2[0:250], '-g')
#plt.plot(t[0:250], R3[0:250], '-b')
#plt.plot(t[0:100], R1[0:100], 'magenta')
#plt.plot(t[0:100], R2[0:100], '-g')
#plt.plot(t[0:100], R3[0:100], '-b')

plt.xlabel('time step')
plt.ylabel('R = R3/R2')
p1 = plt.Rectangle((0, 0), 1, 1, fc='magenta')
p2 = plt.Rectangle((0, 0), 1, 1, fc='g')
p3 = plt.Rectangle((0, 0), 1, 1, fc='b')
plt.legend([p1, p2, p3], ['C=1e-2','C=1e-3','C=1e-4'])
plt.title('Re = 84: dt=100')
grid(True)
plt.hold('on')

