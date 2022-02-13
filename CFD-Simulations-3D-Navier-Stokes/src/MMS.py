# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 12:34:05 2015

@author: Mo
"""
from sympy import *
x,y,z = symbols('x y z')

X = Matrix([[x],[y],[z]])

psi = zeros(3,1)
psi[0,0] = sin(2*pi*x)*y**2*(1-y)**2*z**2*(1-z)**2
psi[2,0] = x**2*(1-x)**2*y**2*(1-y)**2*sin(2*pi*z)

curl_psi = zeros(3,1)
ux = curl_psi[0] = diff(psi[2],X[1]) - diff(psi[1],X[2])
uy = curl_psi[1] = diff(psi[0],X[2]) - diff(psi[2],X[0])
uw = curl_psi[2] = diff(psi[1],X[0]) - diff(psi[0],X[1])