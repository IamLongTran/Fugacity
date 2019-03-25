# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:28:22 2019

@author: Long Tran
"""
import math as m
import numpy as np
from scipy.optimize import root
from numpy import array as ar

class Pc:
    def __init__(self, pa,bar):
        Pc.bar = bar
        Pc.pa = pa
        
class R:
    def __init__ (self,e,v):
        R.e = e #Pa
        R.v = v # bar*L/mol/K
        

    #critical prop.
Tc = 425.2	 #K
Pc = Pc(3.8e6,38)	 	#bar
Vc = 0.255 		#m3/kmol ~ L/mol 
Zc = 0.274 
w = 0.193 # omega
k = 0.37464 + (1.5422 * w) - (0.26992 * (w**2))

Tbp = 272.7	#K at P = 1 bar ~ atm
R = R(8.314,8.314e-2) 	# J/molK  # L*bar/K/mol
Pt = .7e-5 # bar
Tt = 134.6 #K

    ## Cp*
    # Cp = a + bT + cT^2 + dT^3
a = 3.954
b = 37.12e-2
c = -18.326e-5
d = 34.979 *10 **-9
    ## Antoine Eqn
    #	lg(P) = A - B/(T+C)	[P] = [mmHg], 	 [T] = [‘C]
A = 6.82485
B = 943.453
C = 239.711
class Gas:
    def __init__ (self, Tc, Pc, w):
        
        self.Tc = Tc
        self.Pc = Pc
        self.w = w

butane = Gas(Tc,Pc.bar,w=w)

def Z(butane, T,P):
            Tr = T/butane.Tc
            Pr = P/Pc.bar
            
            a = 0.457235 * R.v**2 * butane.Tc**2 / butane.Pc ## not including alpha
            b = 0.0777961 * R.v * butane.Tc / butane.Pc
            
            alpha = (1 + k* (1 - np.sqrt(Tr)))**2
            A = a * alpha * P / R.v**2 / T**2
            B = b * P / R.v / T
                        
            coeff = [1,-1+B,A-3*B**2-2*B,-A*B+B**2+B**3]
            sol = np.roots(coeff)
            sol = np.real(sol[np.imag(sol)==0])
            ZV = sol[0]
            ZL = np.real(sol[-1]) if np.imag( sol[-1]) == 0 and sol[0] != sol[-1] else 'DOE'
            
            
            def fugacity(Z):
                fuga_c = m.e**(Z-1-m.log(Z-B) - A /m.sqrt(8)/B * m.log((Z + (1+m.sqrt(2)) *B)/(Z + (1 - m.sqrt(2)) *B)))
                return (P*fuga_c) # fugacity
            fV = fugacity(ZV)
            if ZL != 'DOE':
                fL = fugacity(ZL)
            else:
                fL = 'DOE'
                
            
            ## The VLE occurs when Tt < T < Tc.
                # the interval of P @ VLE depends on the chemical
            
            ## Solving for P range 
                # dP/dV = 0
                    # d2P/dV2 > 0 <=> P of liquid @ VLE
                    
                    # d2P/dV2 < 0 <=> P of vapor @ VLE
            ## dP/dV = -R*Tc/(Vc - b)**2  + 2*a/Vc**3
                # reduce dP/dV for simplication
                # Tc**2/Pc = 1 and Tc/Pc = 1 (in a and b)
                
            ## d2P/dV2 = 2*R*Tc/(Vc-b)**3  - 6*a/Vc**4
            
            ## Pc * Vc = Z * R * Tc ## EOS of real gas 
            # To solve for the 2 peaks 
            # => solve the system of 3 equations for 3 variables(P,Z,V) @ fixed T
                # Z = f(P,T)
                # PV = ZRT
                # dP/dV = 0
            if fL != 'DOE':
                delta = abs(fL-fV)
                print('the difference between fL and fV is', abs(fL-fV))
    
            #return ZV,ZL,fV,fL,sol,delta
        
'''
def equations(x,maxiter = 100000):
    #(Z,P,V) = x
    a = 0.457235 * R.v**2 *Tc**2 / Pc.bar ## not including alpha
    b = 0.0777961 * R.v *Tc /Pc.bar
      
    alpha = (1 + k* (1 - np.sqrt(T/Tc)))**2
    A = a * alpha * x[1] / R.v**2 / T**2
    B = b * x[1] / R.v / T
    eqn1 = x[0]**3 + (-1+B) * x[0]**2 + (A-3*B**2-2*B)*x[0] +(-A*B+B**2+B**3)
    eqn2 = x[1]*x[2] - x[0]*R.v*T
    eqn3 = -R.v*T/(x[2]-b)**2 + 2*a/(x[2])**3
    return ar([eqn1,eqn2,eqn3])
sol = root(equations,ar([1.5,5,5]),method = 'broyden2')

'''

'''
pressure = []
T25 = np.arange (150,400,10)
Pmin = 1.5e-3
Pmax = 2*Pmin
for T in T25:
    (ZV,ZL,fV,fL,sol,delta) = Z(butane,T,Pmax)
    di = delta
    P = Pmax
    while delta >=5e-5:
        
        P = (Pmax+Pmin)/2
        (ZV,ZL,fV,fL,sol,delta) = Z(butane,T,P)
        df = delta
        if df > di:
            Pmax = P
            
            P = np.average(Pmax,Pmin)
            (ZV,ZL,fV,fL,sol,delta) = Z(butane,T,P)
            
        elif df<di and delta >= 5.5e-5:
            Pmin = P
            
            P = np.average(Pmax,Pmin)
            (ZV,ZL,fV,fL,sol,delta) = Z(butane,T,P)
        else:
            
            break
    Pmax = 2*P
    P.append(pressure)

'''

        
        
'''
Pmax = Pc
Pmin = 1
for T in T25:
    P = (Pmax + Pmin)/2
    if Z(butane,T,P)[1] == 'DOE':
        P = Pmax
'''


'''
P = int(input('Guess P:\n'))
(ZV,ZL,fV,fL,sol)=Z(butane,T,P)
i = 0
while fL == 'DOE' and i in range(10):
    P /= 2
    (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
while fL == 'DOE' and i in range(10):
    P *= 2
    (ZV,ZL,fV,fL,sol)=Z(butane,T,P)


while abs(fL-fV) >= 1e-4:
    P = P * fL/fV
    (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
    
'''

'''
while fL =='DOE' or abs(fV-fL) >= 1e-4:
        (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
        i = 0
        while i in range(10) and ZL =='DOE':
            P = (Pmax+P)/2
            (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
            i+=1
        if ZL != 'DOE':
            while fV>fL:
                Pmin = P
                P = (Pmax+P)/2
                (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
            while fV<fL:
                Pmax = P
                P = (Pmin+P)/2
                (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
                
        i=0
        while i in range(10) and ZL== 'DOE':
            P = (Pmin+P)/2
            (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
            i+=1
        if ZL != 'DOE':
            while fV>fL:
                Pmin = P
                P = (Pmax+P)/2
                (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
            while fV<fL:
                Pmax = P
                P = (Pmin+P)/2
                (ZV,ZL,fV,fL,sol)=Z(butane,T,P)
            
            
        
'''
    




