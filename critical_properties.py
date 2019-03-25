# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:36:51 2019

@author: Long Tran
"""





import math as m
import numpy as np
from scipy.optimize import root
from numpy import array as ar
from scipy import integrate

class Pc:
    def __init__(self, pa,bar):
        Pc.bar = bar
        Pc.pa = pa
        
class R:
    def __init__ (self,e,v):
        R.e = e #Pa
        R.v = v # bar*L/mol/K
R = R(8.314,8.314e-2) 	# J/molK  # L*bar/K/mol
        
print('Please provide physical properties of the fluids below:')

    #critical prop.
Tc = 425.2	 #K
#Pc = Pc(3.8e6,38)	 	#bar
Vc = 0.255 		#m3/kmol ~ L/mol 
Zc = 0.274 
w = 0.193 # omega
k = 0.37464 + (1.5422 * w) - (0.26992 * (w**2))
Tc = float(input('Please enter the critical temperature in Kelvin:\n'))
Pc_input = float(input('Please enter the critical pressure in bar:\n'))
T = float(input('Enter the temperature:\n'))
P =float(input('Enter the pressure:\n'))

Pc = Pc(Pc_input*1e5, Pc_input)
'''
Tbp = 272.7	#K at P = 1 bar ~ atm
R = R(8.314,8.314e-2) 	# J/molK  # L*bar/K/mol
Pt = .7e-5 # bar
Tt = 134.6 #K
'''
    ## Cp*
    # Cp = a + bT + cT^2 + dT^3
ac = 3.954
bc= 37.12e-2
cc = -18.326e-5
dc= 34.979 *10 **-9
    ## Antoine Eqn
    #	lg(P) = A - B/(T+C)	[P] = [mmHg], 	 [T] = [‘C]
A = 6.82485
B = 943.453
C = 239.711
def Z(T,P):
            Tr = T/Tc
            Pr = P/Pc.bar
            
            a = 0.457235 * R.v**2 * Tc**2 / Pc.bar ## not including alpha
            b = 0.0777961 * R.v * Tc / Pc.bar
            
            alpha = (1 + k* (1 - np.sqrt(Tr)))**2
            A = a * alpha * P / R.v**2 / T**2
            B = b * P / R.v / T
                        
            coeff = [1,-1+B,A-3*B**2-2*B,-A*B+B**2+B**3]
            sol = np.roots(coeff)
            sol = np.real(sol[np.imag(sol)==0])
            ZV = sol[0]
            ZL = np.real(sol[-1]) if np.imag( sol[-1]) == 0 and sol[0] != sol[-1] else 'DOE'
            print('Z is:',ZV)
            
            
            def fugacity(Z):
                fuga_c = m.e**(Z-1-m.log(Z-B) - A /m.sqrt(8)/B * m.log((Z + (1+m.sqrt(2)) *B)/(Z + (1 - m.sqrt(2)) *B)))
                return (P*fuga_c) # fugacity
            print('fugacity is:',fugacity(ZV))
            fV = fugacity(ZV)
            if ZL != 'DOE':
                fL = fugacity(ZL)
            else:
                fL = ZL
            Ti = 298 #K
            Pi = 1 #bar
            derv_aPR = -0.45724 * R.e**2*Tc**2*k*m.sqrt(alpha/T/Tc)/Pc.pa
            aPR = 0.45724 * R.e**2*Tc**2/Pc.pa * alpha
            bPR = 0.0778 * R.e * Tc/Pc.pa
            BPR = bPR * P*1e5 / R.e / T
            
            
            def S(T,P):
                
                #alpha = (1 + k*(1-m.sqrt(Tf/Tc)))**2
                #derv_aPR= -0.45724 * R.e**2*Tc**2*k*m.sqrt(alpha/Tf/Tc)/Pc.pa
                #B = bPR * Pf /(R.e*Tf)
               # Zi = crit_prop(butane,Ti,Pi)[0] * 1e5
                fS = lambda x: (ac + bc*x + cc*x**2 + dc*x**3)/x
                S = integrate.quad(fS,Ti,T)[0] - R.e*m.log(P/Pi)
                SV = S + (R.e)*m.log(ZV- BPR) + derv_aPR*m.log((ZV + (1+m.sqrt(2))*BPR)/(ZV+(1-m.sqrt(2))*BPR))/m.sqrt(8)/bPR
                if ZL != 'DOE':
                    SL = S + (R.e)*m.log(ZL- BPR) + derv_aPR*m.log((ZL + (1+m.sqrt(2))*BPR)/(ZL+(1-m.sqrt(2))*BPR))/m.sqrt(8)/bPR
                    print('S liquid @ T = %.2f and P = %.2f is'%(T,P), SL, 'J/mol')
               
                else:   
                    SL = ZL
                print('S vapor @ T = %.2f and P = %.2f is'%(T,P), SV, 'J/mol')
                return SL,SV
            Sliq,Svap = S(T,P)
           
            def H(T,P):
                fH = lambda x:  ac + bc*x + cc*x**2 + dc*x**3
                H = integrate.quad(fH,Ti,T)[0]
                HV = H + R.e * T *(ZV - 1) + (T * derv_aPR - aPR) * m.log((ZV + (1+m.sqrt(2))*BPR)/(ZV+(1-m.sqrt(2))*BPR))/m.sqrt(8)/bPR
                if ZL != 'DOE':
                    HL = H + R.e * T *(ZL - 1) + (T * derv_aPR - aPR) * m.log((ZL + (1+m.sqrt(2))*BPR)/(ZL+(1-m.sqrt(2))*BPR))/m.sqrt(8)/bPR
                    print('H liquid @ T = %.2f K and P = %.2f bar is' %(T,P), HL, 'J/mol')
                
                
                else:
                        HL = ZL
                print('H vapor @ T = %.2f K and P = %.2f bar is' %(T,P), HV, 'J/mol')
                return HL,HV
            Hliq,Hvap = H(T,P)
            
            return ZL,ZV,fL,fV,Sliq,Svap,Hliq,Hvap
Zliq,Zvap,fliq,fvap,Sliq,Svap,Hliq,Hvap = Z(T,P)
'''
            def equations(x):
                #(Z,P [bar],V [L/mol]) = x
                a = 0.457235 * R.v**2 *Tc**2 / Pc.bar ## not including alpha
                b = 0.0777961 * R.v *Tc /Pc.bar
                  
               
                A = a * alpha * x[1] / R.v**2 / T**2
                B = b * x[1] / R.v / T
                eqn1 = x[0]**3 + (-1+B) * x[0]**2 + (A-3*B**2-2*B)*x[0] +(-A*B+B**2+B**3) # cubic Z
                eqn2 = x[1]*x[2] - x[0]*R.v*T # PV = ZRT
                eqn3 = -R.v*T/(x[2]-b)**2 + 2*a/(x[2])**3  #dP/dV = 0
                return ar([eqn1,eqn2,eqn3])
            sol1 = root(equations,ar([1,1,1]),method = 'broyden2')     
            print('')
            print(sol1)
            sol2 = root(equations, ar([0.05,0.05,0.05]), method ='broyden2')
            print('')
            print(sol2)
            
            def d2P_dV2():
                return
            
                
'''
                
                
                
                
              