



import math as m
import numpy as np
from scipy.optimize import root
from numpy import array as ar
from scipy import integrate
'''
name = input('enter name of the fluid:\n')
Tc = float(input('critical temperature in Kelvin:\n'))
Pcr = float(input('critical pressure in bar:\n'))
w = float(input('acentric factor:\n'))
k = float(input('kappa:\n'))
ac = float(input('heat capacity parameter a:\n'))
bc = float(input('b:\n'))
cc = float(input('c:\n'))
dc = float(input('d:\n'))
'''

class R:
    def __init__ (self,e,v):
        R.e = e #Pa
        R.v = v # bar*L/mol/K
R = R(8.314,8.314e-2) 	# J/molK  # L*bar/K/mol


def Z(T,P):
            class Pc:
                def __init__ (self,bar,pa):
                    Pc.bar = bar
                    Pc.pa = pa
            Pc = Pc(Pcr,Pcr*1e5)
                    

            Tr = T/Tc
            #Pr = P/Pc.bar
            
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
            print('ZV is:',ZV)
            if ZL != 'DOE':
                print('ZL is:',ZL)
                    
            
            def fugacity(Z):
                fuga_c = m.e**(Z-1-m.log(Z-B) - A /m.sqrt(8)/B * m.log((Z + (1+m.sqrt(2)) *B)/(Z + (1 - m.sqrt(2)) *B)))
                return (P*fuga_c) # fugacity
            print('fugacity V is:',fugacity(ZV))
            fV = fugacity(ZV)
            if ZL != 'DOE':
                fL = fugacity(ZL)
                print('fugacity L is:',fugacity(ZL))
                if abs(fL-fV) <= 1e-4:
                    print(name + ' is in SAT CURVE')
            else:
                fL = ZL
                if P > Pc.bar:
                    print(name + ' is in LIQUID PHASE')
                else:
                    print(name + ' is in VAPOR PHASE')
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

if __name__ == '__main__':
    name = input('enter name of the fluid:\n')
    Tc = float(input('critical temperature in Kelvin:\n'))
    Pcr = float(input('critical pressure in bar:\n'))
    w = float(input('acentric factor:\n'))
    k = float(input('kappa:\n'))
    ac = float(input('heat capacity parameter a:\n'))
    bc = float(input('b:\n'))
    cc = float(input('c:\n'))
    dc = float(input('d:\n'))
    T = float(input('T:\n'))
    P = float(input('P:\n'))
    print('')
    Z(T,P)
    
                
                
                
                