---
Calculate compressibility factor Z, fugacity f, Enthalpy H, and entropy S of any real fluids (with known Peng-Robinson parameters)

Peng-Robinson paramenters include critical temperature (Tc), critical pressure (Pc), and acentric factor omega(w). 

This is an ongoing project. The end goal of this project is to determine the equalibrium curve (phase-transition line on P-T diagram) between vapor phase and liquid phase, as well as critical thermodynamics properties of any fluids, at any given P and T.
Ultimately, any real fluid will have its own "steam table."

An intensive use of solver optimization is used to solve complicated system of nonlinear equations. A good guess point (T_guess, P_guess) could reduce the calculation process time.

---
Run script file. Terminal/Console will ask for Tc,Pc, omega, etc. 
Run Z(T,P) with T,P are temperature and pressure at the condition of interest. Printed output is the calculated result.

---
Furture release is expcted to have
- replicated (loops) calculations for full table of properties (matrix of data sets).
- ...
