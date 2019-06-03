# Fugacity

Calculate compressibility factor Z, fugacity f, Enthalpy H, and entropy S of any real fluids (with known Peng-Robinson parameters)

Peng-Robinson paramenters include critical temperature (Tc), critical pressure (Pc), and accentric factor omega(w). 

This is an ongoing project. The end goal of this project is determining the equalibrium curve (phase-transition line on P-T diagram) between vapor phase and liquid phase, as well as critical thermodynamics properties of at any specific P and T.
Ultimately, any chemical will have its own "steam table."

An intensive use of solver optimization is used to solve complicated system of nonlinear equations. A good guess point (T_guess, P_guess) could reduce the calculation process time.

---
Run/debug script file. Terminal/Console will ask for Tc,Pc, omega, etc. 
Run Z(T,P) with T,P are the temperature and pressure at any given condition. Printed output is the calculated result.

---
Furture release is expcted to have
- functions added for loops calculations (property table) (matrix of data sets)
- bugs fixed (if any)
