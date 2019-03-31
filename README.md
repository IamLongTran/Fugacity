# Fugacity
Calculate compressibility factor Z, fugacity f, Enthalpy H, and entropy S of any real fluids (with known Peng-Robinson parameters)
Peng-Robinson paramenters include critical temperature Tc, critical pressure Pc, and accentric factor omega(w). This project was initally developed to calculate butane properties, so omega is already given as butane's. 

This is an on-going project. The end goal of this project is, besides calculating critical properties, determining the equalibrium curve (phase-transition line) on P-T diagram between vapor, liquid, and solid. 
Ultimately, every gas will have its critical properties table that is expected to be as detailed as steam table.

An extensive use of built-in solver optimization would be used to solve for complicated system of nonlinear equations. A good guess point (T_guess, P_guess) might reduce the calculation process time.


Run/debug script file. Terminal/Console will ask for Tc,Pc, omega. 
Run Z(T,P) with T,P are the temperature and pressure at the calculated point. Printed result in the console is the calculated result.

---
Furture release is expcted to have
- calculation in solid phase
- updated functions added for loops calculations (table)
- bugs fixed (if any)
