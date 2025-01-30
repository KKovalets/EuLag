The Eulag is an Eulerian-Lagrangian model to solve equations describing gravitational sinking of organic particles under influence of the size and age of particles, temperature and oxygen concentration on their dynamics and degradation processes. 
EuLag.f90 is the source file. It was tested under Windows 11 using Intel Fortran Compiler with VS2008 shell.

Input files: 
Input.nml - here parameters are defined, 
Temperature.dat - temperature profile (C), 
O2_concentration.dat - oxygen concentration profile (kg/m^3).

Output files: 
if time-independent degradation rate is considered,
Fd(z)_0.dat - total POM flux profile (kg/m^2/s), 
Sp(z)_0.dat - total POM concentration profile (kg/m^3),  
Fd(z)_norm_0.dat, Sp(z)_norm_0.dat - normalized total POM flux and concentration profiles, 
gamma_0.dat - degradation rate profile (1/s), 
M0xCm.dat - the product of M0 and Cm parameters, calculated in program using measurements data. 
If time-dependent degradation rate is considered, then alpha-value in days is added in the end of the filenames:
for example, Fd(z)_30.dat, Sp(z)_30.dat, Fd(z)_norm_30.dat, Sp(z)_norm_30.dat, gamma_30.dat if alpha=30 days)

How to use:
Step 1: input the parameters into input-file. Follow the instructions to input the data correctly.
Step 2: If you want to change temperature and oxygen-concentration profiles, replace the corresponding data-files and parameters n_data1 and n_data2 in the input-file.
Step 3: Compile an executable file using EuLag.f90 and run it.
