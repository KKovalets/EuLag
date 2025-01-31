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
Step 1: input the parameters into input-file. Follow the instructions to input the data correctly. In this release parameters and temperature and oxygen profiles are set as ones that were used in our simulations for the Atlantic Ocean.

Step 2: If you want to model some specific ocean and change temperature and oxygen concentration profiles, replace the corresponding data-files "Temperature.dat" and "O2 concentration.dat" with new ones with the same names. You can found measurements data in https://github.com/KKovalets/EuLag_DataSet repository or in EuLag_DataSet release. Change the corresponding parameters n_data1 and n_data2 in the input-file. n_data1 is a number of lines in "Temperature.dat" and n_data2 is a number of lines in "O2 concentration.dat" file.

Step 3: If you are modeling some specific ocean and you have set new already temperature and oxygen concentration profiles, set new Z_measur and Sp_measur or Fd_measur values in the Input.nml file. To find the corresponding POM concentration (Sp) or flux (Fd) measurments, you can also check https://github.com/KKovalets/EuLag_DataSet repository or an EuLag_DataSet release. 
In our simulations we used POM concentration measurments. To reproduce them for some specific ocean open the corresponding file from EuLag_Data repository or release with POM concentrations profiles in this ocean and take the highest measurment from there. Then set the Z_measur as rounded to the nearest integer Z-coordinate and Sp_measur as corresponding POM concentration value.
In general, user can also use POM flux measurments as a parameter instead of concentration. To do this, set the Sp_measur=0. and assign the Z_measur and Fd_measur parameters the corresponding values of the highest measurement from corresponding POM flux file. We didn 't do it in our simulations.

Step 4: Compile an executable file using EuLag.f90 and run it.
