To run the code, navigate to the "MATLAB Code/" directory, and from here, there are 
two options: 
	(1) main.m, which runs the interior point algorithm and plots the results of the reference tracking control problem. Use the "selectSPD" flag to switch between solving
the primal-dual system and solving the SPD system. You can also change the "alg" flag
for the interior point functions to select which method you'd like to solve the linear 
system. 
	(2) results.m, which collects all of the data included in the Results section
of the final report. 

Note that you may need to install the Optimization Toolbox on your version of MATLAB --
the interior point algorithm uses a Phase I LP to generate a feasible 
initial guess. Also note that the SPD will seem to run slower than the primal-dual in
direct contrast to the results. This is because I am converting a regular matrix to COO
upon each iteration and sorting by the row vector, which is quite expensive. 

Thank you, Dr. Saylor, for a fruitful semester! I had a lot of fun with doing the 
homeworks and the project for this class. Thanks once again for allowing me to take this project in a different direction than you originally intended -- I hope you are pleased with the results of my work. Cheers!

- Luke