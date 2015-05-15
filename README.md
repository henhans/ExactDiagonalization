# ExactDiagonalization
Full exact diagonalization for 1D Hubbard and Anderson model with open boundary condition(OBC). The code is based on the exercise at 2012 bcgs school's material by Dr. Andrew Mitchell, Dr. Priv.-Doz, Dr. Ralf Bulla, and Prof. Dr. Simon Trebst. http://www.thp.uni-koeln.de/~mitchell/day3/day3.html
The code used lapack and blas or the intel-mkl library. Using Make to compile the program.

After compiling, the executable is "ExactDiagonalization.out". The code can gnerate energies, spectralfunction(at 1st site), and green's function at zero temperature(at first site). One can use the command line argument to change the parameters. The commandline arguments are: N, U, mu, t, broaden, model(0 for anderson or 1 for hubbard ). The code use the symmetry of total number and total Sz.

Some implementation notes and references can be found in directory doc.

The next step is to implement the Lanczos algorithm, finite temperature Green's function, and other observables. Then, the cluster methods such as CPT.

Tsung-Han Lee
