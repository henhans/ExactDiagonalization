# ExactDiagonalization
Full exact diagonalization for 1D Hubbard and Anderson model. The code is based on the exercise at 2012 bcgs school's material by Dr. Andrew Mitchell, Dr. Priv.-Doz, Dr. Ralf Bulla, and Prof. Dr. Simon Trebst. http://www.thp.uni-koeln.de/~mitchell/day3/day3.html
The code used lapack and blas or the intel-mkl library.

After compiling, the executable is "ExactDiagonalization.out". The code can gnerate energies, spectralfunction, and green's function at zero temperature. One can use the command line argument to change the parameters.
The commandline arguments are: N, U, mu, t, broaden, model(0 for anderson or 1 for hubbard ).

The next step is to implement the Lanczos algorithm, and finite temperature Green's function. Then, the cluster methods such as CPT.

Tsung-Han Lee
