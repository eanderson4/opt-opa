Must include cplex and concert libraries for the optimization routines, not packaged with this (I am using through Student Intiative by IBM)

Runs example power flow problem

Test cases will be added as I confirm they output same results as matpower, which is great power flow solver in matlab for research

http://www.pserc.cornell.edu/matpower/



To run
===============================

Rename and adjust Makefile.generic to include cplex libraries

./make depend

./make

pow case/30.db



