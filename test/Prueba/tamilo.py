The points abs(mu) = 1 can be solved for exactly and are in this subroutine.
Description of eup.f
--------------------
In this subroutine phi(k,l,m,in) is updated to phi(k,l,m,out) for the energy term by Barton's method of monotonic transport.
Description of muup.f
---------------------
In this subroutine phi(k,l,m,in) is updated to phi(k,l,m,out) for the pitch angle diffusion term. The Crank-Nicholson's implicit method is employed and the subroutine cntridag is called to solve the tridiagonal matrix.
Cray J90 routines
-----------------
In order to improve the execution speed of the code the routines below were rewritten to allow multitasking on the Cray J90. To run the code on the Cray J90 please use the following Makefile:
Makefile.J90
and the subroutines ending with “_par”.
A sample NQS batch file for the Cray J90:
#QSUB -s /bin/csh #QSUB -eo
#QSUB -o fkrplk.log setenv NCPUS 20 cd /tmp/fpcode
set echo
set timestamp ja ja.out ./fkrplk.x
ja -cst ja.outprint("Hola mundo")
