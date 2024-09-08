
      program main
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to simulate the evolution of the electron velocity
c     distribution in a magnetized plasma using the Fokker-Planck
c     kinetic equation.
c
c     Reference: Hamilton, Lu, and Petrosian, ApJ 354, 726 (1990)
c
c     Last modification: 2000/05/12 by J. McTiernan
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
c     Declare variables
c     ...
      real*8 tr(10) ! Times at which to report progress
      integer*4 ntr ! Number of report times
c     ...
c     Set initial parameters and values
c     ...
c     Initialize output file for simulation results
      open(7, file='fkrplk.test', status='unknown')
c     Call subroutine to initialize arrays and set initial conditions
      call init
c     Write initial parameters to output file
c     ...
c     Time-stepping loop to advance the simulation
      do 1 i=1,nt
c       Logic to handle reporting times and adjust time step
c       ...
c       Call subroutines to manage electron injection and update
c       the distribution function in space
        call inject
        call tauup
c       Logic to handle data writing and progression reporting
c       ...
c     End of time-stepping loop
    1 continue
      end
