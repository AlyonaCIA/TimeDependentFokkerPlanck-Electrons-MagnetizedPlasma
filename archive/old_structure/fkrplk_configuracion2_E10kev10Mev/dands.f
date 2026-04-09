      subroutine dands(t,rho,s)
c
c  Returns the electron number density (rho)
c  and spatial coordinate (s) given the "column depth" (t).
c
c  The density function used here has the following properties:
c     the density is a constant (rho0).
c
      parameter(xlmbda=20.0, a=1.0e24, rho0=1.0e11)
c
c  Convert "column depth" (t) to distance s, then compute 
c  the density.
c
      c=xlmbda*rho0/a
      rho = rho0
      s = t/c
      return 
      end