      subroutine bfield(s,b,dbds)
c
c  Returns the logarithmic derivative of the magnetic field 
c  dbds=dlnB/ds and magnetic field strength b 
c  at a given spatial position s.
c
c  The magnetic field used here has the following properties:
c  B = B_o(1 + s^2/L_B^2).  L_B is determined by the length of the
c  corona L_C and the mirror ratio r_m.  L_B = L_C / sqrt(r_m - 1).
c
c  xc = length of corona 
c  rm = mirror ratio
c
      parameter(b0=100.0)
      parameter(xc = 1.0e9 , rm = 2.0)
      b = b0*(1. + (rm-1.)*(s**2)/xc**2 )
      dbds = 2.*s*(rm-1.) / ( xc**2 + (rm-1.)*s**2 )
      return 
      end