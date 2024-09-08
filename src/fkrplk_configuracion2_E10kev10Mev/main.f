c23456789012345678901234567890123456789012345678901234567890123456789012
c fkrplk is a new version of the Russ H, Ed L, cpuhog program, at      c
c present (7-feb-1997), the only difference between this and the       c
c original version is that temporary arrays are used for Phi in the    c
c *up routines, this speeds things up by about 20%.                    c
c jmm, 7-feb-1997                                                      c
c Reference: Hamilton, Lu and Petrosian, Apj 354, 726 (1990)           c
c                                                                      c
c Changes by Leon Ofman:                                               c
c 3-1-2001 /univ/ and /univ1/ were moved to fkrplk.h                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fkrplk
c
c  This is it!!
c
      include 'fkrplk.h'
      common /ctau/ dtau(ntaumax), rlmbda2(ntaumax)
      common /bmag/ bm(0:ntaumax)
      common /ce/ de(nemax), beta2(0:nemax), e(-1:nemax+1)
      common /cmu/ dmu(-nmumax+1:nmumax), onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax), rb3g2(0:nemax), bconv(0:ntaumax)


      real soft(0:ntaumax), tau(0:ntaumax)
      real tr(0:99)
      logical symloop, impulsive

      data tr /0.00,0.11,0.21,  0.32,  0.42,  0.53,  0.64,  0.74,  0.85,
     +     0.96,  1.06,  1.17,  1.27,  1.38,  1.49,  1.59,  1.70,  1.80,
     +     1.91,  2.02,  2.12,  2.23,  2.33,  2.44,  2.55,  2.65,  2.76,
     +     2.87,  2.97,  3.08,  3.18,  3.29,  3.40,  3.50,  3.61,  3.71,
     +     3.82,  3.93,  4.03,  4.14,  4.25,  4.35,  4.46,  4.56,  4.67,
     +     4.78,  4.88,  4.99,  5.09,  5.20,  5.31,  5.41,  5.52,  5.62,
     +     5.73,  5.84,  5.94,  6.05,  6.16,  6.26,  6.37,  6.47,  6.58,
     +     6.69,  6.79,  6.90,  7.00,  7.11,  7.22,  7.32,  7.43,  7.54,
     +     7.64,  7.75,  7.85,  7.96,  8.07,  8.17,  8.28,  8.38,  8.49,
     +     8.60,  8.70,  8.81,  8.92,  9.02,  9.13,  9.23,  9.34,  9.45,
     +     9.55,  9.66,  9.76,  9.87,  9.98, 10.08, 10.19, 10.29, 10.40,
     +     10.51/    
c
c  Number of times to report (ntr) and the maximum time of 
c  the simulation tmax.
c
      ntr = 99
      tmax = tr(ntr)
      itr = 0
c
c  Number of spatial grid points (ntauwr) and energy grid
c  points (newr) to write to output file.
c
      ntauwr = 1
      newr = 120 
c
c  logical variables
c
      symloop = .true.
      impulsive = .false.
c
c  Output file 
c
      open(unit=10,file='fkrplk.conf2_E10kev10Mev')
c
c  Intialize the arrays ...
c
      call init(rhomax,soft,tau)
c
c  Convert dimensionless time "y" to real time "t".
c  Assuming that ln(Lambda) = 20 and fully ionized hydrogen.
c
      ytot = 1.6666667e12/rhomax
      ttoy = 1./ytot
      dt = ytot*dy
      t = 0.0
      i = 0
c
 202  format(1p6e13.5)
 204  format(15i5)
      write(10,202) dt, tmax, dy, rhomax
      write(10,204) newr,nmu,ntau/ntauwr,ntr
      write(10,202) (511.*e(k), k=0,newr)
      write(10,202) (amu(l), l=-nmu,nmu)
      write(10,202) (dtau(m), m=ntauwr,ntau,ntauwr)
      write(10,202) (tau(m), m=ntauwr,ntau,ntauwr)
      write(10,202) (soft(m), m=ntauwr,ntau,ntauwr)
      write(10,202) (bm(m), m=0,ntau,ntauwr)
c
c  Step forward in time and compute the solution.
c
      ii = 1
      jj = 0
 1    continue
         if(impulsive .and. i .gt. 0)goto 4
         call inject(ii,jj,t)
 4       i = i+1
         t = i*dt
c
c  Check to see if we need a report.  If time step
c  steps over the report time readjust the time step.
c
         if(t .ge. tr(itr))then
            dtnew = tr(itr) - (i-1)*dt
            dyold = dy
            dy = ttoy*dtnew
c
c  update with smaller step to the reporting time
c
            call tauup(jj,ii,symloop)
            call eup(ii,jj)
            call bcup(jj,ii)
            call muup(ii,jj)
c
c  output the distribution
c
            write(10,*) tr(itr), i
            do 2 m=0,ntau,ntauwr
               do 3 l=-nmu,nmu
                  write(10,202) (phi(k,l,m,jj), k=0,newr)
 3             continue
 2          continue
c
c  get back on step
c
            dtnew = i*dt - tr(itr)
            dy = ttoy*dtnew
c
c  update
c
            call tauup(jj,ii,symloop)
            call eup(ii,jj)
            call bcup(jj,ii)
            call muup(ii,jj)
c
c  advance to next reporting time
c  reset dy
c
            itr = itr + 1
            dy = dyold
         else
c
c  update
c
            call tauup(jj,ii,symloop)
            call eup(ii,jj)
            call bcup(jj,ii)
            call muup(ii,jj)
         endif
      if(t .lt. tmax)goto 1
      close(unit=10)
      end