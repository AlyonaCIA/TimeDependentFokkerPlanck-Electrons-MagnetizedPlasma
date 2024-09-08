      subroutine init(rhomax,soft,tau)
c
c  This is the initialization subroutine.  The values for emin, emax,
c  ne, amu, nmu, taumin, taumax, and ntau are either read in or 
c  assume the default values set here in the data statement. 
c  This subroutine also sets dy the time step variable.  The arrays 
c  de, dtau, beta, beta2, rb3g2, rlmbda, rlmbda2, dmu, onmu2, and 
c  bconv are filled.  The following external subroutines  are needed:
c  dands(tau,rho,s)-- returns the electron number density rho and 
c  spatial postion s at a given column depth tau.
c  bfield(s,b,dbds)-- returns the value of dln(B)/ds == dbds and b 
c  at postion s.
c
c  Energy steps increase linearly from e(0)/4 up to emax so that 
c  ne steps are taken.
C  Added calculation tau as an output parameter
c Changes by Leon Ofman:
c 5-2-2001 check that dimensions do not exceed maximal dimensions
c
      parameter(a = 1.e24, aloglmda = 20.)
      include 'fkrplk.h'
      common /ctau/ dtau(ntaumax),rlmbda2(ntaumax)
      common /bmag/ bm(0:ntaumax)
      common /ce/ de(nemax),beta2(0:nemax),e(-1:nemax+1)
      common /cmu/ dmu(-nmumax+1:nmumax),onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax), rb3g2(0:nemax), bconv(0:ntaumax)


      real soft(0:ntaumax), tau(0:ntaumax)
c
c  Set values for grid parameters.
c
      data emin,emax,ned,nmud,taumin,taumax,ntaud  
     +     / 10.,10000.,120,60,1.e-5,2.e-3,120 /
      nmu = nmud
      if(nmu.gt.nmumax) then
        print*,'ERROR: nmu>nmumax! increase nmumax in fkrplk.h'
        stop
      endif
c
      ntau = ntaud
      if(ntau.gt.ntaumax) then
        print*,'ERROR: ntau>ntaumax! increase ntaumax in fkrplk.h'
        stop
      endif
c
      ne = ned
      if(ne.gt.nemax) then
        print*,'ERROR: ne>nemax! increase nemax in fkrplk.h'
        stop
      endif
c
c  Set up "energy" dependent arrays.
c
      e(0)=emin/511.
      deo=e(0)/4.0
      g=e(0)+1.
      beta(0)=sqrt(1.0 - 1.0/(g**2))
      rb3g2(0)=1./((beta(0)**3)*g*g)
      e(-1) = e(0)-deo
      g = 0.5*(e(0)+e(-1)) + 1.
      beta2(0) = sqrt(1.0-1.0/(g*g))
      des = (emax-emin)/(511.*ne)
      if(des .le. deo)then
         do 1 k=1,ne
            e(k) = e(0) + k*des
            de(k)=e(k)-e(k-1)
            g = e(k) + 1.0 
            beta(k)=sqrt(1.0 - 1.0/(g*g))
            rb3g2(k)=1.0/( (beta(k)**3)*g*g )
            g= 0.5*( e(k)+e(k-1) ) + 1.0
            beta2(k)=sqrt(1.0-1.0/(g*g))
 1       continue
         e(ne+1) = e(ne) + des
      else
         slope = 0.000001
         r = 10.
         delt = emax/511. - e(0) - deo*ne
 11      if(511.*delt/emax .gt. 1.e-4)then
            slope = slope*r
            sum = 1.
            do 12 k=1,ne-1
               sum = sum+(slope+1.)**k
 12         continue
            delt = emax/511. - e(0) - deo*sum
            goto 11
         elseif(delt .lt. 0.)then
            slope = slope/r
            r = sqrt(r)
            delt = e(0)
            goto 11
         endif
         do 13 k=1,ne
            de(k) = deo + (e(k-1)-e(0))*slope
            e(k) = e(k-1) + de(k)
            g = e(k) + 1.0 
            beta(k)=sqrt(1.0 - 1.0/(g*g))
            rb3g2(k)=1.0/( (beta(k)**3)*g*g )
            g= 0.5*( e(k)+e(k-1) ) + 1.0
            beta2(k)=sqrt(1.0-1.0/(g*g))
 13      continue
         e(ne+1) = e(ne) + de(ne)
      endif
c 
c  Set up "mu" dependent arrays.
c
      dr=1./nmu
      amu(-nmu)=-1.0
      amu(nmu) = 1.0
      onmu2(-nmu)=0.0
      onmu2(nmu) =0.0
      do 2 l=-nmu+1,nmu-1
c         amu(l)=cos(1.570796327*(1.-dr*l))
         amu(l) = real(l)/real(nmu)
         onmu2(l)=1.-amu(l)*amu(l)
         onmu22(l) = 1.-0.25*(amu(l)+amu(l-1))**2
         dmu(l)=amu(l)-amu(l-1)
 2    continue
      dmu(nmu)=1.-amu(nmu-1)
      onmu22(nmu) = 1.-0.25*(amu(nmu)+amu(nmu-1))**2
      amu(0) = 0.
      onmu2(0) = 1.
c
c  Set up "tau" dependent arrays.
c
c    For logarithmic grid spacing: taus=taumax/taumin
      taus=taumax/ntau
      to = 0.0
      soft(0) = to
      tau(0) = to
      call dands(taumax,rhomax,s)
      call dands(to,rho,s)
      bfvar = a/(aloglmda*rhomax)
      rlmbda(0)=rho/rhomax
      call bfield(to,b,dbds)
      bconv(0) = bfvar*dbds
      bm(0) = b
      do 3 m=1,ntau
c    For logarithmic grid spacing:
c                  tn=taumin*taus**(real(m-1)/real(ntau-1))
         tn=taus*m
         dtau(m)=tn-to
         call dands(tn,rho,s)
         soft(m) = s
         tau(m) = tn
         rlmbda(m)=rho/rhomax
         call bfield(s,b,dbds)
         bconv(m) = bfvar*dbds
         bm(m) = b
         call dands( 0.5*(tn+to),rho,s)
         rlmbda2(m)=rho/rhomax
         to=tn
 3    continue


c
c  Set dy to the the value required for stability of 
c  all subroutines.
c
c  collision terms
c
      dycoul = 0.9*amin1( 0.5*dmu(nmu)/rb3g2(0), beta(0)*de(1) )
c 
c  advective term
c     
      dytau = dtau(1) / rlmbda2(1)
      do 4 m=1,ntau
         dytau = amin1( dytau, dtau(m)/rlmbda2(m) )
 4    continue
c
c  mirroring term
c
      dymir = 1.e20
      if(bm(ntau) .eq. bm(0))goto 7
      do 5 l=0,nmu-1
         dymir = amin1(dymir,dmu(l+1)/onmu2(l))
 5    continue
      dymiro = dymir
      if(bconv(ntau) .ne. 0.)dymir = dymiro/bconv(ntau)
      do 6 m=0,ntau-1
         if(bconv(m) .ne. 0.)dymir = amin1(dymir,dymiro/bconv(m))
 6    continue
      dymir = 0.5*dymir/beta(ne)
 7    dy = amin1(dymir,dytau)
      dy = amin1(dy,dycoul)
c
      return
      end