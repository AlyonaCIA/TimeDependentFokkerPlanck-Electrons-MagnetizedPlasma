      subroutine tauup(in,out,symloop)


c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the tau term using
c  monotonic transport.  Boundary conditions for positive
c  velocity at m=0 are either symetric or no flux comes in from m < 0;
c  for negative velocity the boundary condition at m=ntau is that
c  no flux comes in.
c  Note:  Nonuniform step sizes in tau are not allowed. 
c
      include 'fkrplk.h'
      common /ctau/ dtau(ntaumax), rlmbda2(ntaumax)
      real tphi(-nmumax:nmumax,0:ntaumax),
     +     tphi1(-nmumax:nmumax,0:ntaumax)
      logical symloop
      integer in,out
 
      a0 = dy/dtau(1)
      do 1 k=0,ne
         a1 = a0*beta(k)
         do 7 m=0,ntau
            do 6 l=-nmu,nmu
               tphi(l,m) = phi(k,l,m,in)
               tphi1(l,m) = tphi(l,m)
 6          continue
 7       continue
c
c  NEGATIVE VELOCITIES
c
         do 2 l=-nmu,-1
            a = a1*amu(l)
            a1pa = a*(1.+a)
c
c  m=ntau:  For this point we use the boundary condition that 
c           phi = 0. for m > ntau
c
            phic = tphi(l,ntau)
            phin = tphi(l,ntau-1)
            dphi = phic - phin
            delt = -dphi*phic
            if(delt .gt. 0.)then
               delt = -delt/phin
               tphi1(l,ntau) = phic*(1.+a) - a1pa*delt
            else
               delt = 0.
               tphi1(l,ntau) = phic*(1.+a) 
            endif


            do 3 m=ntau-1,1,-1
               phip = phic
               phic = phin
               phin = tphi(l,m-1)
               dphip = dphi
               dphi = phic - phin
               deltp = delt
               delt = dphi*dphip
               if(delt .gt. 0.)then
                  delt=delt/(phip - phin)
               else
                  delt = 0.
               endif
               tphi1(l,m) = phic - a*dphip + a1pa*(deltp - delt)
 3          continue
c
c  m=0:  no problem for negative velocities, but the correction
c        to the upwind method depends on the symmetry of the loop.
c
            phip = phic
            phic = phin
            dphip = dphi
            deltp = delt
            if(symloop)then
               dphi = phic - tphi(-l,1)
               delt = dphi*dphip
               if(delt .gt. 0.)then
                  delt = delt/(phip - tphi(-l,1))
               else
                  delt = 0.
               endif
            else
               delt = tphi(l,0)*dphip
               if(delt .gt. 0.)then
                  delt = delt/tphi(l,1)
               else
                  delt = 0.
               endif
            endif
            tphi1(l,0) = phic - a*dphip + a1pa*(deltp - delt)
 2       continue


c
c  POSITIVE VELOCITIES
c
         do 8 l=1,nmu
            a = a1*amu(l)
            a1ma = a*(1.-a)
c
c  m=0:  The treatment of this point depends on the loop symmety.
c        If symmetric, we use phi(k,l,-m)=phi(k,-l,m) for m > 0.
c        If not symmetric, then the we use phi = 0. for m < 0.
c
            phic = tphi(l,0)
            phin = tphi(l,1)
            if(symloop)then
               dphi = phic - tphi(-l,1)
               dphip = phin - phic
               delt = dphi*dphip
               if(delt .gt. 0.)then
                  delt = delt/(phin - tphi(-l,1))
               else
                  delt = 0.
               endif
               delto = dphi*(tphi(-l,1)-tphi(-l,2))
               if(delto .gt. 0.)then
                  delto = delto/(phic - tphi(-l,2))
               else
                  delto = 0.
               endif
               tphi1(l,0) = phic - a*dphi - a1ma*(delt - delto)
            else
               dphip = phin - phic
               delt = phic*dphip
               if(delt .gt. 0.)then
                  delt = delt/phin
               else
                  delt = 0.
               endif
               tphi1(l,0) = phic*(1.-a) - a1ma*delt
            endif


            do 9 m=1,ntau-1
               phip = phic
               phic = phin
               phin = tphi(l,m+1)
               delto = delt
               dphi = dphip
               dphip = phin - phic
               delt = dphi*dphip
               if(delt .gt. 0.)then
                  delt=delt/(phin - phip)
               else
                  delt = 0.
               endif
               tphi1(l,m) = phic - a*dphi - a1ma*(delt - delto)
 9          continue
            phip = phic
            phic = phin
            delto = delt
            dphi = dphip
            delt = -phic*dphi
            if(delt .gt. 0.)then
               delt = -delt/phip
            else
               delt = 0.
            endif
            tphi1(l,ntau)=phic- a*dphi - a1ma*(delt-delto)
 8       continue
         do 17 m=0,ntau
            do 16 l=-nmu,nmu
               phi(k,l,m,out) = tphi1(l,m)
 16         continue
 17      continue
 1    continue
      return
      end 