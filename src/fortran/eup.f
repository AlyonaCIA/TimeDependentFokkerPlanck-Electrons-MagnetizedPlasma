      subroutine eup(in,out)
c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the energy term by
c  Barton's method of monotonic transport.
c
      include 'fkrplk.h'
      common /ce/ de(nemax), beta2(0:nemax), e(-1:nemax+1)
      integer in,out
      real tphi(0:nemax),tphi1(0:nemax)
      do 1 m=0,ntau
         a = 2.*dy*rlmbda(m)
         do 2 l=-nmu,nmu
            do 6 k=0,ne
               tphi(k)=phi(k,l,m,in)
 6          continue
c     
c     k=ne:  For this point the best value of d2 or d3
c            is used as the transported value; d1 is not 
c            computable.
c
            k = ne
            km1 = k-1
            d3 = tphi(k)
            d2 = 0.5*(tphi(k)+tphi(km1))
            if(tphi(k) .le. tphi(km1)) then
               dk = amax1(d3,d2)
            else
               dk = amin1(d3,d2)
            endif
            f1 = dk/beta2(k)
            tphi1(k) = d3 - a*f1/(e(ne+1)-e(km1))
            do 3 k=ne-1,1,-1
               kp1 = k+1
               km1 = k-1
               f1p = f1
               d3 = tphi(k)
               d2 = 0.5*(d3 + tphi(km1))
               delt = 0.5*(e(k)-e(km1))/de(kp1)
               d1 = d3*(1.+delt)-tphi(kp1)*delt
               if(d3 .le. tphi(km1)) then
                  dk=amax1(d3,amin1(d1,d2))
               else
                  dk=amin1(d3,amax1(d1,d2))
               endif
               f1 = dk/beta2(k)
               tphi1(k) = d3 + a*(f1p-f1)/(e(kp1)-e(km1))
 3          continue
c
c     k=0:  For this point the best value between the upwind 
c           and d1 value is used.  Note that phi(0) is always
c           greater than phi(-1).
c
            k = 0
            kp1 = k+1
            km1 = k-1
            f1p = f1
            d3 = tphi(k)
            delt = 0.5*(e(k)-e(km1))/de(kp1)
            d1 = d3*(1.+delt)-tphi(kp1)*delt
            dk = amin1(d3,d1)
            f1 = dk/beta2(0)
            tphi1(k) = d3 + a*(f1p-f1)/(e(kp1)-e(km1))
            do 7 k=0,ne
               phi(k,l,m,out)=tphi1(k)
 7          continue
 2       continue
 1    continue
      return 
      end