      subroutine bcup(in,out)
c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the converging
c  magnetic field term by Barton's monotonic transport method.  
c
c  The points abs(mu) = 1 can be solved for exactly and are in 
c  this subroutine.
c
      include 'fkrplk.h'


      common /ce/ de(nemax), beta2(0:nemax), e(-1:nemax+1)
      common /cmu/ dmu(-nmumax+1:nmumax), onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax),rb3g2(0:nemax), bconv(0:ntaumax)
      integer in,out
      real tphi(-nmumax:nmumax),tphi1(-nmumax:nmumax)
      do 1 k=0,ne
         a = dy*beta(k)
         do 2 m=0,ntau
            do 6 l=-nmu,nmu
               tphi(l)=phi(k,l,m,in)
 6          continue
            aa = a*bconv(m)
            if(aa .eq. 0.)then
               do 5 l=-nmu,nmu
                  tphi1(l) = tphi(l)
 5             continue
            else
               psf = exp(-aa)
               tphi1(nmu) = psf*tphi(nmu)
               tphi1(-nmu) = tphi(-nmu)/psf
               if(aa .gt. 0.)then
                  fl = tphi(nmu)*onmu22(nmu)
                  do 3 l=nmu-1,-nmu+1,-1
                     lm1 = l-1
                     lp1 = l+1
                     flp = fl
                     d3 = tphi(l)
                     d2 = 0.5*(d3+tphi(lm1))
                     delt = 0.5*dmu(l)/dmu(lp1)
                     d1 = d3*(1.+delt)-tphi(lp1)*delt
                     if(d3 .le. tphi(lm1)) then
                        dk = amax1(d3,amin1(d1,d2))
                     else
                        dk = amin1(d3,amax1(d1,d2))
                     endif
                     fl = dk*onmu22(l)
                     tphi1(l)=
     +                    d3+aa*(flp-fl)/(amu(lp1)-amu(lm1))
 3                continue
               else
                  fl = tphi(-nmu)*onmu22(-nmu+1)
                  do 4 l=-nmu+1,nmu-1
                     lm1 = l-1
                     lp1 = l+1
                     flp = fl
                     d3 = tphi(l)
                     d2 = 0.5*(tphi(lp1)+d3)
                     delt = 0.5*dmu(lp1)/dmu(l)
                     d1 = d3*(1.+delt)-tphi(lm1)*delt
                     if(d3 .le. tphi(lm1)) then
                        dk = amin1(d3,amax1(d1,d2))
                     else
                        dk = amax1(d3,amin1(d1,d2))
                     endif
                     fl = dk*onmu22(lp1)
                     tphi1(l)=
     +                    d3+aa*(fl-flp)/(amu(lp1)-amu(lm1))
 4                continue
               endif
            endif
            do 7 l=-nmu,nmu
               phi(k,l,m,out)=tphi1(l)
 7          continue
 2       continue
 1    continue


      return
      end
