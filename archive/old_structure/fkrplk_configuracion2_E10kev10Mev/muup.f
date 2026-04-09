      subroutine muup(in,out)


c     phi(k,l,m,in) is updated to phi(k,l,m,out) for the pitch angle 
c     diffusion term.  The Crank-Nicholson method is employed and 
c     subroutine cntridag is called to solve the tridiagonal matrix.  


      include 'fkrplk.h'
      common /cmu/ dmu(-nmumax+1:nmumax), onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax), rb3g2(0:nemax), bconv(0:ntaumax)


      real a(-nmumax:nmumax),b(-nmumax:nmumax),c(-nmumax:nmumax),
     +     r(-nmumax:nmumax)
      real tphi(-nmumax:nmumax),tphi1(-nmumax:nmumax)
      integer in,out
      do 10 k=0,ne
         b1=beta(k)*.25
         do 20 m=0,ntau
            b4=rlmbda(m)*rb3g2(k)
            do 25 l=-nmu,nmu
               tphi(l)=phi(k,l,m,in)
 25         continue
            do 30 l=-nmu+1,nmu-1
               lp1=l+1
               lm1=l-1
               dydsum=dy/(dmu(lp1)+dmu(l))
               b5=b4*onmu2(l)
               aa=b5*dydsum
               bb1=(-b4*amu(l))*dydsum
c     a(l), b(l), c(l) are the elements of the tridiagonal matrix. 
c     r(l) is the vector product of the matrix and the vector 
c     of updated phi(k,l,m,out).  
               a(l)=bb1-aa/dmu(l)
               b(l)=1.+b5*dy/(dmu(l)*dmu(lp1))
               c(l)=-bb1-aa/dmu(lp1)
               r(l)=aa*((tphi(lp1)-tphi(l))/dmu(lp1)-
     +              (tphi(l)-tphi(lm1))/dmu(l))+bb1*
     +              (tphi(lp1)-tphi(lm1))+tphi(l)
 30            continue
c     at amu=-1, the only term remaining is the first derivative 
c     from the pitch angle diffusion term which we approximate by 
c     a forward difference
            c(-nmu)=-b4*dy/dmu(-nmu+1)
            b(-nmu)=1.-c(-nmu)
            r(-nmu)=-c(-nmu)*(tphi(-nmu+1)-tphi(-nmu))
     +           +tphi(-nmu)
c     at amu=1, use the backward difference
            a(nmu)=-b4*dy/dmu(nmu)
            b(nmu)=1.-a(nmu)
            r(nmu)=a(nmu)*(tphi(nmu)-tphi(nmu-1))
     +           +tphi(nmu)
            call cntridag(a,b,c,r,nmu,tphi1)
            do 35 l=-nmu,nmu
               phi(k,l,m,out)=tphi1(l)
 35         continue
 20         continue
 10      continue
      return
      end