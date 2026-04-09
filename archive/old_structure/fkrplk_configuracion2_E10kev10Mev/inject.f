      subroutine inject(in,out,t)
c
c  This subroutine injects an electron distribution which is:
c     initially a gaussian in space, but becomes a delta function in space 
c     isotropic in pitch-angle
c     power-law in energy
c     a triangular-shaped pulse, 2 seconds long
c
      include 'fkrplk.h'
      common /ctau/ dtau(ntaumax), rlmbda2(ntaumax)
      integer out
c
      delta = 5.
      e0 = 10./511.
      sa = 5.0e-2
c
      if(t .le. 1.0) then
         ft0 = 0.1+0.9*t
      else
         ft0 = 1.9-0.9*t
      endif
      if(ft0 .le. 0) ft0 = 0.0
c
      if(t .le. 0.)then
         do k=0,ne
            do m=0,ntau
               do l=-nmu,nmu
                 e = 1./sqrt(1.-beta(k)**2)-1.
                 g = (e0/e)**delta
                 s = real(m)/real(ntau)
                 f = exp(-(s/sa)**2)
                 phi(k,l,m,out) = g*f*ft0
               enddo
            enddo
         enddo
      else
         do k=0,ne
            do l=-nmu,nmu
               e = 1./sqrt(1.-beta(k)**2)-1.
               g = (e0/e)**delta
               phi(k,l,0,out) = g*ft0
            enddo
         enddo
      endif
      return 
      end