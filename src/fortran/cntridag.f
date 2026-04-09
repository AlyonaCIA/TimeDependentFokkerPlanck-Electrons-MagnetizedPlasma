      subroutine cntridag(a,b,c,r,n,tphi1)
c
c  Solves the tridiagonal matrix equation which results from the 
c  Crank-Nicholson method used in subroutine muup for phi(k,l,m,out) 
c  for given k, m, and out at all l.  
c  Adapted from Press, W.H., Flannery, B.P., Teukolsky, S.A., and 
c  Vetterling, W.T. 1986, "Numerical Recipes", (New York: Cambridge
c  University Press), p. 40.
c
      include 'fkrplk.h'
      real gam(-nmumax:nmumax),a(-nmumax:nmumax),b(-nmumax:nmumax),
     +     c(-nmumax:nmumax),r(-nmumax:nmumax)
      real tphi1(-nmumax:nmumax)


      if (b(-n).eq.0) stop
      bet=b(-n)
      tphi1(-n)=r(-n)/bet
      do 11 j=-n+1,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
c--      if (bet.eq.0) pause
         tphi1(j)=(r(j)-a(j)*tphi1(j-1))/bet
 11   continue
      do 12 j=n-1, -n,-1
         tphi1(j)=tphi1(j)-gam(j+1)*tphi1(j+1)
 12   continue
      return
      end