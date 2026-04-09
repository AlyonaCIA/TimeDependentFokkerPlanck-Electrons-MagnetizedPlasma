      parameter(nemax=120,nmumax=60,ntaumax=120)
      common /univ/ phi(0:nemax,-nmumax:nmumax,0:ntaumax,0:1)
      common /univ1/ amu(-nmumax:nmumax), ne, nmu, ntau, 
     +       beta(0:nemax), rlmbda(0:ntaumax), dy