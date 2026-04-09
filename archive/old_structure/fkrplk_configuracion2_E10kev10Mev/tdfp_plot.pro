pro tdfp_plot


;+
;Reads and plots output from the time-dependent Fokker-Planck code
;
;Electron distribution function
; phi(Energy, pitch angle, position, time)
;
;Revised 10 September 2001 by Gordon Holman
;-


;Read and store data from output file

openr, ten, 'fkrplk.test', /get_lun
readf, ten, dt, tmax, dy, rhomax
readf, ten, newr, nmu, ntau, ntr
newr = fix(newr) & nmu = fix(nmu) & ntau = fix(ntau) & ntr = fix(ntr)
ekev = fltarr(newr+1)
amu = fltarr(2*nmu+1)
dtau = fltarr(ntau)
tau = dtau
soft = dtau & bm = fltarr(ntau+1)
tr = fltarr(ntr+1)


readf, ten, ekev
eta = double(ekev/510.99)
beta = sqrt(1.0-1.0/(eta+1.0)^2)
eta = eta^2/(eta+1.0)


readf, ten, amu
pad = 180.0*acos(amu)/!Pi
readf, ten, dtau
readf, ten, tau
tau = [0.0, tau]
readf, ten, soft
soft = [0.0,soft]
readf, ten, bm


phi = fltarr(newr+1, 2*nmu+1, ntau+1, ntr+1)
tphi = phi(*, *, *, 0)
i = fltarr(ntr+1)
temp = fltarr(2)
FOR itr = 0, ntr DO BEGIN
   readf, ten, temp
   tr(itr) = temp(0)
   i(itr) = temp(1)
   readf, ten, tphi
   phi(*, *, *, itr) = tphi
ENDFOR
free_lun, ten


;Specify the values of the indices that remain fixed

paindex=0
sindex=10
eindex=10

;Plot data

phiE0 = Reform(phi(*,paindex,sindex,0))
phiE1 = Reform(phi(*,paindex,sindex,1))
phiE2 = Reform(phi(*,paindex,sindex,2))
phiE3 = Reform(phi(*,paindex,sindex,3))
phiE4 = Reform(phi(*,paindex,sindex,4))


Window,0,Title = 'Energy Dependence'
plot, ekev, phiE0,/XLog,/YLog,XTitle = 'Energy (keV)', $
  YTitle = 'Normalized Density', $
  Title = 'Electron Distribution Function Energy Dependence',$
  XRange=[10.,1000.],YRange=[1.E-30,1.E00],XStyle=1,YStyle=1
oplot, ekev, phiE1, LineStyle = 1, color=112
oplot, ekev, phiE2, LineStyle = 2, color=208
oplot, ekev, phiE3, LineStyle = 3, color=165
oplot, ekev, phiE4, LineStyle = 4, color=254

;Create and position Legend for electron distribution function vs. energy graph
XYOutS, 1.0E+02, 1.0E-02, 'Pitch Angle = '+ strtrim((pad(paindex)),2)+' degrees', Size = 1.3
XYOutS, 1.0E+02, 1.0E-04, 'Position = '+ strtrim((soft(sindex)),2)+' km', Size = 1.3


phiPA0 = Reform(phi(eindex,*,sindex,0))
phiPA1 = Reform(phi(eindex,*,sindex,1))
phiPA2 = Reform(phi(eindex,*,sindex,2))
phiPA3 = Reform(phi(eindex,*,sindex,3))
phiPA4 = Reform(phi(eindex,*,sindex,4))


Window,1,Title = 'Pitch Angle Dependence'
plot, pad, phiPA0,/YLog,XTitle = 'Pitch Angle (degrees)', $
  YTitle = 'Normalized Density', $
  Title = 'Electron Distribution Function Pitch Angle Dependence',$
  XRange=[0.,180.],YRange=[1.E-30,1.E00],XStyle=1,YStyle=1
oplot, pad, phiPA1, LineStyle = 1, color=112
oplot, pad, phiPA2, LineStyle = 2, color=208
oplot, pad, phiPA3, LineStyle = 3, color=165
oplot, pad, phiPA4, LineStyle = 4, color=254

;Create and position Legend for electron distribution function vs. pitch angle graph
XYOutS, 90., 1.0E-02, 'Energy = '+ strtrim((ekev(eindex)),2)+' keV', Size = 1.3
XYOutS, 90., 1.0E-04, 'Position = '+ strtrim((soft(sindex)),2)+' km', Size = 1.3


phiS0 = Reform(phi(eindex,paindex,*,0))
phiS1 = Reform(phi(eindex,paindex,*,1))
phiS2 = Reform(phi(eindex,paindex,*,2))
phiS3 = Reform(phi(eindex,paindex,*,3))
phiS4 = Reform(phi(eindex,paindex,*,4))


Window,2,Title = 'Spatial Dependence'
plot, soft, phiS0,/YLog,XTitle = 'Distance (cm)', $
  YTitle = 'Normalized Density', $
  Title = 'Electron Distribution Function Spatial Dependence',$
  XRange=[0.,1.E+09],YRange=[1.E-30,1.E00],XStyle=1,YStyle=1
oplot, soft, phiS1, LineStyle = 1, color=112
oplot, soft, phiS2, LineStyle = 2, color=208
oplot, soft, phiS3, LineStyle = 3, color=165
oplot, soft, phiS4, LineStyle = 4, color=254

;Create and position Legend for electron distribution function vs. distance graph
XYOutS, 5.e8, 1.0E-2, 'Energy = '+ strtrim((ekev(eindex)),2)+' keV', Size = 1.3
XYOutS, 5.e8, 1.0E-4, 'Pitch Angle = '+ strtrim((pad(paindex)),2)+' deg', Size = 1.3

;Specify a vector of position indices for the Time Dependence plot
splt = [0,5,10,20,30]

phiT0 = Reform(phi(eindex,paindex,splt(0),*))
phiT1 = Reform(phi(eindex,paindex,splt(1),*))
phiT2 = Reform(phi(eindex,paindex,splt(2),*))
phiT3 = Reform(phi(eindex,paindex,splt(3),*))
phiT4 = Reform(phi(eindex,paindex,splt(4),*))


Window,3,Title = 'Time Dependence'
plot, tr, phiT0,/YLog,XTitle = 'Time (s)', $
  YTitle = 'Normalized Density', $
  Title = 'Electron Distribution Function Time Dependence',$
  XRange=[0.,11.],YRange=[1.E-30,1.E00],XStyle=1,YStyle=1
oplot, tr, phiT1, LineStyle = 1, color=112
oplot, tr, phiT2, LineStyle = 2, color=208
oplot, tr, phiT3, LineStyle = 3, color=165
oplot, tr, phiT4, LineStyle = 4, color=254

;Create and position Legend for electron distribution function vs. time graph
XYOutS, 5.5, 1.0E-02, 'Energy = '+ strtrim((ekev(eindex)),2)+' keV', Size = 1.3
XYOutS, 5.5, 1.0E-04, 'Pitch Angle = '+ strtrim((pad(paindex)),2)+' deg', Size = 1.3

tlx=0.50
tly=0.50
PlotS, [tlx, tlx+0.05], [tly,tly], LineStyle=0, color=255, /Normal
XYOutS, tlx+0.06, tly, 'Position= '+strtrim((soft(splt(0))),2)+' km', Size = 1.3, /Normal
PlotS, [tlx, tlx+.05], [tly-.05,tly-.05], LineStyle=1, color=112, /Normal
XYOutS, tlx+0.06, tly-.05, 'Position= '+strtrim((soft(splt(1))),2)+' km', Size = 1.3, /Normal
PlotS, [tlx, tlx+0.05], [tly-0.1,tly-0.1], LineStyle=2, color=208, /Normal
XYOutS, tlx+0.06, tly-0.1, 'Position= '+strtrim((soft(splt(2))),2)+' km', Size = 1.3, /Normal
PlotS, [tlx, tlx+0.05], [tly-0.15,tly-0.15], LineStyle=3, color=165, /Normal
XYOutS, tlx+0.06, tly-0.15, 'Position= '+strtrim((soft(splt(3))),2)+' km', Size = 1.3, /Normal
PlotS, [tlx, tlx+0.05], [tly-0.2,tly-0.2], LineStyle=4, color=254, /Normal
XYOutS, tlx+0.06, tly-0.2, 'Position= '+strtrim((soft(splt(4))),2)+' km', Size = 1.3, /Normal


;Optional plot of optical depth vs. distance
; along magnetic field line


;Window,4,Title = 'Tau vs. Distance',retain=2
;plot,soft,tau,Xtitle='Distance (cm)',$
;  YTitle='Optical Depth Tau', $
;  Title = 'Optical Depth vs. Distance',$
;  XRange=[0.,1.E+09], YRange=[0.,0.002],XStyle=1,YStyle=1


END