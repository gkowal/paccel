part = 'p'

tm = readfits(part + '_tim.fits.gz')
en = readfits(part + '_ene.fits.gz')
xp = readfits(part + '_pos.fits.gz')
vp = readfits(part + '_vel.fits.gz')

vv = sqrt(total(vp[*,*]^2,1))

;spawn, "grep 'c  ' ./params.in", res
;c = (float(strmid(res, strpos(res, "=")+1)))[0]
c = 10.0

spawn, "grep 'tunit ' ./params.in", res
tunit = (strmid(res, strpos(res, "=")+1))[0]
case strtrim(tunit) of
" 's'": tunit = 'secs'
" 'm'": tunit = 'mins'
" 'h'": tunit = 'hours'
" 'd'": tunit = 'days'
" 'y'": tunit = 'yrs'
else:
endcase

cs=2.5

;goto, anim

xr=[0,1.0]*1.e1+0.0
xr=[69.0,69.5]
xr=[30,40]
xr=[0.0,100.0]

!p.multi = [4,1,4]
erase
plot, tm, vv/c, xr=xr, charsize=cs, ytit='Particle speed [c]', xs=1, yr=[0.0,1], ys=1
plot, tm, deriv(tm,vv/c)/sqrt(1.0-(vv/c)^2), xr=xr, charsize=cs, ytit='Particle acceleration', xs=1, yr=[-10,10], ys=1
oplot, xr, [0.0,0.0], color=150
;plot, tm, 2.*((xp[1,*]+1.)/2 - floor((xp[1,*]+1.)/2)) - 1.0, xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3
;plot, tm, ((xp[2,*]+0.5) - floor((xp[2,*]+0.5))) - 0.5, xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3
plot, tm, xp[1,*], xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3, ys=1, yr=[-1,1]*5
oplot, tm, xp[2,*]
oplot, tm, xp[0,*]
plot, tm, en, xtit='Time [' + tunit + ']', ytit='E [MeV]', charsize=cs, /yl, xr=xr, xs=1, yr=[1,1e14], ys=1, psym=1, symsize=0.5

goto, fine
anim:
for i = 0, 995 do begin

xr=[0,5]+0.1*i

!p.multi = [4,1,4]
erase
plot, tm, vv/c, xr=xr, charsize=cs, ytit='Particle speed [c]', xs=1, yr=[0,1], ys=1, /yl
plot, tm, deriv(tm,vv/c)/sqrt(1.0-(vv/c)^2), xr=xr, charsize=cs, ytit='Particle acceleration', xs=1, yr=[-10,10], ys=1
oplot, xr, [0.0,0.0], color=150
;plot, tm, 2.*((xp[1,*]+1.)/2 - floor((xp[1,*]+1.)/2)) - 1.0, xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3
;plot, tm, ((xp[2,*]+0.5) - floor((xp[2,*]+0.5))) - 0.5, xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3
plot, tm, xp[1,*], xr=xr, charsize=cs, ytit='Y position', xs=1, psym=3, ys=1, yr=[-1,1]
plot, tm, en, xtit='Time [' + tunit + ']', ytit='E [MeV]', charsize=cs, /yl, xr=xr, xs=1, yr=[1,1e14], ys=1

wait, 0.1
endfor

fine:

end
