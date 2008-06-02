part = 'e'

tm = readfits(part + '_tim.fits.gz')
xp = readfits(part + '_pos.fits.gz')
vp = readfits(part + '_vel.fits.gz')

vv = sqrt(total(vp[*,*]^2,1))

spawn, "grep 'c  ' ./params.in", res
c = (float(strmid(res, strpos(res, "=")+1)))[0]

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

plot, tm, vv/c, xtit='Time [' + tunit + ']', ytit='v/c', charsize=1.4

end
