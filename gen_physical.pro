pro gen_physical
	T = findgen(16) * 500.d0 + 2500.d0
	Pe = 10.d0^(-3+findgen(7))

	openw,2,'physical.dat'
	for i = 0, 15 do begin
		for j = 0, 6 do begin
			res = kappa_c(Pe[j],T[i],5000.,htoverv)
			printf,2,T[i],Pe[j],htoverv
		endfor
	endfor
	close,2

	spawn, './backopa'

	dat = ddread('opacity.dat',offset=1,/count)
	T = reform(dat[0,*])
	Pe = reform(dat[1,*])
	opa = reform(dat[2,*])

	triangulate, T, Pe, triangles, boundary

 	grid = trigrid(T, Pe, opa, triangles)
	
	stop
end

pro gen_physical2
	Tmin = 5000.d0
	Tmax = 20000.d0
	logPgmin = -2.d0
	logPgmax = 5.d0

	n = 100
	T = findgen(n)/(n-1.d0) * (Tmax-Tmin) + Tmin
	Pg = 10.d0^(findgen(n)/(n-1.d0) * (logPgmax-logPgmin) + logPgmin)

	openw,2,'physical.dat'
	printf,2,'T      Pe        P_T'
	printf,2,n_elements(T)*n_elements(Pg)
	for i = 0, n-1 do begin
		for j = 0, n-1 do begin
			printf,2,T[i],1.d0,Pg[j]
		endfor
	endfor
	close,2

	stop
; 	spawn, './backopa'

	dat = ddread('opacity.dat',offset=1,/count)
	T = reform(dat[0,*])
	Pg = reform(dat[1,*])
	Pe = reform(dat[2,*])

	Pe = reform(Pe,n,n)

	stop
end