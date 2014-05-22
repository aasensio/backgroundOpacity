program test
use chemical_eq, only : calculate_abundance_Pe_from_T_Pg, calculate_abundance_Pg_from_T_Pe, calculate_abundance_from_T_Pe_Pg
use background_opacity_module, only : background_opacity

	integer :: ngrid, i, n
	real(kind=8), dimension(1) :: h, T, PT, ne, PH_out, PHminus_out, PHplus_out, &
		PH2_out, PH2plus_out, res, opa1, opa2, Pe, n_htot
	real(kind=8) :: lambda
	character(len=40) :: file

	open(unit=12,file='physical.out',action='write',status='replace')
	write(12,*) '   T         Pe        P_T'
	open(unit=13,file='physical.dat',action='read',status='old')
	read(13,*)

	read(13,*) n

	do i = 1, n
		read(13,*) T, Pe, PT

		print *
		print *, 'INPUT'
		print *, 'T = ', T
		print *, 'Pe = ', Pe
		print *, 'P_T = ', PT
		
		ngrid = 1
		h = 1.d0
		file = 'test'

! n_htot is the total number of hydrogen nuclei
! 		res = calculate_abundance_Pe_from_T_Pg(15, ngrid, h, T, PT, file, PH_out, PHminus_out, &
! 			PHplus_out, PH2_out, PH2plus_out, Pe)
! 
! 		n_htot = (PH_out(1)+PHminus_out(1)+PHplus_out(1)+PH2_out(1)+PH2plus_out(1)) / (1.3806503d-16 * T(1))
! 		
! 		print *, 'OUTPUT'
! 		print *, 'n(CO)/nHtot = ', res(1) / n_htot(1)
! 		print *, 'nH/nHtot   = ', PH_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
! 		print *, 'nH-/nHtot  = ', PHminus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
! 		print *, 'nH+/nHtot  = ', PHplus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
! 		print *, 'nH2/nHtot  = ', PH2_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
! 		print *, 'nH2+/nHtot = ', PH2plus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
!  		print *, 'Pe = ', Pe(1)
!  		print *, 'ne/nHtot = ', Pe(1) / (n_htot(1)*1.3806503d-16 * T(1))
!  		pause


! First compute the partial pressures of hydrogen
 		res = calculate_abundance_from_T_Pe_Pg(15, ngrid, h, T, Pe, PT, file, PH_out, PHminus_out, &
			PHplus_out, PH2_out, PH2plus_out)

		n_htot = (PH_out(1)+PHminus_out(1)+PHplus_out(1)+PH2_out(1)+PH2plus_out(1)) / (1.3806503d-16 * T(1))

		print *, 'OUTPUT'		
		print *, 'nH/nHtot   = ', PH_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
		print *, 'nH-/nHtot  = ', PHminus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
		print *, 'nH+/nHtot  = ', PHplus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
		print *, 'nH2/nHtot  = ', PH2_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
		print *, 'nH2+/nHtot = ', PH2plus_out(1) / (n_htot(1)*1.3806503d-16 * T(1))
 		print *, 'Pe = ', Pe(1)
 		print *, 'ne/nHtot = ', Pe(1) / (n_htot(1)*1.3806503d-16 * T(1))

 		
! Compute the opacities at 5000 and 15000 A
		opa1 = background_opacity(T(1), Pe(1), PH_out(1), PHminus_out(1), PHplus_out(1), PH2_out(1), PH2plus_out(1), 5000.d0)
		opa2 = background_opacity(T(1), Pe(1), PH_out(1), PHminus_out(1), PHplus_out(1), PH2_out(1), PH2plus_out(1), 15000.d0)

		write(12,FMT='(F8.2,2X,E12.5,2X,E15.9,2X,E15.9)') T, Pe, opa1 / n_htot, opa2 / n_htot		

	enddo
	close(12)
	close(13)

end program test
