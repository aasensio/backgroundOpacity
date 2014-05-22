!*******************************************************************
! Global variables
!*******************************************************************
module variables
implicit none

! equil(9,273) : coefficients of the series expansion of the equilibrium constant
! molec(273) : name of the molecule as it appears in the file
! nombre(273) : name of the molecule
! elements(21) : name of the 21 atomic elements included
! abund(21) : abundance of the 21 atomic elements 
! pot_ion(21) : ionization potential of the 21 species included
! afinidad(21) : electronic affinity of the 21 species included
! estequio(i,j) : stequiometric coefficient of the species j into molecule i
! charge(273) : charge of the molecule
! composicion(i,j) : atoms which constitute molecule i (j=1..4)
! includ(273) : which molecules are included
! n_atoms_mol(273) : number of atoms in each molecule for its units transformation
! atomic_partition(i,j,k) : coefficient i of the atomic species k. The charge is indicated by j
!								(j = 0:neutral, 1:positive, 2:negative)
! equilibrium(i) : equilibrium constant at a given temperature for the whole set of molecules
! equilibrium_atomic(i,j) : equilibrium constant for the atomic species j for the positive (j=1)
!								ionization (A->A+ + e-) or the negative ionization (j=2) (A+e- ->A-)
! x_sol(21) : the abundance of every species during the iteration
! x0_sol(21) : the abundance of every species at the beggining of the iteration

	real(kind=8) :: equil(9,273)
	integer :: n_elements, print_abund, n_grid
	character(len=16) :: molec(273), nombre_mol(273)
	character(len=60) :: output_model
	character(len=2) :: elements(21)
	real(kind=8) :: abund_atom(21), pot_ion(21), afinidad(21)
	integer :: which_index_element(21) = (/1,2,6,7,8,9,11,12,13,14,15,16,17,19,&
		20,22,24,25,26,28,29/)
	integer :: estequio(273,21), charge(273), composicion(273,4), n_iters, includ(273), n_included
	integer :: output_species
	real(kind=8) :: n_atoms_mol(273)
	real(kind=8), allocatable :: which_included(:)
	real(kind=8) :: atomic_partition(7,3,21)
	real(kind=8) :: equilibrium(273), nh, temper, n_e, equilibrium_atomic(2,21)
	real(kind=8), parameter :: NA_ME = -4.89104734208d0, PK_CH = 1.3806503d-16
	real(kind=8) :: x_sol(22), x0_sol(22), P_total, P_elec
end module variables

!*******************************************************************
! Nonlinear system of equations solver
!*******************************************************************
module equations
use variables
implicit none
contains

!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
! The equations are valid for obtaining Pg from Pe and T
! P(22) is the total pressure
! P(1:21) are the pressures of each element
!-------------------------------------------------------------------
	subroutine funcv_Pg_from_T_Pe(x,n,fvec,gmat)
	integer :: n
	real(kind=8) :: x(n), fvec(n), gmat(n,n), components(4)
	integer :: i, j, k, l, m, ind, minim, ind2
	real(kind=8) :: P(n), salida, interna, P_e, minim_ioniz, dinterna(n)
		
		P = x

		fvec = 0.d0
		gmat = 0.d0
		
! We know that P_T = P(H)+P(He)+P(C)+...+Pe
! where P(i)=a(i)*P(H), with a(i) the abundance of element i with respect to H
! Therefore, P_T=K*P(H)+Pe, where K=sum(a(i)). Consequently, P(H) can be isolated and
! be P(H)=(P_T-Pe)/K
! Then, P(i) = a(i) * (P_T-Pe)/K
		x0_sol(1:21) = abund_atom * (P(22) - P_elec) / sum(abund_atom)

		do i = 1, 21
			salida = 0.d0

! If temperature is very high, do not take into account molecular formation
			if (temper < 1.d5) then
				do j = 1, n_included
					k = which_included(j)
					interna = 0.d0
					dinterna = 0.d0
					interna = 1.d0
					
					minim = 0
					minim_ioniz = 100.d0
					
! Product of the partial pressures of the components of molecule k
					do l = 1, 4
						ind = composicion(k,l)

						if (ind /= 0.d0) then						
							if (pot_ion(ind) < minim_ioniz) then
								minim = ind
								minim_ioniz = pot_ion(ind)
							endif

							interna = interna*P(ind)**estequio(k,ind)

! Derivatives
							do m = 1, 4
								ind2 = composicion(k,m)
																
								if (ind2 /= 0.d0) then

									if (ind2 == ind) then

! Include the derivative of P(i) on df(i)/dP(i)
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										else
											dinterna(ind2) = dinterna(ind2) * estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										endif
										
									else

! And then, the rest of partial pressures in the molecule
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = P(ind2)**estequio(k,ind2)
										else
											dinterna(ind) = dinterna(ind) * P(ind2)**estequio(k,ind2)
										endif
										
									endif

								endif
								
							enddo
							
						endif
					enddo


					if (equilibrium(k) == 0.d0) then
						salida = 0.d0
						dinterna = 0.d0
					else
						if (charge(k) == 1) then
							salida = salida + estequio(k,i)*interna / (equilibrium(k) * P_elec ) * &
								equilibrium_atomic(1,minim)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / (equilibrium(k) * P_elec ) * &
								equilibrium_atomic(1,minim)
						else
							salida = salida + estequio(k,i)*interna / equilibrium(k)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / equilibrium(k)
						endif
					endif
																							
				enddo
			endif

! Positive ions partial pressure
			salida = salida + equilibrium_atomic(1,i) * P(i) / P_elec

! Derivatives of positive ions partial pressure
			gmat(i,i) = gmat(i,i) + equilibrium_atomic(1,i) / P_elec

! Negative ions partial pressure
			salida = salida + 1.d0 / equilibrium_atomic(2,i) * P(i) * P_elec

! Derivatives of negative ions partial pressure
			gmat(i,i) = gmat(i,i) + 1.d0 / equilibrium_atomic(2,i) * P_elec
			
! P(i) = Pi + P+ + P- + Pmolec 			
			fvec(i) = P(i) + salida - x0_sol(i)
			
			gmat(i,i) = gmat(i,i) + 1.d0
			gmat(i,22) = gmat(i,22) - abund_atom(i) / sum(abund_atom)

! Contribution to P_e of positive atomic ions
			fvec(22) = fvec(22) + equilibrium_atomic(1,i) * P(i) / P_elec

 			gmat(22,i) = gmat(22,i) + equilibrium_atomic(1,i) / P_elec

! Contribution to P_e of negative atomic ions
			fvec(22) = fvec(22) - 1.d0 / equilibrium_atomic(2,i) * P(i) * P_elec

 			gmat(22,i) = gmat(22,i) - 1.d0 / equilibrium_atomic(2,i) * P_elec
			
		enddo
		
  		fvec(22) = fvec(22) - P_elec
										
	end subroutine funcv_Pg_from_T_Pe

!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
! The equations are valid for obtaining Pe from Pgas and T
! P(22) is the electron pressure
! P(1:21) are the pressures of each element
!-------------------------------------------------------------------
	subroutine funcv_Pe_from_T_Pg(x,n,fvec,gmat)
	integer :: n
	real(kind=8) :: x(n), fvec(n), gmat(n,n)
	integer :: i, j, k, l, m, ind, minim, ind2
	real(kind=8) :: P(n), salida, interna, P_e, minim_ioniz, dinterna(n)
		
		P = x

		fvec = 0.d0
		gmat = 0.d0

! We know that P_T = P(H)+P(He)+P(C)+...+Pe
! where P(i)=a(i)*P(H), with a(i) the abundance of element i with respect to H
! Therefore, P_T=K*P(H)+Pe, where K=sum(a(i)). Consequently, P(H) can be isolated to
! be P(H)=(P_T-Pe)/K
! Then, P(i) = a(i) * (P_T-Pe)/K
		x0_sol(1:21) = abund_atom * (P_total - P(22)) / sum(abund_atom)
		
		do i = 1, 21
			salida = 0.d0

! If temperature is very high, do not take into account molecular formation
			if (temper < 1.d5) then

				do j = 1, n_included
					k = which_included(j)
					interna = 0.d0
					dinterna = 0.d0
					interna = 1.d0
					
					minim = 0
					minim_ioniz = 100.d0
					
! Product of the partial pressures of the components of molecule k
					do l = 1, 4
						ind = composicion(k,l)

						if (ind /= 0.d0) then						
							if (pot_ion(ind) < minim_ioniz) then
								minim = ind
								minim_ioniz = pot_ion(ind)
							endif
							interna = interna*P(ind)**estequio(k,ind)

! Derivatives
							do m = 1, 4
								ind2 = composicion(k,m)

								if (ind2 /= 0.d0) then

									if (ind2 == ind) then

! Include the derivative of P(i) on df(i)/dP(i)
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										else
											dinterna(ind2) = dinterna(ind2) * estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										endif

									else

! And then, the rest of partial pressures in the molecule
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = P(ind2)**estequio(k,ind2)
										else
											dinterna(ind) = dinterna(ind) * P(ind2)**estequio(k,ind2)
										endif

									endif

								endif

							enddo

						endif
					enddo
					
					
					if (equilibrium(k) == 0.d0) then
						salida = 0.d0
					else
						if (charge(k) == 1) then
							salida = salida + estequio(k,i)*interna / (equilibrium(k) * P(22) ) * &
								equilibrium_atomic(1,minim)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / (equilibrium(k) * P(22) ) * &
								equilibrium_atomic(1,minim)
							gmat(i,22) = gmat(i,22) - estequio(k,i)*interna / (equilibrium(k) * P(22)**2 ) * &
								equilibrium_atomic(1,minim)
						else
							salida = salida + estequio(k,i)*interna / equilibrium(k)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / equilibrium(k)
						endif
					endif
																							
				enddo
			endif
			

! Positive ions partial pressure
			salida = salida + equilibrium_atomic(1,i) * P(i) / P(22)

! Derivatives of positive ions partial pressure
			gmat(i,i) = gmat(i,i) + equilibrium_atomic(1,i) / P(22)
			gmat(i,22) = gmat(i,22) - equilibrium_atomic(1,i) * P(i) / P(22)**2
			

! Negative ions partial pressure
			salida = salida + 1.d0 / equilibrium_atomic(2,i) * P(i) * P(22)

! Derivatives of negative ions partial pressure
			gmat(i,i) = gmat(i,i) + 1.d0 / equilibrium_atomic(2,i) * P(22)
			gmat(i,22) = gmat(i,22) - 1.d0 / equilibrium_atomic(2,i) * P(i) / P(22)**2
			
! P(i) = Pi + P+ + P- + Pmolec 			
			fvec(i) = P(i) + salida - x0_sol(i)

			gmat(i,i) = gmat(i,i) + 1.d0
			gmat(i,22) = -abund_atom(i) / sum(abund_atom)

! Contribution to P_e of positive atomic ions
			fvec(22) = fvec(22) + equilibrium_atomic(1,i) * P(i) / P(22)

			gmat(22,i) = gmat(22,i) + equilibrium_atomic(1,i) / P(22)
			gmat(22,22) = gmat(22,22) - equilibrium_atomic(1,i) * P(i) / P(22)**2

! Contribution to P_e of negative atomic ions
			fvec(22) = fvec(22) - 1.d0 / equilibrium_atomic(2,i) * P(i) * P(22)
			
			gmat(22,i) = gmat(22,i) - 1.d0 / equilibrium_atomic(2,i) * P(22)
			gmat(22,22) = gmat(22,22) - 1.d0 / equilibrium_atomic(2,i) * P(i)

		enddo

! Independent term: P_e = all contributors to electrons
		fvec(22) = fvec(22) - P(22)
				
		gmat(22,22) = gmat(22,22) - 1.d0
				
	end subroutine funcv_Pe_from_T_Pg

!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
! The equations are valid for obtaining Pe from Pgas and T
! P(22) is the electron pressure
! P(1:21) are the pressures of each element
!-------------------------------------------------------------------
	subroutine funcv_from_T_Pg_Pe(x,n,fvec,gmat)
	integer :: n
	real(kind=8) :: x(n), fvec(n), gmat(n,n)
	integer :: i, j, k, l, m, ind, minim, ind2
	real(kind=8) :: P(n), salida, interna, minim_ioniz, dinterna(n)

		P = x

		fvec = 0.d0
		gmat = 0.d0

! We know that P_T = P(H)+P(He)+P(C)+...+Pe
! where P(i)=a(i)*P(H), with a(i) the abundance of element i with respect to H
! Therefore, P_T=K*P(H)+Pe, where K=sum(a(i)). Consequently, P(H) can be isolated to
! be P(H)=(P_T-Pe)/K
! Then, P(i) = a(i) * (P_T-Pe)/K
		x0_sol(1:21) = abund_atom * (P_total - P_elec) / sum(abund_atom)

		do i = 1, 21
			salida = 0.d0

! If temperature is very high, do not take into account molecular formation
			if (temper < 1.d5) then

				do j = 1, n_included
					k = which_included(j)
					interna = 0.d0
					dinterna = 0.d0
					interna = 1.d0

					minim = 0
					minim_ioniz = 100.d0

! Product of the partial pressures of the components of molecule k
					do l = 1, 4
						ind = composicion(k,l)

						if (ind /= 0.d0) then
							if (pot_ion(ind) < minim_ioniz) then
								minim = ind
								minim_ioniz = pot_ion(ind)
							endif
							interna = interna*P(ind)**estequio(k,ind)

! Derivatives
							do m = 1, 4
								ind2 = composicion(k,m)

								if (ind2 /= 0.d0) then

									if (ind2 == ind) then

! Include the derivative of P(i) on df(i)/dP(i)
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										else
											dinterna(ind2) = dinterna(ind2) * estequio(k,ind) * P(ind)**(estequio(k,ind)-1.d0)
										endif

									else

! And then, the rest of partial pressures in the molecule
										if (dinterna(ind) == 0.d0) then
											dinterna(ind) = P(ind2)**estequio(k,ind2)
										else
											dinterna(ind) = dinterna(ind) * P(ind2)**estequio(k,ind2)
										endif

									endif

								endif

							enddo

						endif
					enddo


					if (equilibrium(k) == 0.d0) then
						salida = 0.d0
					else
						if (charge(k) == 1) then
							salida = salida + estequio(k,i)*interna / (equilibrium(k) * P_elec ) * &
								equilibrium_atomic(1,minim)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / (equilibrium(k) * P_elec ) * &
								equilibrium_atomic(1,minim)
						else
							salida = salida + estequio(k,i)*interna / equilibrium(k)
							gmat(i,:) = gmat(i,:) + estequio(k,i)*dinterna / equilibrium(k)
						endif
					endif

				enddo
			endif


! Positive ions partial pressure
			salida = salida + equilibrium_atomic(1,i) * P(i) / P_elec

! Derivatives of positive ions partial pressure
			gmat(i,i) = gmat(i,i) + equilibrium_atomic(1,i) / P_elec


! Negative ions partial pressure
			salida = salida + 1.d0 / equilibrium_atomic(2,i) * P(i) * P_elec

! Derivatives of negative ions partial pressure
			gmat(i,i) = gmat(i,i) + 1.d0 / equilibrium_atomic(2,i) * P_elec

! P(i) = Pi + P+ + P- + Pmolec
			fvec(i) = P(i) + salida - x0_sol(i)

			gmat(i,i) = gmat(i,i) + 1.d0			

		enddo

	end subroutine funcv_from_T_Pg_Pe


! ! !-------------------------------------------------------------------
! ! ! Calculates the Jacobian using forward-differences
! ! !-------------------------------------------------------------------	
! 	subroutine fdjac_Pe_from_T_Pg(n,x,fvec,df)
! 	integer :: n,np
! 	real(kind=8) :: df(n,n),fvec(n),x(n),EPS,gmat(n,n)
! 	PARAMETER (EPS=1.d-4)
! 	integer :: i,j
! 	real(kind=8) :: h,temp,f(n)
! 	do j=1,n
! 		temp=x(j)
! 		h=EPS*dabs(temp)
! 		if(h.eq.0.d0)h=EPS
! 		x(j)=temp+h
! 		h=x(j)-temp
! 		call funcv_Pe_from_T_Pg(x,n,f,gmat)
! 		x(j)=temp
! 		do i=1,n
! 			df(i,j)=(f(i)-fvec(i))/h
! 		enddo
! 	enddo
! 
! 	end subroutine fdjac_Pe_from_T_Pg
! ! 
! ! !-------------------------------------------------------------------
! ! ! Calculates the Jacobian using forward-differences
! ! !-------------------------------------------------------------------	
! 	subroutine fdjac_Pg_from_T_Pe(n,x,fvec,df)
! 	integer :: n,np
! 	real(kind=8) :: df(n,n),fvec(n),x(n),EPS,gmat(n,n)
! 	PARAMETER (EPS=1.d-6)
! 	integer :: i,j
! 	real(kind=8) :: h,temp,f(n)
! 	do j=1,n
! 		temp=x(j)
! 		h=EPS*dabs(temp)
! 		if(h.eq.0.d0)h=EPS
! 		x(j)=temp+h
! 		h=x(j)-temp
! 		call funcv_Pg_from_T_Pe(x,n,f,gmat)
! 		x(j)=temp
! 		do i=1,n
! 			df(i,j)=(f(i)-fvec(i))/h
! 		enddo
! 	enddo
! 
! 	end subroutine fdjac_Pg_from_T_Pe

! !-------------------------------------------------------------------
! ! Returns the equations and the Jacobian of the nonlinear set to be solved at point x
! !-------------------------------------------------------------------
! 	subroutine usrfun_Pe_from_T_Pg(x,n,fvec,fjac)
! 	integer :: n, i, j
! 	real(kind=8) :: x(n), fvec(n), fjac(n,n), gmat(n,n)
! 
! 		call funcv_Pe_from_T_Pg(x,n,fvec,gmat)
! 		call fdjac_Pe_from_T_Pg(n,x,fvec,fjac)
! 
! 	end subroutine usrfun_Pe_from_T_Pg
! 	
! !-------------------------------------------------------------------
! ! Returns the equations and the Jacobian of the nonlinear set to be solved at point x
! !-------------------------------------------------------------------
! 	subroutine usrfun_Pg_from_T_Pe(x,n,fvec,fjac)
! 	integer :: n, i, j
! 	real(kind=8) :: x(n), fvec(n), fjac(n,n), gmat(n,n)
! 
! 		call funcv_Pg_from_T_Pe(x,n,fvec,fjac)
! ! 		call fdjac_Pg_from_T_Pe(n,x,fvec,fjac)
! 
! ! 		do i = 1, 22
! ! 			do j = 1, 22
! ! 				print *, fjac(i,j), gmat(i,j)
! ! 			enddo
! ! 			pause
! ! 		enddo
! 
! 	end subroutine usrfun_Pg_from_T_Pe
	
end module equations

module maths_chemical
use equations
implicit none
contains

!-------------------------------------------------------------------
! Solves a system of nonlinear equations using the Newton mthod
!-------------------------------------------------------------------
	subroutine mnewt(which,ntrial,x,n,tolx,tolf)
	integer :: which, n,ntrial
	real(kind=8) :: tolf,tolx,x(n)
	integer :: i,k,indx(n)
	real(kind=8) :: d,errf,errx,fjac(n,n),fvec(n),p(n)
	do k=1,ntrial
	
		if (which == 1) then
			call funcv_Pe_from_T_Pg(x,n,fvec,fjac)   !User subroutine supplies function values at x in fvec
		endif

		if (which == 2) then
			call funcv_Pg_from_T_Pe(x,n,fvec,fjac)   !User subroutine supplies function values at x in fvec
		endif

		if (which == 3) then
			call funcv_from_T_Pg_Pe(x,n,fvec,fjac)   !User subroutine supplies function values at x in fvec
		endif

		errf = 0.d0
		do i=1,n  !Check function convergence.
			errf=errf+dabs(fvec(i))
		enddo 
		if (errf <= tolf) then
			n_iters = k
			return
		endif
		p = -fvec
		
		call ludcmp(fjac,n,n,indx,d) !Solve linear equations using LU decomposition.
		
		call lubksb(fjac,n,n,indx,p)
		
		errx=0.d0  ! Check root convergence.
		x = x + p
		
		do i=1,n   !Update solution.
			errx=errx+dabs(p(i))			
		enddo 
		if(errx <= tolx) then
			n_iters = k
			return
		endif
	enddo 
	if (k == ntrial) then
		x = x0_sol
	endif
	
	end subroutine mnewt
	
!-------------------------------------------------------------------
! Performs the LU decomposition of a matrix
!-------------------------------------------------------------------	
	subroutine ludcmp(a,n,np,indx,d)
	integer :: n,np,indx(n),NMAX
	real(kind=8) :: d,a(np,np),TINY
	PARAMETER (NMAX=500,TINY=1.0d-20) !Largest expected n, and a small number.
	integer :: i,imax,j,k
	real(kind=8) :: aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
	d=1.d0 
	do i=1,n 
		aamax = 0.d0
		do j=1,n
			if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
		enddo 
		if (aamax.eq.0.) then
		  print *, 'singular matrix in ludcmp' 
		  x_sol = x0_sol
		  exit
	   endif
		vv(i)=1./aamax 
	enddo 
	do j=1,n 
		do i=1,j-1 
			sum=a(i,j)
			do k=1,i-1
				sum=sum-a(i,k)*a(k,j)
			enddo 
			a(i,j)=sum
		enddo 
		aamax=0.d0
		do i=j,n 
			sum = a(i,j)
			do k=1,j-1
				sum=sum-a(i,k)*a(k,j)
			enddo 
			a(i,j)=sum
			dum=vv(i)*dabs(sum) 
			if (dum.ge.aamax) then 
				imax=i
				aamax=dum
			endif
		enddo 
		if (j.ne.imax)then 
			do k=1,n 
				dum=a(imax,k)
				a(imax,k)=a(j,k)
				a(j,k)=dum
			enddo 
			d=-d 
			vv(imax)=vv(j) 
		endif
		indx(j)=imax
		if(a(j,j).eq.0.)a(j,j)=TINY
		if(j.ne.n)then 
			dum=1.d0/a(j,j)
			do i=j+1,n
				a(i,j)=a(i,j)*dum
			enddo 
		endif
	enddo 

	end subroutine ludcmp

!-------------------------------------------------------------------
! Solves a system of linear equations using the LU decomposition
!-------------------------------------------------------------------	
	subroutine lubksb(a,n,np,indx,b)
	integer :: n,np,indx(n)
	real(kind=8) :: a(np,np),b(n)

	integer :: i,ii,j,ll
	real(kind=8) :: sum
	ii=0 
	do i=1,n
		ll=indx(i)
		sum=b(ll)
		b(ll)=b(i)
		if (ii.ne.0)then
			do j=ii,i-1
				sum=sum-a(i,j)*b(j)
			enddo 
		else if (sum.ne.0.) then
			ii=i  
		endif
		b(i)=sum
	enddo 
	do i=n,1,-1 
		sum=b(i)
		do j=i+1,n
			sum=sum-a(i,j)*b(j)
		enddo 
		b(i)=sum/a(i,i) 
	enddo 
	end subroutine lubksb
	


end module maths_chemical

!*******************************************************************
! Constants reading module
!*******************************************************************
module lectura
use variables
implicit none

contains

!-------------------------------------------------------------------
! Read the equilibrium constants and the molecular constitution
!-------------------------------------------------------------------
	subroutine read_equil_cte
	integer :: i, j

		open(unit=42,file='TE_data/equil_cte.dat',status='old',action='read')
		read(42,*)
		do i = 1, 273
			read(42,FMT='(A16,5X,F8.5,2X,F10.5,3X,F10.5,5X,F10.5,3X,F10.5,4X,F10.5,4X,F10.5,6X,F9.5,1X,F9.5)')&
				molec(i),(equil(j,i),j=1,9)
		enddo
	
		close(42)
		
	end subroutine read_equil_cte
	
!-------------------------------------------------------------------
! Read the elements included in the molecules
!-------------------------------------------------------------------	
	subroutine read_elements
	integer :: i
		open(unit=42,file='TE_data/elements.dat',status='old',action='read')
		read(42,*) 
		
		do i = 1, 21
			read(42,*) elements(i), abund_atom(i), pot_ion(i), afinidad(i)
! OLD			
			abund_atom(i) = 10.d0**(abund_atom(i)-12.d0)
! Now use the element abundance read in the main program and enhanced as desired
! 			abund_atom(i) = abundance(which_index_element(i))
		enddo
	
		close(42)
	end subroutine read_elements
	
!-------------------------------------------------------------------
! Parse all the molecules with the stequiometric coefficients
!-------------------------------------------------------------------	
	subroutine read_estequio
	integer :: i, j, from, step, found, found_bak, step_valencia, temp_int
	character(len=16) :: temp, temp2
	character(len=1) :: valencia
	character(len=2) :: carga
		charge = 0.d0
		do i = 1, 273
			nombre_mol(i) = ''
			temp2 = ' '
			temp = molec(i)
			found = 1
			found_bak = 1
			do j = 1, 21
				from = index(molec(i),elements(j))				
				if (from /= 0) then															
					
					if (temp(from+1:from+1) == '_') then 
						step = 0
						step_valencia = 2
					else
						step = 1
						step_valencia = 2
					endif
					
					found_bak = found 
					found = found + 2
					valencia = temp(from+step_valencia:from+step_valencia)
					
					read(valencia,*) temp_int
					estequio(i,j) = temp_int
										
					if (valencia == '1') then
						valencia = ''
						found = found - 1
					endif
					temp2(found_bak:found) = temp(from:from+step)//valencia					
											
				endif
			enddo
			
			from = index(temp,'/')
			read (temp(from+1:from+2),*) temp_int
			if (temp_int /= 0) then
				if (temp_int < 10) then
					charge(i) = -temp_int	
				else
					charge(i) = temp_int / 10
				endif
			endif
			
			nombre_mol(i) = temp2
						
		enddo
		
	end subroutine read_estequio
	
!-------------------------------------------------------------------
! Read the elements present in each molecule
!-------------------------------------------------------------------	
	subroutine read_composition
	integer :: i, j, k, temp
				
		n_atoms_mol = 0.d0
		
		do i = 1, 273
			k = 1
			do j = 1, 21
				temp = index(molec(i),elements(j))
				if (temp /= 0) then
					composicion(i,k) = j
					k = k + 1
					n_atoms_mol(i) = n_atoms_mol(i) + estequio(i,j)
				endif
			enddo			
		enddo				
	
	end subroutine read_composition
					
!-------------------------------------------------------------------	
! Calculates the equilibrium constants for a given temperature
!-------------------------------------------------------------------		
	subroutine calc_equil(T)
	real(kind=8) :: T, temp, logar
	integer :: i, j
		
		logar = dlog10(5040.d0 / T)
		do j = 1, 273
			temp = 0.d0
			do i = 0, 8
				temp = temp + equil(i+1,j) * (logar)**i
			enddo			
			
			equilibrium(j) = 10.d0**temp
						
		enddo
		
! Transform the units from SI to cgs multipliying by 10 the necessary times depending on the
! units of the equilibrium constant
		equilibrium = equilibrium * 10.d0 ** (n_atoms_mol - 1.d0)
		
	end subroutine calc_equil

!-------------------------------------------------------------------
! Calculates the equilibrium constants for a given temperature
! Do it only for the molecules that we use
!-------------------------------------------------------------------
	subroutine calc_equil_onlyusedmols(T)
	real(kind=8) :: T, temp, logar
	integer :: i, j, k


		logar = dlog10(5040.d0 / T)
		do j = 1, n_included
			k = which_included(j)
			temp = 0.d0
			do i = 0, 8
				temp = temp + equil(i+1,k) * (logar)**i
			enddo

			equilibrium(k) = 10.d0**temp
		enddo

! Transform the units from SI to cgs multipliying by 10 the necessary times depending on the
! units of the equilibrium constant
		equilibrium = equilibrium * 10.d0 ** (n_atoms_mol - 1.d0)

	end subroutine calc_equil_onlyusedmols
	
!-------------------------------------------------------------------
! Read the equilibrium constants for the atomic and ionic species
!-------------------------------------------------------------------
	subroutine read_partition_cte_atomic
	integer :: i, j, codigo, cod_specie, carga
	real(kind=8) :: data(7)

		open(unit=42,file='TE_data/atomic.dat',status='old',action='read')
		
		atomic_partition = 0.d0
		do i = 1, 59
			read(42,*) codigo, (data(j),j=1,7)
			cod_specie = codigo / 10
			carga = codigo - cod_specie * 10
			atomic_partition(:,carga+1,cod_specie) = data
		enddo
	
		close(42)
		
	end subroutine read_partition_cte_atomic	
	
!-------------------------------------------------------------------	
! Calculates the atomic equilibrium constants for a given temperature
!-------------------------------------------------------------------		
	subroutine calc_equil_atomic(T)
	real(kind=8) :: T, temp, logar, phi(3)
	integer :: i, j, k, l
		
		logar = dlog10(5040.d0 / T)		
	
! All the species
		do l = 1, 21
			
! Calculate the partition function of the neutral, + and - ion of the species			
			do k = 1, 3
				temp = 0.d0
				do i = 0, 6
					temp = temp + atomic_partition(i+1,k,l) * (logar)**i
				enddo		
				phi(k) = 10.d0**temp				
			enddo
									
! Positive ion constant
			equilibrium_atomic(1,l) = 3.41405d0 + NA_ME + 2.5d0 * dlog10(T) - 5039.9d0 * pot_ion(l) / T + &
				dlog10( 2.d0*phi(2)/phi(1) )
			equilibrium_atomic(1,l) = 10.d0**(equilibrium_atomic(1,l))
			
! Negative ion constant			
			equilibrium_atomic(2,l) = 3.41405d0 + NA_ME + 2.5d0 * dlog10(T) - 5039.9d0 * afinidad(l) / T + &
				dlog10( 2.d0*phi(1)/phi(3) )		
			equilibrium_atomic(2,l) = 10.d0**(equilibrium_atomic(2,l))
						
		enddo
		
! Transform the units from N/m^2 to dyn/cm^2
		equilibrium_atomic = equilibrium_atomic * 10.d0
		
	end subroutine calc_equil_atomic
	
!-------------------------------------------------------------------	
! Read what species are going to be included
!-------------------------------------------------------------------			
	subroutine read_what_species
	integer :: i, j
	character(len=16) :: name
		
		includ = 0
		
		open(unit=43,file='TE_data/species.dat',status='old',action='read')
		
		read(43,*) n_included
		allocate(which_included(n_included))
		do i = 1, n_included
			read(43,FMT='(A16)') name
					
			do j = 1, 273
				if (molec(j) == name) then
					includ(j) = 1
					which_included(i) = j
!					print *, 'Including molecule number ',i,' : ',molec(j), ' with index ',j
				endif				
				
			enddo
			
		enddo
		
		close(43)
	end subroutine read_what_species
	
!-------------------------------------------------------------------	
! Write the output model with the corresponding molecular abundance
!-------------------------------------------------------------------
	subroutine sacar(which)
	integer :: n_depths, j, k, which, i
	real(kind=8) :: molecu(273), atomic(21), height
		!open (unit=12,file=output_file,status='old',action='read')
		open (unit=43,file=output_model,status='replace',action='write')
        
		read(42,*) n_depths
        
		do j = 1, n_depths              
			read(42,*) height, (atomic(i),i=1,21), (molecu(k), k=1,273)               
        
			if (which < 0) then
				write(43,*) height, atomic(-which)
			else
				write(43,*) height, molecu(which)
			endif
		enddo
        
!		print *, 'File ',output_model,' created with the abundance of molecule number ',which
        
		close(43)
!		close(12)
	
	end subroutine sacar

end module lectura


!*******************************************************************
!*******************************************************************
! Calculates the molecular and atomic abundances in thermodynamic 
! equilibrium
!*******************************************************************
!*******************************************************************

module chemical_eq
use variables
use lectura
use maths_chemical
implicit none
contains

!-----------------------------------------------------------------
! Calculates the abundance of a given molecule in chemical equilibrium
! INPUT :
!	mol_code : integer to identify the molecule (see file TE_data/equil_cte.dat)
!	n_grid : number of grid points in the T, nH and ne arrays
!	height_in : vector of heights (it is used for nothing indeed)
!	temper_in : array of temperatures in K
!  PT_in : array of total pressure in dyn*cm^2
!  abundance_file : file in which the molecular abundances will be saved
! OUTPUT : 
!  PH_out : partial pressure of H atoms (dyn/cm^2)
!  PHminus_out : partial pressure of H- atoms (dyn/cm^2)
!  PHplus_out : partial pressure of H+ atoms (dyn/cm^2)
!  PH2_out : partial pressure of H2 molecules (dyn/cm^2)
!  PH2plus_out : partial pressure of H2+ molecules (dyn/cm^2)
!-----------------------------------------------------------------
	function calculate_abundance_Pe_from_T_Pg(mol_code, n_grid, height_in, temper_in, PT_in, &
		abundance_file, PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out, P_elec)
	integer :: n_grid, mol_code
	real(kind=8) :: height_in(n_grid), temper_in(n_grid), PT_in(n_grid), n_e_in(n_grid), initial_value
	real(kind=8) :: calculate_abundance_Pe_from_T_Pg(n_grid), abun_out(n_grid)
	real(kind=8), dimension(n_grid) :: PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out, P_elec
	character(len=40) :: abundance_file
	real(kind=8) :: mole(273), height, minim_ioniz
	integer :: i, ind, l, loop	, minim, k
	
		open(unit=45, file=abundance_file,action='write',status='replace')
		write(45,*) n_grid

! Reading equilibrium constants of the 273 molecules included
		call read_equil_cte
! Reading equilibrium constants of atomic and ionic species
		call read_partition_cte_atomic
! Reading 21 elements
		call read_elements
! Reading estequiometer values
		call read_estequio
! Reading composition of the 273 molecules
		call read_composition
! Reading what molecules are included
		call read_what_species
		
		do loop = 1, n_grid
		
			height = height_in(loop)
			temper = temper_in(loop)
			P_total = PT_in(loop)

! Calculating molecular equilibrium constants
			call calc_equil_onlyusedmols(temper)
! Calculating ionic equilibrium constants
			call calc_equil_atomic(temper)		

! Initial conditions
! Initialize assuming that the H pressure equals the total gas pressure			
			x_sol(1:21) = P_total * abund_atom
			initial_value = 1
			x_sol(22) = initial_value*P_total

			call mnewt(1,100,x_sol,22,1.d-4,1.d-4)
			
			do while (minval(x_sol(1:22)) < 0.d0)
				print *, 'SOLVING AGAIN...'
				initial_value = initial_value / 10.d0
				x_sol(22) = initial_value*P_total
				call mnewt(1,100,x_sol,22,1.d-4,1.d-4)
			enddo
									
			P_elec(loop) = x_sol(22)

	 		mole = 1.d0					
			i = mol_code
			minim = 0
			minim_ioniz = 100.d0
			do l = 1, 4
				ind = composicion(i,l)
				if (ind /= 0.d0) then						
					if (pot_ion(ind) < minim_ioniz) then
						minim = ind
						minim_ioniz = pot_ion(ind)
					endif
					mole(i) = mole(i) * x_sol(ind)**estequio(i,ind)
				endif
			enddo	
			if (equilibrium(i) == 0.d0) then
				mole(i) = 0.d0
			else
				if (charge(i) == 1) then
					mole(i) = (mole(i) / (equilibrium(i) * P_elec(loop)) * &
						equilibrium_atomic(1,minim)) / (PK_CH * temper)
				else
					mole(i) = (mole(i) / equilibrium(i)) / (PK_CH * temper)
				endif
			endif
			
			if (.not.(mole(i)>0) .and. .not.(mole(i)<=0) ) mole(i) = 0.d0
			write(45,*) height, mole(i)
			abun_out(loop) = mole(i)
			
! Now extract also the partial pressure from H, H-, H+, H2 and H2+			
			
			PH_out(loop) = x_sol(1)			
			PHminus_out(loop) = x_sol(1) * P_elec(loop) / equilibrium_atomic(2,1)
			PHplus_out(loop) = x_sol(1) / P_elec(loop) * equilibrium_atomic(1,1)
			if (temper < 1.d5) then
				PH2_out(loop) = x_sol(1)**2 / equilibrium(1)
				PH2plus_out(loop) = x_sol(1)**2 / P_elec(loop) * (equilibrium_atomic(1,1) / equilibrium(2))
			else
				PH2_out(loop) = 0.d0
				PH2plus_out(loop) = 0.d0
			endif

		enddo

		close(45)
		print *, 'Calculated abundance for molecule : ', molec(mol_code)		
		calculate_abundance_Pe_from_T_Pg = abun_out
		
		deallocate(which_included)

	end function calculate_abundance_Pe_from_T_Pg
	
!-----------------------------------------------------------------
! Calculates the abundance of a given molecule in chemical equilibrium
! INPUT :
!	mol_code : integer to identify the molecule (see file TE_data/equil_cte.dat)
!	n_grid : number of grid points in the T, nH and ne arrays
!	height_in : vector of heights (it is used for nothing indeed)
!	temper_in : array of temperatures in K
!  PT_in : array of total pressure in dyn*cm^2
!  abundance_file : file in which the molecular abundances will be saved
! OUTPUT : 
!  PH_out : partial pressure of H atoms (dyn/cm^2)
!  PHminus_out : partial pressure of H- atoms (dyn/cm^2)
!  PHplus_out : partial pressure of H+ atoms (dyn/cm^2)
!  PH2_out : partial pressure of H2 molecules (dyn/cm^2)
!  PH2plus_out : partial pressure of H2+ molecules (dyn/cm^2)
!-----------------------------------------------------------------
	function calculate_abundance_Pg_from_T_Pe(mol_code, n_grid, height_in, temper_in, Pe_in, &
		abundance_file, PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out, P_total)
	integer :: n_grid, mol_code
	real(kind=8) :: height_in(n_grid), temper_in(n_grid), Pe_in(n_grid), n_e_in(n_grid)
	real(kind=8) :: calculate_abundance_Pg_from_T_Pe(n_grid), abun_out(n_grid), initial_values(22)
	real(kind=8), dimension(n_grid) :: PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out, P_total
	character(len=40) :: abundance_file
	real(kind=8) :: mole(273), height, minim_ioniz
	integer :: i, ind, l, loop	, minim
	
		open(unit=45, file=abundance_file,action='write',status='replace')
		write(45,*) n_grid

! Reading equilibrium constants of the 273 molecules included
		call read_equil_cte
! Reading equilibrium constants of atomic and ionic species
		call read_partition_cte_atomic
! Reading 21 elements
		call read_elements
! Reading estequiometer values
		call read_estequio
! Reading composition of the 273 molecules
		call read_composition
! Reading what molecules are included
		call read_what_species
		
		do loop = 1, n_grid
		
			height = height_in(loop)
			temper = temper_in(loop)
			P_elec = Pe_in(loop)

! Calculating molecular equilibrium constants
			call calc_equil_onlyusedmols(temper)
! Calculating ionic equilibrium constants
			call calc_equil_atomic(temper)		

! Initial conditions
! Initialize assuming that the total gas pressure is given by H
			initial_values = 1.d-3
 			x_sol(1:22) = initial_values !P_elec * abund_atom * 1.d0

			call mnewt(2,100,x_sol(1:22),22,1.d-6,1.d-6)

			do while (minval(x_sol(1:22)) < 0.d0)
				print *, 'SOLVING AGAIN...'
				initial_values = initial_values / 10.d0
				x_sol(1:22) = initial_values
				call mnewt(2,100,x_sol(1:22),22,1.d-6,1.d-6)
			enddo

			P_total(loop) = x_sol(22)
		
	 		mole = 1.d0
			i = mol_code
			minim = 0
			minim_ioniz = 100.d0
			do l = 1, 4
				ind = composicion(i,l)
				if (ind /= 0.d0) then						
					if (pot_ion(ind) < minim_ioniz) then
						minim = ind
						minim_ioniz = pot_ion(ind)
					endif
					mole(i) = mole(i) * x_sol(ind)**estequio(i,ind)
				endif
			enddo	
			if (equilibrium(i) == 0.d0) then
				mole(i) = 0.d0
			else
				if (charge(i) == 1) then
					P_elec = n_e * PK_CH * temper
					mole(i) = (mole(i) / (equilibrium(i) * P_elec) * &
						equilibrium_atomic(1,minim)) / (PK_CH * temper)
				else
					mole(i) = (mole(i) / equilibrium(i)) / (PK_CH * temper)
				endif
			endif
			
			if (.not.(mole(i)>0) .and. .not.(mole(i)<=0) ) mole(i) = 0.d0
			write(45,*) height, mole(i)
			abun_out(loop) = mole(i)
			
! Now extract also the partial pressure from H, H-, H+, H2 and H2+
			PH_out(loop) = x_sol(1)			
			PHminus_out(loop) = x_sol(1) * P_elec / equilibrium_atomic(2,1)
			PHplus_out(loop) = x_sol(1) / P_elec * equilibrium_atomic(1,1)
			if (temper < 1.d5) then
				PH2_out(loop) = x_sol(1)**2 / equilibrium(1)
				PH2plus_out(loop) = x_sol(1)**2 / P_elec * (equilibrium_atomic(1,1) / equilibrium(2))
			else
				PH2_out(loop) = 0.d0
				PH2plus_out(loop) = 0.d0
			endif

		enddo

		close(45)
		print *, 'Calculated abundance for molecule : ', molec(mol_code)		
		calculate_abundance_Pg_from_T_Pe = abun_out
		
		deallocate(which_included)

	end function calculate_abundance_Pg_from_T_Pe


!-----------------------------------------------------------------
! Calculates the abundance of a given molecule in chemical equilibrium
! INPUT :
!	mol_code : integer to identify the molecule (see file TE_data/equil_cte.dat)
!	n_grid : number of grid points in the T, nH and ne arrays
!	height_in : vector of heights (it is used for nothing indeed)
!	temper_in : array of temperatures in K
!  PT_in : array of total pressure in dyn*cm^2
!  abundance_file : file in which the molecular abundances will be saved
! OUTPUT :
!  PH_out : partial pressure of H atoms (dyn/cm^2)
!  PHminus_out : partial pressure of H- atoms (dyn/cm^2)
!  PHplus_out : partial pressure of H+ atoms (dyn/cm^2)
!  PH2_out : partial pressure of H2 molecules (dyn/cm^2)
!  PH2plus_out : partial pressure of H2+ molecules (dyn/cm^2)
!-----------------------------------------------------------------
	function calculate_abundance_from_T_Pe_Pg(mol_code, n_grid, height_in, temper_in, Pe_in, Pt_in,&
		abundance_file, PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out)
	integer :: n_grid, mol_code
	real(kind=8) :: height_in(n_grid), temper_in(n_grid), Pe_in(n_grid), n_e_in(n_grid), Pt_in(n_grid)
	real(kind=8) :: calculate_abundance_from_T_Pe_Pg(n_grid), abun_out(n_grid), initial_values(21)
	real(kind=8), dimension(n_grid) :: PH_out, PHminus_out, PHplus_out, PH2_out, PH2plus_out
	character(len=40) :: abundance_file
	real(kind=8) :: mole(273), height, minim_ioniz
	integer :: i, ind, l, loop	, minim

		open(unit=45, file=abundance_file,action='write',status='replace')
		write(45,*) n_grid

! Reading equilibrium constants of the 273 molecules included
		call read_equil_cte
! Reading equilibrium constants of atomic and ionic species
		call read_partition_cte_atomic
! Reading 21 elements
		call read_elements
! Reading estequiometer values
		call read_estequio
! Reading composition of the 273 molecules
		call read_composition
! Reading what molecules are included
		call read_what_species

		do loop = 1, n_grid

			height = height_in(loop)
			temper = temper_in(loop)
			P_elec = Pe_in(loop)
			P_total = Pt_in(loop)

! Calculating molecular equilibrium constants
! 			call calc_equil(temper)
			call calc_equil_onlyusedmols(temper)
! Calculating ionic equilibrium constants
			call calc_equil_atomic(temper)

! Initial conditions
! Initialize assuming that the total gas pressure is given by H
			initial_values = 1.d-3
 			x_sol(1:21) = initial_values !P_elec * abund_atom * 1.d0

			call mnewt(3,100,x_sol(1:21),21,1.d-6,1.d-6)

			do while (minval(x_sol(1:21)) < 0.d0)
				print *, 'SOLVING AGAIN...'
				initial_values = initial_values / 10.d0
				x_sol(1:21) = initial_values
				call mnewt(3,100,x_sol(1:21),21,1.d-6,1.d-6)
			enddo

	 		mole = 1.d0
			i = mol_code
			minim = 0
			minim_ioniz = 100.d0
			do l = 1, 4
				ind = composicion(i,l)
				if (ind /= 0.d0) then
					if (pot_ion(ind) < minim_ioniz) then
						minim = ind
						minim_ioniz = pot_ion(ind)
					endif
					mole(i) = mole(i) * x_sol(ind)**estequio(i,ind)
				endif
			enddo
			if (equilibrium(i) == 0.d0) then
				mole(i) = 0.d0
			else
				if (charge(i) == 1) then
					P_elec = n_e * PK_CH * temper
					mole(i) = (mole(i) / (equilibrium(i) * P_elec) * &
						equilibrium_atomic(1,minim)) / (PK_CH * temper)
				else
					mole(i) = (mole(i) / equilibrium(i)) / (PK_CH * temper)
				endif
			endif

			if (.not.(mole(i)>0) .and. .not.(mole(i)<=0) ) mole(i) = 0.d0
			write(45,*) height, mole(i)
			abun_out(loop) = mole(i)

! Now extract also the partial pressure from H, H-, H+, H2 and H2+
			PH_out(loop) = x_sol(1)
			PHminus_out(loop) = x_sol(1) * P_elec / equilibrium_atomic(2,1)
			PHplus_out(loop) = x_sol(1) / P_elec * equilibrium_atomic(1,1)
			if (temper < 1.d5) then
				PH2_out(loop) = x_sol(1)**2 / equilibrium(1)
				PH2plus_out(loop) = x_sol(1)**2 / P_elec * (equilibrium_atomic(1,1) / equilibrium(2))
			else
				PH2_out(loop) = 0.d0
				PH2plus_out(loop) = 0.d0
			endif

		enddo

		close(45)
		print *, 'Calculated abundance for molecule : ', molec(mol_code)
		calculate_abundance_from_T_Pe_Pg = abun_out

		deallocate(which_included)

	end function calculate_abundance_from_T_Pe_Pg


end module chemical_eq