module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline
    real(dp)     :: dx, S_func, j_func, z, eta0, x0, x_min, x_max, d, e, f, dk
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand, etaa, transfer1, transfer2, transfer3, transfer4, transfer5, transfer6 
    real(dp), 	  allocatable, dimension(:)	  :: spectrum1, spectrum2, spectrum3, spectrum4, spectrum5, spectrum6
    real(dp),     allocatable, dimension(:,:)     :: j_l, j_l2
    real(dp),     allocatable, dimension(:)       :: x_arg, int_arg, cls, cls2
    real(dp),     pointer, dimension(:)           :: k, x
    real(dp),     pointer, dimension(:,:,:,:)     :: S_coeff
    real(dp),     pointer, dimension(:,:)         :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2, final
    real(dp),     allocatable, dimension(:)       :: x_hires, k_hires

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    call get_hires_source_function(k, x, S)
    allocate(k_hires(5000))
    allocate(x_hires(5000))
    k_hires = k
    x_hires = x

    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    allocate(j_l_spline(n_spline))
    allocate(j_l_spline2(n_spline))
    allocate(integrand(l_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    ! arrays for printing out integrands to file for 6 different l values
    allocate(transfer1(5000))
    allocate(transfer2(5000))
    allocate(transfer3(5000))
    allocate(transfer4(5000))
    allocate(transfer5(5000))
    allocate(transfer6(5000))
    allocate(spectrum1(5000))
    allocate(spectrum2(5000))
    allocate(spectrum3(5000))
    allocate(spectrum4(5000))
    allocate(spectrum5(5000))
    allocate(spectrum6(5000))

    ! stores the Bessels to file so only computed once
    inquire(file='clmod.unf', exist=exist)
    if (exist) then
       open(96,file='clmod.unf', form='unformatted')
       read(96) j_l
       read(96) j_l2
       read(96) z_spline
       close(96)
    else

    ! Initializes spherical Bessel functions for each l; use 5400 sampled points between 
    ! z = 0 and 3500. Each function is then splined
    dx = (3500.d0 - 0.1d0) / real(n_spline - 1.d0)
    do j = 1, l_num
       do i = 1, n_spline
    	  z_spline(i) = 0.1 + (i-1) * dx
          call sphbes(ls(j), z_spline(i), j_l(i,j))
       end do 
	
       call spline(z_spline, j_l(:,j), 1.d30, 1.d30, j_l2(:,j))
    end do 

      open(96,file='clmod.unf', form='unformatted')
       write(96) j_l
       write(96) j_l2
       write(96) z_spline
       close(96)
   end if 

    ! Computes the C_l's for each given l
    ! pre-compute eta grid and z to save calculation time
    allocate(etaa(5000))
    eta0 = get_eta(x_hires(5000))
    do i = 1, 5000
	etaa(i) = eta0 - get_eta(x_hires(i))
    end do 

    open(96, file = 'bes.dat')
    do l= 60, 120
      call sphbes(l,100.d0,x_min)
      write(96,*) l, x_min*x_min
    end do
    close(96)

    integrand(:) = 0.d0
    cls(:) = 0.d0
    do l = 1, l_num !l loop

	do j = 1, 5000 !k loop
	   !reset integrand array after each integration over x 
	   integrand(l) = 0.d0
	   do i = 1, 5000 !x loop  
	      ! set up integrand of transfer function, Theta_l(k) for calculation over inner loop of x   
	      dx = 1.d0 / 4999.d0  
	      z = k_hires(j) * etaa(i)
	      j_func = splint(z_spline, j_l(:,l), j_l2(:,l), z) 
	      S_func = S(j,i)
	      integrand(l) = integrand(l) + S_func * j_func * dx

            if(l==1 .and. j==1700)then
	    open(96, file = 'integrand.dat')
     	    write(96,*) x_hires(i),  S_func*j_func
	    endif 

	   end do 

	   ! Integrates P(k) * (Theta_l^2 / k) looped over k to find un-normalized C_l's
	   dk = 1.d0 / 4999.d0
	   d = n_s - 1.d0
	   e = (c * k_hires(j) / H_0) ** d
	   f = integrand(l) * integrand(l)
	   cls(l) = cls(l) + e * f * dk / k_hires(j)
	
           !store integrand arrays for ls(l) = (12, 68, 260, 510, 830, 1150)
	   if(l==7)then
	   transfer1(j) = integrand(l)
           spectrum1(j) = f / k_hires(j)
	   endif

	   if(l==14)then
	   transfer2(j) = integrand(l)
           spectrum2(j) = f / k_hires(j)
	   endif
	  
	   if(l==25)then
	   transfer3(j) = integrand(l)
           spectrum3(j) = f / k_hires(j)
	   endif

	   if(l==31)then
	   transfer4(j) = integrand(l)
           spectrum4(j) = f / k_hires(j)
	   endif

	   if(l==37)then
	   transfer5(j) = integrand(l)
           spectrum5(j) = f / k_hires(j)
	   endif

	   if(l==43)then
	   transfer6(j) = integrand(l)
           spectrum6(j) = f / k_hires(j)
	   endif

	end do
	! integrated Cl stored for each l in ls. 
        cls(l) = cls(l) * real(ls(l)*(ls(l) + 1.d0), dp) / (2.d0* pi)
        write(*,*) cls(l)
    end do
 close(96)

    open(96, file = 'transferl12.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer1(j)
    end do
    close(96)

    open(96, file = 'transferl70.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer2(j)
    end do
    close(96)

    open(96, file = 'transferl275.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer3(j)
    end do
    close(96)

    open(96, file = 'transferl550.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer4(j)
    end do
    close(96)

    open(96, file = 'transferl850.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer5(j)
    end do
    close(96)

    open(96, file = 'transferl1150.dat')
    do j = 1, 5000
       write(96,*) c*k_hires(j)/H_0,  transfer6(j)
    end do
    close(96)
   
    open(96, file = 'spectintegrandl12.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  real(ls(7)*(ls(7)+1.d0),dp)*spectrum1(j)*H_0/c
    end do
    close(96)

    open(96, file = 'spectintegrandl70.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  4970.d0*spectrum2(j)*H_0/c
    end do
    close(96)

    open(96, file = 'spectintegrandl275.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  75900.d0*spectrum3(j)*H_0/c
    end do
    close(96)

    open(96, file = 'spectintegrandl550.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  303050.d0*spectrum4(j)*H_0/c
    end do
    close(96)

    open(96, file = 'spectintegrandl850.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  723350.d0*spectrum5(j)*H_0/c
    end do
    close(96)

    open(96, file = 'spectintegrandl1150.dat')
    do j = 1,5000
       write(96,*) c*k_hires(j)/H_0,  1323650.d0*spectrum6(j)*H_0/c
    end do
    close(96)
	    

    ! Spline C_l's found above, and output smooth C_l curve for each integer l
    allocate(final(1200))
    
    call spline(real(ls, dp), cls, 1.d30, 1.d30, cls2)
    
    ! writes splined power spectrum to file 
    open(96, file = 'powspec.dat')
    do l= 1, 1200
      final(l) = splint(real(ls, dp), cls, cls2, real(l, dp))
      write(96,*) l, 5775.d0 * final(l) / splint(real(ls, dp), cls, cls2, 225.d0)
    end do
    close(96)

  end subroutine compute_cls
  
end module cl_mod
