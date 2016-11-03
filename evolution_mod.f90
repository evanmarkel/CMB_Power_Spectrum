module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  use spline_1D_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init      = 1.d-8
  real(dp),     parameter, private :: k_min       = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max       = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k         = 100
  integer(i4b), parameter          :: n_brec      = 500  ! number of x-grid points before recombination
  integer(i4b), parameter          :: n_xtot      = 1000 ! total x-grid points to be used for parameter integration
  integer(i4b), parameter, private :: lmax_int    = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list and time array for x_init until recombination
  real(dp), allocatable, dimension(:) :: ks
  real(dp), allocatable, dimension(:) :: x_brec
  real(dp), allocatable, dimension(:) :: x_tot

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

contains

  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j, n1, n2, n_hitot
    real(dp)     :: a, g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp)	 :: dx, x_start_rec, x_end_rec, x_0, one, two, three, four
    real(dp), allocatable, dimension(:,:) :: S_lores
    
    ! splined coeffiecients
    real(dp), allocatable, dimension(:,:,:,:) :: coeff

    ! allocate low res source function 
    allocate(S_lores(n_k, n_t))
    allocate(coeff(4,4,n_k,n_t))

    do i = 0, n_t-1

    ! get values for k-independent parameters
    a      = exp(x_t(i+1))
    g      = get_g(x_t(i+1))
    dg     = get_dg(x_t(i+1))
    ddg    = get_ddg(x_t(i+1))
    tau    = get_tau(x_t(i+1))
    dt     = get_dtau(x_t(i+1))
    ddt    = get_ddtau(x_t(i+1))
    H_p    = get_H_p(x_t(i+1))
    dH_p   = get_dH_p(x_t(i+1))
    ddHH_p = 5.d-1 * H_0 * H_0 * ((Omega_b + Omega_m) * a ** (-1) + 4.d0 * (Omega_r) * a ** (-2) + &
             & 4.d0 * (Omega_lambda) * a**2)

       do j = 0, n_k-1
       ! get values for k, x dependent parameters
       Pi   = Theta(n_brec+i-1,2,j)
       dPi  = dTheta(i,2,j)
       ddPi = 4.d-1 * c*ks(j+1) / H_p * (-dH_p / H_p * Theta(n_brec+i-1,1,j) + dTheta(i,1,j)) + &
	      & 3.d-1 * (ddt * Pi + dt * dPi) - (6.d-1 * c*ks(j+1)/ H_p) * & 
	      & ((-dH_p/H_p)*Theta(n_brec+i-1,3,j) + dTheta(i,3,j))

       ! Compute Source function for set k, x grids
       one   = g * (Theta(n_brec+i-1,0,j) + Psi(n_brec+i-1,j) + Pi / 4.d0)
       two   = exp(-tau) * (dPsi(i,j) - dPhi(i,j))
       three = -(1.d0 / (c*ks(j+1))) * (g * v_b(n_brec+i-1,j) * dH_p + g * H_p * dv_b(i,j) + H_p * v_b(n_brec+i-1,j) * dg)
       four  = 75.d-2 / (c*c*ks(j+1)*ks(j+1)) * (ddHH_p * g * Pi + 3.d0 * H_p * dH_p * (dg * Pi + g * dPi) + H_p * H_p * (ddg * Pi + 2.d0 * dg * dPi + g * ddPi))

       S_lores(j+1,i+1) = one + two + three + four

       end do 
    end do 

    ! 2D spline the source function over existing k,x grids
    call splie2_full_precomp(ks, x_t, S_lores, coeff)

    ! Variables from time_mod needed to make high_res k,x grids
    x_start_rec = -log(1631.4d0)  ! x of start of recombination
    x_end_rec   = -log(615.2d0)   ! x of end of recombination
    n1          = 2000		  ! grid points during recombination
    n2	        = 3000		  ! grid points after recombination til present
    n_hitot     = n1 + n2	  ! total high resolution grid point total
    x_0         = 0.d0            ! x today

    ! Allocate x,k arrays for high resolution
    allocate(x(n_hitot))
    allocate(k(n_hitot))
    allocate(S(n_hitot,n_hitot))

    ! x values in recombination era. 200 grid points in the range z = [1630.4, 614].
    dx = (x_end_rec - x_start_rec) / (n1 - 1)
    do i = 1 , n1
       x(i) = x_start_rec + (i - 1) * dx
    end do

    ! x va lues in post-recombination era. 300 grid points in the range z = (614, 0]. 
    dx = (x_0 - x_end_rec) / (n2)
    do i = 1 , n2
	x(n1 + i) = x_end_rec + i * dx
    end do

    ! Initializes k-grid, ks; quadratic between k_min and k_max
 
    do i = 1 , n_hitot
    k(i) = k_min + (k_max - k_min) * real(i - 1.d0) / real(n_hitot - 1.d0)
    end do

    ! Outputs a pre-computed high-res 2D array (over k and x) for the
    ! source function, S(k,x). 
    do j = 1, n_hitot
 	do i = 1, n_hitot
		S(j,i) = splin2_full_precomp(ks, x_t, coeff, k(j),x(i))	 
	end do 
    end do

  end subroutine get_hires_source_function


  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i
    real(dp)     :: x_init, x_rc, H_p, dtau, epsil, dx_brec

    ! initialize x=log(a_initial) and beginning of recombination to calculate IC values of dtau and H_p
    x_init = log(a_init)
    x_rc   = x_t(1)
    H_p    = get_H_p(x_init)
    dtau   = get_dtau(x_init)

    allocate(ks(n_k))

    ! Allocate x_brec array which will store 500 x values from [x_init, x_start_rec]
    allocate (x_brec(n_brec))
    allocate (x_tot(n_xtot-1))

    ! Initializes k-grid, ks; quadratic between k_min and k_max
    do i = 1 , n_k
    ks(i) = k_min + (k_max - k_min) * (real(i-1.d0) / real(n_k - 1.d0))**2
    end do

    ! x values before recombination
    dx_brec = (x_rc - x_init) / (n_brec - 1)
    do i = 1 , n_brec
       x_brec(i) = x_init + (i - 1) * dx_brec
    end do

    ! combine two arrays into total x-array
    x_tot(1:n_brec) = x_brec(1:n_brec)
    x_tot(n_brec:n_xtot) = x_t(1:n_t)

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_xtot-1, 0:lmax_int, 0:n_k-1))
    allocate(delta(0:n_xtot-1, 0:n_k-1))
    allocate(delta_b(0:n_xtot-1, 0:n_k-1))
    allocate(v(0:n_xtot-1, 0:n_k-1))
    allocate(v_b(0:n_xtot-1, 0:n_k-1))
    allocate(Phi(0:n_xtot-1, 0:n_k-1))
    allocate(Psi(0:n_xtot-1, 0:n_k-1))
    allocate(dPhi(0:n_t-1, 0:n_k-1))
    allocate(dPsi(0:n_t-1, 0:n_k-1))
    allocate(dv_b(0:n_t-1, 0:n_k-1))
    allocate(dTheta(0:n_t-1, 0:lmax_int, 0:n_k-1))

    ! Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,0:n_k-1)     = 1.d0
    delta(0,0:n_k-1)   = 1.5d0*Phi(0,0:n_k-1)
    delta_b(0,0:n_k-1) = 1.5d0*Phi(0,0:n_k-1)

    do i = 0, n_k-1
    epsil = (c * ks(i+1)) / (H_p * dtau)

       v(0,i)       = (c * ks(i+1) * Phi(0,i)) / (2.d0 * H_p)
       v_b(0,i)     = v(0,i)
       Theta(0,0,i) = 0.5d0 * Phi(0,i)
       Theta(0,1,i) = -(c * ks(i+1) * Phi(0,i)) / (6.d0 * H_p)
       Theta(0,2,i) = -(20.d0 * c * ks(i+1) * Theta(0,1,i)) / (45.d0 * H_p * dtau)
       do l = 3, lmax_int
      Theta(0,l,i) = -(real(l) /(2.d0*real(l) + 1.d0)) * epsil * Theta(0,l-1,i)
       end do
       Psi(0,i)     = -Phi(0,i) - 12.d0*(H_0 / (c * ks(i+1)* a_init))**2 * Omega_r * Theta(0,2,i)
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l, i_tc_end
    real(dp)     :: xb1, xb2, x1, x2, x_init, x_tc, epsil
    real(dp)     :: eps, hmin, h1, H_p, dtau, t1, t2, a_current, R
    logical(lgt) :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    x_init = log(a_init)
    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

!    inquire(file='cmbdata.unf', exist=exist)
 !   if (exist) then
  !     open(58,file='cmbdata.unf', form='unformatted')
   !    read(58) delta
    !   read(58) delta_b
     !  read(58) v
      ! read(58) v_b
 !      read(58) Theta
  !     read(58) Phi
   !    read(58) Psi 
    !   read(58) dPhi
     !  read(58) dPsi
      ! read(58) dTheta
       !read(58) dv_b
     !  close(58)
    !else

    ! Propagate each k-mode independently
    do k = 0, n_k-1

       k_current = ks(k+1)  ! Store k_current as a global module variable
       h1        = 1.d-5

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed,
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
       if (abs(x_tc) - abs(x_t(1)) < 1.d-3) then
          x_tc = x_t(1)
       end if

    ! Loop over x-before recombination grid and integrate tight coupling
    ! equations and regular regime equations between x_tc and x_start_rec
    ! for higher k-modes 
    do j = 1, n_brec-1
    xb1 = x_brec(j)
    xb2 = x_brec(j+1)
    a_current = exp(xb2)

         ! Integrates from x_init until the end of tight coupling, using
         ! the tight coupling equations
         call odeint(y_tight_coupling, xb1, xb2, eps, h1, hmin, derivs_tight_coupling, bsstep, output_tight_coupling)

         ! Defines H_prime and dtau during tight coupling for global variable calcs
         H_p = get_H_p(xb2)
         dtau  = get_dtau(xb2)
         epsil = (c * k_current) / (H_p * dtau)
  
         ! Store tight coupling values in global parameter arrays
         delta(j,k)   = y_tight_coupling(1)
         delta_b(j,k) = y_tight_coupling(2)
         v(j,k)       = y_tight_coupling(3)
         v_b(j,k)     = y_tight_coupling(4)
         Phi(j,k)     = y_tight_coupling(5)
         Theta(j,0,k) = y_tight_coupling(6)
         Theta(j,1,k) = y_tight_coupling(7)
         Theta(j,2,k) = -(20.d0 * c * k_current * Theta(j,1,k)) / (45.d0 * H_p * dtau)
         do l = 3, lmax_int
            Theta(j,l,k) = - (real(l) /(2.d0*real(l) + 1.d0)) * epsil * Theta(j,l-1,k)
         end do
         Psi(j,k)     = -Phi(j,k) - 12.d0 *(H_0 / ( c * k_current * a_current))**2 * &
			&Omega_r * Theta(j,2,k)

         if (xb2 >= x_tc) then
         i_tc_end = j+1
         ! Sets up variables for integration from the end of tight coupling
         ! until today
         y(1:7) = y_tight_coupling(1:7)
         y(8)   = -(20.d0 * c * k_current * y(7)) / (45.d0 * H_p * dtau)
         do l = 3, lmax_int
            y(6+l) = - (real(l) /(2.d0*real(l) + 1.d0)) * epsil * y(5+l)
         end do
         !set to exit loop once switched to after tc regime
         exit
         end if

    end do

         ! Integrates equations from the end of tight coupling until the start of recombination
         if ( x_tc /= x_t(1)) then
            do i = i_tc_end, n_brec-1
	       
	       call odeint(y, x_brec(i), x_brec(i+1), eps, h1, hmin, derivs_be, bsstep, output_be)

	       ! Stores first non-tc values at time step i_tc_end+1 in global variables
	       a_current = exp(x_brec(i+1))         
               delta(i,k)   = y(1)
               delta_b(i,k) = y(2)
               v(i,k)       = y(3)
               v_b(i,k)     = y(4)
               Phi(i,k)     = y(5)
               do l = 0, lmax_int
                  Theta(i,l,k) = y(6 + l)
               end do
               Psi(i,k)     = -Phi(i,k) - 12.d0*(H_0 /(c * k_current * a_current))**2 * Omega_r * &
			     & Theta(i,2,k)           
            end do
	
        end if
    
         ! Stores derivatives that are required for C_l estimation
	 a_current = exp(x_brec(n_brec))
         call derivs_be(x_brec(n_brec),y,dydx)
         dPhi(0,k)     = dydx(5)
         dv_b(0,k)     = dydx(4)
         dTheta(0,:,k) = dydx(6:6+lmax_int)
         dPsi(0,k)     = -dPhi(0,k) + 12.d0*(H_0 /(c * k_current * a_current))**2 * &
			 &Omega_r * (2.d0 * y(8)-dTheta(0,2,k))

        do i = 1, n_t-1
         ! Defines the values of x, a, and R at the start of recombination
         x1 = x_t(i)
         x2 = x_t(i+1)
         a_current = exp(x2)
         R = (4.d0 * Omega_r) / (3.d0 * Omega_b * a_current)

	 ! Integrates equations from Recombination to today
         call odeint(y, x1, x2, eps, h1, hmin, derivs_be, bsstep, output_be)

         ! Store variables at time step i in global variables. 
         ! Two x-arrays overlap by 1 so n_xtot-1 is total number in array.
         ! Values of global parameters per k-mode.
         delta(n_brec + i-1,k)   = y(1)
         delta_b(n_brec + i-1,k) = y(2)
         v(n_brec + i-1,k)       = y(3)
         v_b(n_brec + i-1,k)     = y(4)
         Phi(n_brec + i-1,k)     = y(5)
         do l = 0, lmax_int
             Theta(n_brec + i-1,l,k) = y(6 + l)
         end do
         Psi(n_brec + i-1,k)     = -Phi(n_brec + i-1,k) - 12.d0*(H_0 / (c * k_current * a_current))**2 &
                    & * Omega_r * Theta(n_brec + i-1,2,k)

         ! Stores derivatives that are required for C_l estimation
         call derivs_be(x2,y,dydx)
         dPhi(i,k)     = dydx(5)
         dv_b(i,k)     = dydx(4)
         dTheta(i,:,k) = dydx(6:6+lmax_int)
         dPsi(i,k)     = -dPhi(i,k) + 12.d0 * (H_0 /(c * k_current * a_current))**2 * Omega_r * &
			 &(2.d0 * y(8) - dTheta(i,2,k))

      end do
	write(*,*) 'k value and x_tc is ', k, x_tc, i_tc_end

!      open(58,file='cmbdata.unf', form='unformatted')
 !      write(58) delta
  !     write(58) delta_b
   !    write(58) v
    !   write(58) v_b
     !  write(58) Theta
      ! write(58) Phi
!       write(58) Psi 
 !      write(58) dPhi
  !     write(58) dPsi
   !    write(58) dTheta
    !   write(58) dv_b
     !  close(58)

    end do
!end if

    ! write parameters to file for given k modes
    open(96, file = 'phi.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Phi(i-1,k)
      end do
    end do

    open(96, file = 'psi.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Psi(i-1,k)
      end do
    end do

    open(96, file = 'delta.dat')
    do k = 0, n_k-1
    do i = 1, n_xtot-1
    write (96,*) x_tot(i), delta(i-1,k)
        end do
    end do
    close(96)
 
    open(96, file = 'delta_b.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), delta_b(i-1,k)
      end do
    end do

    open(96, file = 'v.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), v(i-1,k)
      end do
    end do

    open(96, file = 'v_b.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), v_b(i-1,k)
      end do
    end do

    open(96, file = 'Theta0.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,0,k)
      end do
    end do

    open(96, file = 'Theta1.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,1,k)
      end do
    end do

    open(96, file = 'Theta2.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,2,k)
      end do
    end do

    open(96, file = 'Theta3.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,3,k)
      end do
    end do

    open(96, file = 'Theta4.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,4,k)
      end do
    end do

    open(96, file = 'Theta5.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,5,k)
      end do
    end do

    open(96, file = 'Theta6.dat')
    do k= 0, n_k-1
      do i = 1, n_xtot-1
      write(96,*) x_tot(i), Theta(i-1,6,k)
      end do
    end do

    open(96, file = 'dphi.dat')
    do k= 0, n_k-1
      do i = 1, n_t
      write(96,*) x_t(i), dPhi(i-1,k)
      end do
    end do

    open(96, file = 'dpsi.dat')
    do k= 0, n_k-1
      do i = 1, n_t
      write(96,*) x_t(i), dPsi(i-1,k)
      end do
    end do
 
    open(96, file = 'dv_b.dat')
    do k= 0, n_k-1
      do i = 1, n_t
      write(96,*) x_t(i), dv_b(i-1,k)
      end do
    end do

    open(96, file = 'dTheta0.dat')
    do k= 0, n_k-1
      do i = 1, n_t
      write(96,*) x_t(i), dTheta(i-1,2,k)
      end do
    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  ! Returns the time at which
  ! tight coupling ends. In this project, we define this as either when
  ! x_tc if: dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    real(dp)              :: x_init, x_start_rec, x, x_tc, H_p, dtau, eps, z_start_rec

    x_init 	= log(a_init)
    x    	= x_init
    z_start_rec = 1630.4d0
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination

    ! while loop ends and sets value of x_tc when one of criteria are met
    do while (x < x_start_rec)
    x = x + 1.d-5
      H_p = get_H_p(x)
    dtau = get_dtau(x)
    eps = (c * k) / (H_p * dtau)

      ! sets x_tc if: dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
      if (abs(dtau) < 10.d0 .or. abs(eps) > 0.1d0 .or. x > x_start_rec) then
        x_tc = x
        exit
      end if     
    end do
 
   get_tight_coupling_time = x_tc

  end function get_tight_coupling_time

  ! subroutines for y_tight_coupling differential solver and implementation of equation
  subroutine derivs_tight_coupling(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),                  intent(in)  :: x
    real(dp), dimension(:),    intent(in)  :: y
    real(dp), dimension(:),    intent(out) :: dydx
    real(dp)                		   :: a, n, q, R, k, H_p, dH_p, dtau, ddtau, epsil
    real(dp)            	           :: delta, delta_b, v, v_b, Phi, Theta_0, Theta_1, Theta_2, Psi
    real(dp)                               :: dPhi, dTheta_0, dTheta_1, dv_b, ddelta, ddelta_b, dv

    ! assignation of current values of equations for differential calculations
    delta     = y(1)
    delta_b   = y(2)
    v         = y(3)
    v_b       = y(4)
    Phi       = y(5)
    Theta_0   = y(6)
    Theta_1   = y(7)

    ! input values for current time value (x)
    a     = exp(x)
    k     = k_current
    H_p   = get_H_p(x)
    dH_p  = get_dH_p(x)
    dtau  = get_dtau(x)
    ddtau = get_ddtau(x)
    epsil = (c*k)/(H_p)

    ! Tight Coupling Einstein-Boltzmann equations
    R = (4.d0 * Omega_r) / (3.d0 * Omega_b * a)
    Theta_2 = -(20.d0 * epsil * Theta_1) / (45.d0 * dtau)
    Psi = -Phi  - 12.d0 * (H_0 / (c * k * a))**2 * Omega_r * Theta_2

    ! Use Calculated Psi and Phi values to calculate the derivative of Phi, dPhi
    dPhi = Psi - (epsil * epsil * Phi) / (3.d0) + &
           & ((H_0 * H_0) / (2.d0 * H_p * H_p)) * &
           & (((Omega_m * delta) / a) + ((Omega_b * delta_b) / a) + ((4.d0 * Omega_r * Theta_0) &
	   & / (a * a)))

    !Use the previous to calculate dTheta_0 and then q and dv_b. Then dTheta_1, ddelta, dv, ddelta_b
    dTheta_0 = -(epsil * Theta_1) - dPhi
    q        = ( -( (1.d0 - 2.d0 * R) * dtau + (1.d0 + R) * ddtau ) * (3.d0 * Theta_1 + v_b) - &
               & epsil * Psi +  (1.d0 - dH_p/H_p) * epsil * (-Theta_0 + 2.0 * Theta_2) - &
		 epsil * dTheta_0) / ((1.d0 + R) * dtau + (dH_p/H_p) - 1.d0)
    dv_b     = (1.d0 / (1.d0 + R)) * (-v_b - epsil * Psi + R * (q + epsil * &
             & (-Theta_0 + 2.d0 * Theta_2) - epsil * Psi))
    dTheta_1 = (q - dv_b) / 3.d0
    dv       = -v - epsil * Psi
    ddelta   = epsil * v - 3.d0 * dPhi
    ddelta_b = epsil * v_b - 3.d0 * dPhi

    dydx(1)  = ddelta
    dydx(2)  = ddelta_b
    dydx(3)  = dv
    dydx(4)  = dv_b
    dydx(5)  = dPhi
    dydx(6)  = dTheta_0
    dydx(7)  = dTheta_1

  end subroutine derivs_tight_coupling
 
  subroutine output_tight_coupling(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y

  end subroutine output_tight_coupling

  ! subroutines for Einstein-Boltzmann differential solver and implementation of equation for post-tight coupling regime
  subroutine derivs_be(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    integer(i4b)              	        :: i, l
    real(dp), allocatable, dimension(:) :: Theta, dTheta
    real(dp)                            :: R, a, k, epsil, eta, H_p, dH_p, dtau, ddtau, delta 
    real(dp)                            :: dPhi, dv_b, ddelta, ddelta_b, dv, delta_b, v, v_b, Phi, Psi

    ! create the length of Theta array based on number on lmax_int
    allocate(Theta(0:lmax_int))
    allocate(dTheta(0:lmax_int))

    ! assignation of current values of equations for differential calculations
    delta    = y(1)
    delta_b  = y(2)
    v        = y(3)
    v_b      = y(4)
    Phi      = y(5)
    do i = 0, lmax_int
       Theta(i) = y(6+i)
    end do

    ! input values for current time value (x)
    a        = exp(x)
    k     = k_current
    eta   = get_eta(x)
    H_p   = get_H_p(x)
    dH_p  = get_dH_p(x)
    dtau  = get_dtau(x)
    ddtau = get_ddtau(x)
    R     = (4.d0 * Omega_r) / (3.d0 * Omega_b * a)
    epsil = (c*k)/(H_p)

    ! Calculate Einstein-Boltzmann equations   
    Psi = -Phi  - 12.d0*(H_0/ (c * k * a))**2 * Omega_r * Theta(2)

    ! Use Calculated Psi and Phi values to calculate the derivative of Phi, dPhi
    dPhi = Psi - (epsil * epsil * Phi) / (3.d0) + &
           & ((H_0 * H_0) / (2.d0 * H_p * H_p)) * &
           & (((Omega_m * delta) / a) + ((Omega_b * delta_b) / a) + ((4.d0 * Omega_r * Theta(0))&
	   & / (a * a)))

    ! Use Psi and dPhi to calculate ddelta, ddelta_b, dv, dv_b, and dTheta(0:lmax_int)
    ddelta    = epsil * v - 3.d0 * dPhi
    ddelta_b  = epsil * v_b - 3.d0 * dPhi
    dv        = -v - epsil * Psi
    dv_b      = -v_b - epsil * Psi + dtau * R * (3.d0 * Theta(1) + v_b)
    dTheta(0) = -epsil * Theta(1) - dPhi
    dTheta(1) = (epsil / 3.d0) * Theta(0) - ( 2.d0 * epsil / 3.d0) * Theta(2) + &  
          & (epsil / 3.d0) * Psi + dtau * (Theta(1) + (1.d0/3.d0) * v_b)
    dTheta(2)= (2.d0 * epsil / 5.d0) * Theta(1) - (3.d0/5.d0) * epsil* Theta(3) &
          &  + dtau * (Theta(2) - .1d0 * Theta(2))
    do l = 3, lmax_int-1
     dTheta(l)= ((real(l) * epsil * Theta(l-1)) / (2.d0*real(l) + 1.d0)) - (((real(l) + 1.d0) * &
		& epsil * Theta(l+1)) / ((2.d0*real(l) + 1.d0))) + (dtau * Theta(l))
    end do
    dTheta(lmax_int) = epsil * Theta(lmax_int - 1) - ((c * (lmax_int + 1) * Theta(lmax_int)) / &
             & (H_p * eta)) + (dtau * Theta(lmax_int))

    dydx(1)  = ddelta
    dydx(2)  = ddelta_b
    dydx(3)  = dv
    dydx(4)  = dv_b
    dydx(5)  = dPhi
    dydx(6:6 + lmax_int) = dTheta(:)

  end subroutine derivs_be
 
  subroutine output_be(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y

  end subroutine output_be

end module evolution_mod
