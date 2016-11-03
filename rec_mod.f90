module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                  	     ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             	     ! Grid
  real(dp), allocatable, dimension(:), private :: tau, log_tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, log_n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22                 ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y_saha, eps, h1, hmin, qa, qb, qc, qd, T_b, n_b, xmin, xmax, dx, xstart, xstop, a
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H
    real(dp), dimension(1) :: y, x

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstop

    eps   	= 1.d-5       ! step parameters for odeint solver
    h1    	= 1.d-10		    
    hmin  	= 1.d-30

    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(log_tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(log_n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! x_rec grid for n points ranging from xtart to xstop 
    dx = (xstop - xstart) / (n - 1)
    do i = 1 , n
	x_rec(i) = xstart + (i - 1) * dx
    end do 


    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n
       if (use_saha) then
          ! Use the Saha equation. Define a and terms dependent on a. 
	  a = exp(x_rec(i))
	  n_b = (Omega_b * rho_c) / (m_H * a**3)
 	  T_b = T_0 / a

	  ! quadratic formula terms. takes the positive root for electron density. 
	  qa =    1.d0
	  qb =   (1.d0 / n_b) * ((m_e * T_b * k_b) / (2 * pi * hbar * hbar)) ** &
                & 1.5d0 * (exp(-epsilon_0 / (T_b * k_b)))
	  qc = - qb
	  qd = sqrt(qb * qb - 4.d0 * qc)
		!y_saha = (-qb + sqrt(qb * qb - 4.d0 * qc)) / (2.d0*qa)
		! Taylor approximation is used to estimate solution for Saha equation 
		! as was found to be more stable
		y_saha = 1.d0 - 1.d0 / qb + 2.d0 / (qb**2)
	
	  X_e(i) = y_saha
          n_e(i) = y_saha * n_b
 
          if (X_e(i) < saha_limit) then
		use_saha = .false.
		  !Initial conditions for the Peebles equation. 
		  !The last Saha values are used as the first integration point for odeint
		  x(1) = x_rec(i)
	          y(1) = X_e(i)
	 	  X_e(i) = y(1)
		  n_e(i) = X_e(i) * n_b
	  endif 
       else
          ! The Peebles equation is used to calculate electron density 
	  ! First define a, n_b, and T_b 
	  a = exp(x_rec(i))
	  n_b = (Omega_b * rho_c) / (m_H * a**3)
 	  T_b = T_0 / a

          ! Integrate dX_e/dx using the Burlisch-Storr method 
	  ! set value of fractional electron density X_e 
    	  call odeint(y, x_rec(i-1), x_rec(i), eps, h1, hmin, derivs_peebles, bsstep, output_peebles)  
	  X_e(i) = y(1)

	  !set the value of the electron density n_e for calculation of tau
	  n_e(i) = X_e(i) * n_b
       end if
    end do
   
    ! Computes splined (log of) electron density function
    log_n_e = log(n_e)
    call spline(x_rec, log_n_e, 1.d30, 1.d30, n_e2)

    ! Computes optical depth tau at all grid times. 
    ! tau = 0 at n so iterative limit is (1:n-1) to avoid log(0.d0)
    ! Set initial conditions
    y(1) = 0.d0 
    tau(n) = y(1)

    ! tau = 0 at n so iterative limit is (1:n-1) to avoid log(0.d0)
    do i = 1 , n-1
        !Integrate d(tau)/dx using the Burlisch-Storr method.  
  	call odeint(y, x_rec(n+1-i), x_rec(n-i), eps, h1, hmin, derivs_tau, bsstep, output_tau)  
	tau(n-i) = y(1)
    end do

    ! Here we compute the splined (log of) optical depth. 
    ! We only evaluate the log where the value is non-zero.  
    ! Therefore we need to manually set tau at the present day to zero. 
    log_tau(n) = 0.d0
    log_tau(1:n-1) = log(tau(1:n-1))
    call spline(x_rec(1:n-1), log_tau(1:n-1), 1.d30, 1.d30, tau2(1:n-1))

    ! Computes splined second derivative of (log of) optical depth
    call spline(x_rec(1:n-1), tau2(1:n-1), 1.d30, 1.d30, tau22(1:n-1))

    ! Computes splined visibility function
    do i = 1, n
    g(i) = -get_dtau(x_rec(i)) * exp(-get_tau(x_rec(i)))
    end do
    call spline(x_rec, g, 1.d30, 1.d30, g2)

    ! Computes splined second derivative of visibility function
    call spline(x_rec, g2, 1.d30, 1.d30, g22)

    ! writes z and X(e) to file
    open(96, file = 'elecdensity.dat')
    do i = 1 , n
    write (96,*) 1.d0 / exp(x_rec(i)) - 1.d0 , X_e(i)
    end do 
    close(96)

    ! writes x and tau(x) to file
    open(96, file = 'tau.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_tau(x_rec(i))
    end do 
    close(96)

    ! writes x, tau(x), dtau(x), ddtau(x) to file
    open(96, file = 'taudtau.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_tau(x_rec(i)), abs(get_dtau(x_rec(i))), get_ddtau(x_rec(i))
    end do 
    close(96)

    ! writes x and dtau(x) to file
    open(96, file = 'dtau.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_dtau(x_rec(i))
    end do 
    close(96)

    ! writes x and ddtau(x) to file
    open(96, file = 'ddtau.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_ddtau(x_rec(i))
    end do 
    close(96)

    ! writes x and g~(x) to file
    open(96, file = 'g.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_g(x_rec(i))
    end do 
    close(96)

    ! writes x and dg~(x) to file
    open(96, file = 'dg.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_dg(x_rec(i))
    end do 
    close(96)

    ! writes x and ddg~(x) to file
    open(96, file = 'ddg.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_ddg(x_rec(i))
    end do 
    close(96)

    ! writes x, g~(x), dg~(x), and ddg~(x) to file
    open(96, file = 'allg.dat')
    do i = 1 , n
    write (96,*) x_rec(i), get_g(x_rec(i)), get_dg(x_rec(i)) / 10.d0, get_ddg(x_rec(i)) / 300.d0
    end do 
    close(96)

  end subroutine initialize_rec_mod

  !subroutines for peebles integration solver and implementation of Peebles equation 
  subroutine derivs_peebles(x, y, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y
        real(dp), dimension(:), intent(out) :: dydx

        real(dp)     :: phi2, alpha2, beta, beta2, n1s, la, l2s1s, cr, a, n_b, T_b

	a = exp(x)
	n_b = (Omega_b * rho_c) / (m_H * a**3)
 	T_b = T_0 / a

        phi2   = 0.448d0 * log(epsilon_0 / (T_b * k_b))
	alpha2 = ((64.d0 * pi * hbar * hbar) / (sqrt(27.d0 * pi) * c) * &
	        & ((alpha) ** 2) / (m_e) ** 2) * sqrt(epsilon_0 / (T_b * k_b)) * phi2
	beta   = alpha2 * (((m_e * k_b * T_b) / (2.d0 * pi * hbar *hbar)) ** 1.5d0)
	beta2  = beta * exp( -(epsilon_0) / (4.d0 * k_b * T_b) )
	n1s    = (1.d0 - y(1)) * n_b
        la     = (get_H(x) * ((3.d0 * epsilon_0) ** 3)) / (((8.d0 * pi) ** 2) &
	         & * n1s * (hbar * c) **3)
 	l2s1s  = 8.227d0
 	cr 	 = (l2s1s + la) / (l2s1s + la + beta2) 

        dydx(1) = (cr / get_H(x)) * (beta * exp( -(epsilon_0) / (k_b * T_b) ) * &
		 &(1.d0 - y(1)) - n_b * alpha2 * y(1) * y(1))

  end subroutine derivs_peebles
  
  subroutine output_peebles(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y

  end subroutine output_peebles
  
  ! subroutines for tau differential solver and implementation of equation
  subroutine derivs_tau(x, y, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y
        real(dp), dimension(:), intent(out) :: dydx
	real(dp) 			    :: a, n

	a = exp(x)
        dydx(1) = - (get_n_e(x) * sigma_T * a * c) / get_H_p(x)

  end subroutine derivs_tau
  
  subroutine output_tau(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y

  end subroutine output_tau

  ! Completes routine for computing n_e at arbitrary x, using precomputed information
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e

    get_n_e = exp(splint(x_rec, log_n_e, n_e2, x))

  end function get_n_e

  ! Completes routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

  if(x > x_rec(n-1)) then 
	get_tau = 0.d0
    else
    get_tau = exp(splint(x_rec(1:n-1), log_tau(1:n-1), tau2(1:n-1), x))
    endif

  end function get_tau

  ! Completes routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

    if(x > x_rec(n-1)) then 
	get_dtau = 0.d0
    else
    get_dtau = get_tau(x) * splint_deriv(x_rec(1:n-1), log_tau(1:n-1), tau2(1:n-1), x)
    endif

  end function get_dtau

  ! Completes routine for computing the second derivative of tau at arbitrary x, using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau, t, dt

    t  = get_tau(x)
    dt = get_dtau(x)
    if(x > x_rec(n-1)) then 
	get_ddtau = 0.d0
    else 
    get_ddtau = splint(x_rec(1:n-1), tau2(1:n-1), tau22(1:n-1), x) * t + ((dt * dt) / t)
    endif

  end function get_ddtau

  ! Completes routine for computing the visibility function, g, at arbitrary x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(x_rec, g, g2, x)

  end function get_g

  ! Completes routine for computing the derivative of the visibility function, g, at arbitrary x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(x_rec, g, g2, x)
  end function get_dg

  ! Completes routine for computing the second derivative of the visibility function, g, at arbitrary x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg

end module rec_mod
