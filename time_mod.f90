module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver 
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
  real(dp),    allocatable, dimension(:) :: z_t                ! Grid of relevant z-values

  integer(i4b)                           :: n_eta              ! Number of eta grid points
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2	       ! Eta and eta'' at each grid point

contains
  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2 
    real(dp)	 :: a_init, eps, h1, hmin, x_splint, y_splint, omega_total
    real(dp), dimension(1) :: y, x

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Total number of grid points
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today 
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    eps   	= 1.d-10	            ! initial values for odeint solver
    h1    	= 1.d-5			    
    hmin  	= 0.d0

    ! Allocate x array
    allocate(x_t(n_t))

    ! x values in recombination era. 200 grid points in the range z = [1630.4, 614].
    dx = (x_end_rec - x_start_rec) / (n1 - 1)
    do i = 1 , n1
       x_t(i) = x_start_rec + (i - 1) * dx
    end do

    ! x va lues in post-recombination era. 300 grid points in the range z = (614, 0]. 
    dx = (x_0 - x_end_rec) / (n2)
    do i = 1 , n2
	x_t(n1 + i) = x_end_rec + i * dx
    end do

    ! universal expansion scale factor a = exp(x) grid. 
    allocate(a_t(n_t))
    allocate(z_t(n_t))

       a_t = exp(x_t)
       z_t = 1.d0 / a_t - 1.d0

    !allocate the arrays for computing the conformal time eta(x_eta). x_eta contains n_eta=1000 points. 
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))
    
    !set initial values for eta, corresponding to a = 10.d-10
    x(1) = x_eta1
    y(1) = c / ((H_0 * sqrt(Omega_r))) * a_init
    eta(1) = y(1)

    !determine step length dx and use to integrate eta, solving at each x_eta grid point
    dx = (x_eta2 - x_eta1) / (n_eta - 1)
    do i = 1 , n_eta
	x_eta(i) = x_eta1 + (i - 1) * dx
        !Integrate d(eta)/dx using the burlisch-storr method.  
  	call odeint(y, x_eta(i), x_eta(i) + dx, eps, h1, hmin, derivs_eta, bsstep, output_eta)  
	eta(i+1) = y(1)
    end do

    ! spline the integrated function used for spline to calculate eta''. 
    ! yp1, yp2 are large so eta'' vanishes at end points. eta'' to be used by splint function to 
    ! accurately calculate values of eta in between grid points
    call spline(x_eta, eta, 1.d30, 1.d30, eta2)

    ! create dat files and print array data to graph

    ! eta values for x_eta steps
    open(96, file = 'eta.dat')
    do i = 1 , n_eta
    write (96,*) x_eta(i) , eta(i) / Mpc
    end do 
    close(96)

    ! hubble constant for x_t in (km/s) / Mpc
    open(96, file = 'H(x).dat')
    do i = 1 , n_eta
    write(96,*) x_eta(i), get_H(x_eta(i)) * Mpc / 1.d3
    end do 
    close(96)	

    ! hubble constant function of z and a in (km/s) / Mpc 
    open(96, file = 'H(z).dat')
    do i = 1 , n_t
    write(96,*) z_t(i), a_t(i), get_H(x_t(i)) * Mpc / 1.d3
    end do 
    close(96)

    ! hubble prime constant for x_t in (km/s) / Mpc
    open(96, file = 'H_p(x).dat')
    do i = 1 , n_eta
    write(96,*) x_eta(i), get_H_p(x_eta(i)) * Mpc / 1.d3
    end do
    close(96)

    ! derivative of hubble prime constant for x_t in (km/s) / Mpc
    open(96, file = 'dH_p(x).dat')
    do i = 1 , n_eta
    write(96,*) x_eta(i), get_dH_p(x_eta(i)) * Mpc / 1.d3
    end do
    close(96)

    ! interpolate eta over the grid range by a factor of ten
    open(96, file = 'interpolate_eta.dat')
    do i = 1 , 10*n_eta
    x_splint = x_eta1 + 1.d-1 * (i - 1) * dx
    y_splint = get_eta(x_splint)
    write(96,*) x_splint, y_splint / Mpc
    end do
    close(96)

    ! plot the relative energy density parameters for baryons, dark matter, radiation, and lambda as a function of x
    ! omega_total normalizes the densities to 1 assuming a flat universe
    open(96, file = 'energy.dat')
    do i = 1 , n_eta
    omega_total =   (Omega_b * (exp(x_eta(i))**(-3))) / rho_c + & 
                  & (Omega_m * (exp(x_eta(i))**(-3))) / rho_c + &
                  & (Omega_r * (exp(x_eta(i))**(-4))) / rho_c + &
		  & (Omega_lambda) / rho_c
    write(96,*) x_eta(i), &
                & ((Omega_b * (exp(x_eta(i))**(-3))) / rho_c) / omega_total, &
	        & ((Omega_m * (exp(x_eta(i))**(-3))) / rho_c) / omega_total, &
	        & ((Omega_r * (exp(x_eta(i))**(-4))) / rho_c) / omega_total, &
		& (Omega_lambda / rho_c) / omega_total
    end do
    close(96)

  end subroutine initialize_time_mod

  subroutine derivs_eta(x, y, detada)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y
        real(dp), dimension(:), intent(out) :: detada

        detada(1) = c / get_H_p(x)

  end subroutine derivs_eta
  
  subroutine output_eta(x, y)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: y

  end subroutine output_eta

  ! Here is a function that computes H at given x, using the Friedmann equation using params module
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H, a
 
    a = exp(x)
    get_H = H_0 * sqrt( (Omega_b + Omega_m) * (a ** (-3)) + (Omega_r) * (a ** (-4)) + Omega_lambda)

  end function get_H

  ! This function computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p, a

    a = exp(x)

    get_H_p = H_0 * sqrt( (Omega_b + Omega_m) / a + (Omega_r) * (a**(-2)) + & 
       & Omega_lambda * (a**(2)))

  end function get_H_p

  ! Task: Write a function that computes dH_p/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    real(dp)             :: a

    a = exp(x)
    get_dH_p = 0.5d0 * a * H_0 * ( - (Omega_b + Omega_m) * a ** (-2) - 2.d0 * (Omega_r) * a ** (-3) + &
             & 2.d0 * (Omega_lambda) * a) / sqrt( ((Omega_b + Omega_m) / a) + (Omega_r) * a ** (-2) + & 
             & (Omega_lambda) * a ** (2))

  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    get_eta = splint(x_eta, eta, eta2, x_in)

  end function get_eta

end module time_mod
