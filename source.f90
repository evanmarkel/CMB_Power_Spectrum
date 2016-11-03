
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

