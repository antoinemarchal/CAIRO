!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  use mod_inout

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa(data, std_data, grid_params, fileout, timeout, n_gauss, n_gauss_add, lambda_amp, lambda_mu, &
       lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, amp_fact_init, sig_init, lb_sig_init, &
       ub_sig_init, lb_sig, ub_sig, maxiter_init, maxiter, m, noise, regul, descent, lstd, ustd, init_option, &
       iprint, iprint_init, save_grid, lym, params_init, init_spec, norm_var, lambda_sig_corr_narrow, &
       lambda_sig_corr_broad, lambda_mu_corr_narrow, lambda_mu_corr_broad, ib_sig, lambda_r, lb_amp, ub_amp, params_lsf)
    
    implicit none
    
    logical, intent(in) :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
    logical, intent(in) :: regul           !! if true --> activate regulation
    logical, intent(in) :: descent         !! if true --> activate hierarchical descent to initiate the optimization
    logical, intent(in) :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
    logical, intent(in) :: lym             !! if true --> activate 2-Gaussian decomposition for Lyman alpha nebula emissio
    logical, intent(in) :: init_spec       !! if true --> use params mean spectrum with input
    logical, intent(in) :: norm_var        !! if true --> normalize the var sig energy term

    integer, intent(in) :: n_gauss_add     !! number of gaussian to add at each step
    integer, intent(in) :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: iprint          !! print option 
    integer, intent(in) :: iprint_init     !! print option init
    integer, intent(in) :: maxiter         !! max iteration for L-BFGS-B alogorithm
    integer, intent(in) :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

    real(xp), intent(in) :: lambda_amp     !! lambda for amplitude parameter
    real(xp), intent(in) :: lambda_mu      !! lamnda for mean position parameter
    real(xp), intent(in) :: lambda_sig     !! lambda for dispersion parameter

    real(xp), intent(in) :: lambda_var_amp !! lambda for amp dispersion parameter
    real(xp), intent(in) :: lambda_var_mu  !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_var_sig !! lambda for variance dispersion parameter

    real(xp), intent(in) :: lambda_lym_sig !! lambda for difference dispersion parameter (2-gaussaian)

    real(xp), intent(in) :: lambda_sig_corr_narrow !! lambda for coorelation between narrow sigma field
    real(xp), intent(in) :: lambda_sig_corr_broad  !! lambda for coorelation between broad sigma field

    real(xp), intent(in) :: lambda_mu_corr_narrow !! lambda for coorelation between narrow position field
    real(xp), intent(in) :: lambda_mu_corr_broad  !! lambda for coorelation between broad position field

    real(xp), intent(in) :: lambda_r !! lambda for ratio 1/3 NII or OIII lines

    real(xp), intent(in) :: amp_fact_init  !! times max amplitude of additional Gaussian
    real(xp), intent(in) :: sig_init       !! dispersion of additional Gaussian
    real(xp), intent(in) :: ub_sig_init    !! upper bound sigma init
    real(xp), intent(in) :: lb_sig_init    !! lower bound sigma init
    real(xp), intent(in) :: lb_sig         !! lower bound sigma
    real(xp), intent(in) :: ub_sig         !! upper bound sigma
    real(xp), intent(in) :: ib_sig         !! intermediate bound sigma
    real(xp), intent(in) :: lb_amp         !! lower bound amplitude
    real(xp), intent(in) :: ub_amp         !! upper bound amplitude

    character(len=8), intent(in)   :: init_option !!Init ROHSA with the mean or the std spectrum    
    character(len=512), intent(in) :: fileout   !! name of the output result
    character(len=512), intent(in) :: timeout   !! name of the output result

    integer :: n_gauss      !! number of gaussian to fit
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    real(xp), intent(in), dimension(:,:,:), allocatable :: data        !! initial fits data
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_data    !! standard deviation map fo the cube is given by the user 
    real(xp), intent(in), dimension(:), allocatable     :: params_init !!
    real(xp), intent(in), dimension(:), allocatable     :: params_lsf  !!

    real(xp), intent(inout), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cub

    real(xp), dimension(:,:,:), allocatable :: cube            !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_mean       !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: fit_params      !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:,:), allocatable :: std_cube          !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:,:,:), allocatable :: std_map           !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:,:,:), allocatable :: std_map_2           !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:), allocatable :: b_params            !! unknow average sigma
    real(xp), dimension(:), allocatable :: std_spect           !! std spectrum of the observation
    real(xp), dimension(:), allocatable :: max_spect           !! max spectrum of the observation
    real(xp), dimension(:), allocatable :: max_spect_norm      !! max spectrum of the observation normalized by the max of the mean spectrum
    real(xp), dimension(:), allocatable :: mean_spect          !! mean spectrum of the observation
    real(xp), dimension(:), allocatable :: guess_spect         !! params obtain fi the optimization of the std spectrum of the observation

    real(xp) :: c_lym=1._xp !! minimized the variance of the ratio between dispersion 1 and dispersion of a 2-Gaussian model for Lym alpha nebula

    integer, dimension(3) :: dim_data !! dimension of original data
    integer, dimension(3) :: dim_cube !! dimension of reshape cube
    
    real(xp), dimension(:,:), allocatable :: kernel !! convolution kernel 

    real(xp) :: lctime, uctime
    
    integer :: ios=0 !! ios integer
    integer :: i     !! loop index
    integer :: j     !! loop index
    integer :: k     !! loop index
    integer :: l     !! loop index

    integer :: delta_1 ! delta position between center line and line on its left
    integer :: delta_2 ! delta position between center line and line on its right
    integer :: delta_3 ! delta position between center line and line on its right
    integer :: delta_4 ! delta position between center line and line on its right

    integer, dimension(6) :: position = 0._xp
    integer, dimension(6) :: delta_balmer = 0._xp
    integer, dimension(6) :: delta_forbidden = 0._xp
        
    print*, "fileout = '",trim(fileout),"'"
    print*, "timeout = '",trim(timeout),"'"
    
    print*, " "
    print*, "______Parameters_____"
    print*, "n_gauss = ", n_gauss

    print*, "lambda_amp = ", lambda_amp
    print*, "lambda_mu = ", lambda_mu
    print*, "lambda_sig = ", lambda_sig

    print*, "lambda_var_sig = ", lambda_var_sig

    print*, "lambda_sig_corr_narrow = ", lambda_sig_corr_narrow
    print*, "lambda_sig_corr_broad = ", lambda_sig_corr_broad

    print*, "lambda_mu_corr_narrow = ", lambda_mu_corr_narrow
    print*, "lambda_mu_corr_broad = ", lambda_mu_corr_broad

    print*, "lambda_r = ", lambda_r

    ! print*, "amp_fact_init = ", amp_fact_init
    ! print*, "sig_init = ", sig_init
    ! print*, "lb_sig_init = ", lb_sig_init
    ! print*, "ub_sig_init = ", ub_sig_init
    print*, "lb_sig = ", lb_sig
    print*, "ib_sig = ", ib_sig
    print*, "ub_sig = ", ub_sig
    print*, "lb_amp = ", lb_amp
    print*, "ub_amp = ", ub_amp
    ! print*, "init_option = ", init_option
    ! print*, "maxiter_init = ", maxiter_init
    print*, "maxiter = ", maxiter
    ! print*, "lstd = ", lstd
    ! print*, "ustd = ", ustd
    ! print*, "noise = ", noise
    ! print*, "regul = ", regul
    ! print*, "descent = ", descent
    ! print*, "save_grid = ", save_grid
    ! print*, "init_spec = ", init_spec
    ! print*, "norm_var = ", norm_var

    print*, " "

    ! Check n_gauss = 2 for Lym akpha mode
    if (lym .eqv. .true.) then
       ! if (n_gauss .eq. 2) then
          print*, "Lym alpha mode activated"
          stop
       ! else 
       !    print*, "Lym alpha mode is based on a 2-Gaussian model / please select n_gauss = 2"
       !    stop
       ! end if
    end if

    print*, " "
    
    allocate(kernel(3, 3))
    
    kernel(1,1) = 0._xp
    kernel(1,2) = -0.25_xp
    kernel(1,3) = 0._xp
    kernel(2,1) = -0.25_xp
    kernel(2,2) = 1._xp
    kernel(2,3) = -0.25_xp
    kernel(3,1) = 0._xp
    kernel(3,2) = -0.25_xp
    kernel(3,3) = 0._xp
        
    dim_data = shape(data)
    
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    write(*,*) ""
    write(*,*) "number of los = ", dim_data(2)*dim_data(3)
    
    nside = dim2nside(dim_data)
    
    write(*,*) "nside = ", nside
    
    call dim_data2dim_cube(nside, dim_data, dim_cube)
    
    !Allocate memory for cube
    allocate(cube(dim_cube(1), dim_cube(2), dim_cube(3)))
    allocate(std_cube(dim_cube(1), dim_cube(2), dim_cube(3)))
    
    !Reshape the data (new cube of size nside)
    print*, " "
    write(*,*) "Reshape cube, new dimensions :"
    write(*,*) "dim_v, dim_y, dim_x = ", dim_cube
    print*, " "
    
    print*, "Compute mean and std spectrum"
    allocate(std_spect(dim_data(1)))
    allocate(max_spect(dim_data(1)), max_spect_norm(dim_data(1)))
    allocate(mean_spect(dim_data(1)))
    allocate(b_params(n_gauss))

    call std_spectrum(data, std_spect, dim_data(1), dim_data(2), dim_data(3))
    call mean_spectrum(data, mean_spect, dim_data(1), dim_data(2), dim_data(3))
    call max_spectrum(data, max_spect, dim_data(1), dim_data(2), dim_data(3))
    call max_spectrum(data, max_spect_norm, dim_data(1), dim_data(2), dim_data(3), maxval(mean_spect))
    
    call reshape_up(data, cube, dim_data, dim_cube)
    
    !Allocate memory for parameters grids
    if (descent .eqv. .true.) then
       allocate(fit_params(3*n_gauss, 1, 1))
       !Init sigma = 1 to avoid Nan
       do i=1,n_gauss
          fit_params(1+(3*(i-1)),1,1) = 0._xp
          fit_params(2+(3*(i-1)),1,1) = 1._xp
          fit_params(3+(3*(i-1)),1,1) = 1._xp
       end do
    end if
    
    print*, "                    Start iteration"
    print*, " "
    
    if (descent .eqv. .true.) then
       print*, "Start hierarchical descent"

       if (save_grid .eqv. .true.) then
          !Open file time step
          open(unit=11, file=timeout, status='replace', access='append', iostat=ios)
          write(11,fmt=*) "# size grid, Time (s)"
          close(11)
          call cpu_time(lctime)
       end if

       !Start iteration
       do n=0,nside-1
          power = 2**n
          
          allocate(cube_mean(dim_cube(1), power, power))
          
          call mean_array(power, cube, cube_mean)
          
          if (n == 0) then
             if (init_spec .eqv. .true.) then
                print*, "Use user init params"
                fit_params(:,1,1) = params_init
                print*, "Calculate Delta between lines"
                do i=1, 6
                   position(i) = params_init(2+(3*(i-1)))
                end do
                delta_1 = position(5) - position(1)
                delta_2 = position(6) - position(2)
                delta_3 = position(6) - position(3)
                delta_4 = position(6) - position(4)
                print*, "Delta Halpha and Hbeta = ", delta_1
                print*, "Delta Forbidden lines = ", delta_2
                print*, "Delta Forbidden lines = ", delta_3
                print*, "Delta Forbidden lines = ", delta_4
             else
                if (init_option .eq. "mean") then
                   print*, "Init mean spectrum"        
                   print*, "init_spec == .false. not available"
                   stop
                   ! call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), cube_mean(:,1,1), amp_fact_init, sig_init, &
                   !      lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
                else 
                   print*, "init_option keyword should be 'mean'"
                   stop
                end if
             end if

             !Init b_params
             do i=1, n_gauss       
                b_params(i) = fit_params(3+(3*(i-1)),1,1)
             end do
          end if
                    
          if (regul .eqv. .false.) then
             print*, "regul == .false. not available"
             stop
             ! call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), lb_sig, ub_sig, maxiter, m, iprint)
          end if
          
          if (regul .eqv. .true.) then
             if (n == 0) then                
                print*,  "Update level", n
                ! call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), lb_sig, ub_sig, maxiter, m, iprint)
             end if

             if (n > 0 .and. n < nside) then
                allocate(std_map(dim_cube(1), power, power))
                
                if (noise .eqv. .true.) then
                   call reshape_up(std_data, std_cube, dim_data, dim_cube)
                   call sum_array_square(power, std_cube, std_map)
                   std_map = sqrt(std_map) / real(2**(2*(nside-n)),xp)
                else
                   print*, "no noise = .false. in this branch"
                   stop
                end if

                ! Update parameters 
                print*,  "Update level", n, ">", power
                call update(cube_mean, fit_params, b_params, n_gauss, dim_cube(1), power, power, lambda_amp, lambda_mu, &
                     lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, lb_sig, ub_sig, maxiter, &
                     m, kernel, iprint, std_map, lym, c_lym, norm_var, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
                     lambda_mu_corr_narrow, lambda_mu_corr_broad, ib_sig, lambda_r, lb_amp, ub_amp, params_lsf, &
                     delta_1, delta_2, delta_3, delta_4)

                deallocate(std_map)
             end if
          end if
          
          deallocate(cube_mean)

          ! Save grid in file
          if (save_grid .eqv. .true.) then
             print*, "Save grid parameters"
             call save_process(n, n_gauss, fit_params, power, fileout)
             !Save timestep
             if (n .ne. 0) then
                open(unit=11, file=timeout, status='unknown', access='append', iostat=ios)
                if (ios /= 0) stop "opening file error"
                call cpu_time(uctime)
                print*, dim_cube(1)
                ! write(11,fmt='(I6, f6.3)') power, uctime-lctime
                write(11,fmt=*) power, uctime-lctime
                close(11)
             end if
          end if
          
          ! Propagate solution on new grid (higher resolution)
          call go_up_level(fit_params)
          write(*,*) " "
          write(*,*) "Interpolate parameters level ", n!, ">", power

       enddo
       
       print*, " "
       write(*,*) "Reshape cube, restore initial dimensions :"
       write(*,*) "dim_v, dim_y, dim_x = ", dim_data
              
       call reshape_down(fit_params, grid_params,  (/ 3*n_gauss, dim_cube(2), dim_cube(3)/), &
            (/ 3*n_gauss, dim_data(2), dim_data(3)/))       

       else
          allocate(guess_spect(3*(n_gauss+n_gauss_add)))
          if (init_option .eq. "mean") then
             print*, "Use of the mean spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), mean_spect, amp_fact_init, sig_init, &
                  lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
          else if (init_option .eq. "std") then
             print*, "Use of the std spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), std_spect, amp_fact_init, sig_init, &
                  lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
          else if (init_option .eq. "max") then
             print*, "Use of the max spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), max_spect, amp_fact_init, sig_init, &
                  lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
          else if (init_option .eq. "maxnorm") then
             print*, "Use of the std spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), max_spect_norm, amp_fact_init, sig_init, &
                  lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
          else
             print*, "init_option keyword should be 'mean' or 'std' or 'max'"
             stop
          end if
          call init_grid_params(grid_params, guess_spect, dim_data(2), dim_data(3))

          deallocate(guess_spect)      
    end if
    
    !Update last level
    print*, " "
    print*, "Start updating last level."
    print*, " "
    
    allocate(std_map(dim_data(1), dim_data(2), dim_data(3)))
    
    if (noise .eqv. .true.) then
       std_map = std_data
    else   
       print*, "no noise = .false. in this branch"
       stop
    end if
    
    if (regul .eqv. .true.) then
       call update(data, grid_params, b_params, n_gauss, dim_data(1), dim_data(2), dim_data(3), lambda_amp, lambda_mu, &
            lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, lb_sig, ub_sig, maxiter, m, &
            kernel, iprint, std_map, lym, c_lym, norm_var, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
            lambda_mu_corr_narrow, lambda_mu_corr_broad, ib_sig, lambda_r, lb_amp, ub_amp, params_lsf, &
            delta_1, delta_2, delta_3, delta_4)
    end if
        
  end subroutine main_rohsa
  
end module mod_rohsa
