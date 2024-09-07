!! This module contains optimization subroutine and parametric model
module mod_optimize
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array
  use mod_model 

  implicit none
  
  private

  public :: myfunc_spec, myresidual, f_g_cube_fast_norm, f_g_cube_fast_norm_single, &
       gaussian, dG_da, dG_dmu, dG_dsig, gaussian_lsf, dGlsf_da, dGlsf_dmu, dGlsf_dsig
  
contains

  !!MODEL WITH LSF
  pure function gaussian_lsf(x, a, m, s, lsf)
    !! Gaussian function with LSF
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp), intent(in) :: lsf
    real(xp) :: gaussian_lsf

    gaussian_lsf = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * (lsf**2 + s**2)) );
  end function gaussian_lsf

  pure function dGlsf_da(x, m, s, lsf)
    !! Gradient Gaussian function / a with LSF
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: m, s
    real(xp), intent(in) :: lsf
    real(xp) :: dGlsf_da

    dGlsf_da = exp(-( (real(x,xp) - m)**2 ) / (2._xp * (lsf**2 + s**2)) );
  end function dGlsf_da
  
  pure function dGlsf_dmu(x, a, m, s, lsf)
    !! Gradient Gaussian function / mu
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp), intent(in) :: lsf
    real(xp) :: dGlsf_dmu
    
    dGlsf_dmu = a * (real(x,xp) - m) / (lsf**2._xp + s**2._xp) &
         * exp(-( (real(x,xp) - m)**2 ) / (2._xp * (lsf**2 + s**2)) );
  end function dGlsf_dmu

  pure function dGlsf_dsig(x, a, m, s, lsf)
    !! Gradient Gaussian function / sig
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp), intent(in) :: lsf
    real(xp) :: dGlsf_dsig

    dGlsf_dsig = a * s * (real(x,xp) - m)**2._xp / (lsf**2._xp + s**2._xp)**2._xp &
         * exp(-( (real(x,xp) - m)**2 ) / (2._xp * (lsf**2 + s**2)) );
  end function dGlsf_dsig

  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function gaussian

  pure function dG_da(x, m, s)
    !! Gradient Gaussian function / a
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: m, s
    real(xp) :: dG_da

    dG_da = exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function dG_da

  pure function dG_dmu(x, a, m, s)
    !! Gradient Gaussian function / mu
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: dG_dmu
    
    dG_dmu = a * (real(x,xp) - m) / s**2._xp * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function dG_dmu

  pure function dG_dsig(x, a, m, s)
    !! Gradient Gaussian function / sig
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: dG_dsig

    dG_dsig = a * (real(x,xp) - m)**2._xp / s**3._xp * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function dG_dsig
    
  ! Compute the residual between model and data
  subroutine myresidual(params, line, residual, n_gauss, dim_v, lsf)
    implicit none

    integer, intent(in) :: dim_v, n_gauss
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(in), dimension(3*n_gauss) :: params
    real(xp), intent(inout), dimension(:), allocatable :: residual
    real(xp), intent(in), dimension(n_gauss) :: lsf

    integer :: i, k
    real(xp) :: g    
    real(xp), dimension(dim_v) :: model

    integer :: start, finish
    real(xp) :: delta

    g = 0._xp
    model = 0._xp
    
    do i=1, n_gauss
       ! delta = 5._xp * params(3+(3*(i-1)))
       delta = 5._xp * sqrt(params(3+(3*(i-1)))**2._xp + lsf(i)**2._xp)
       start = int(params(2+(3*(i-1))) - delta)
       finish = int(params(2+(3*(i-1))) + delta)
       if (start .le. 1) then 
          start = 1
       end if
       if (finish > dim_v) then 
          finish = dim_v
       end if
       ! print*, start, end
       ! stop
       ! do k=1, dim_v
       do k=start, finish
          ! g = gaussian(k, params(1+(3*(i-1))), params(2+(3*(i-1))), params(3+(3*(i-1))))
          g = gaussian_lsf(k, params(1+(3*(i-1))), params(2+(3*(i-1))), params(3+(3*(i-1))), lsf(i))
          model(k) = model(k) + g
       enddo
    enddo

    residual = model - line
  end subroutine myresidual

  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_fast_norm(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
       lambda_mu_corr_narrow, lambda_mu_corr_broad, lambda_r, params_lsf, delta_1, delta_2, delta_3, delta_4)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    integer, intent(in) :: delta_1, delta_2, delta_3, delta_4
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig, lambda_r
    real(xp), intent(in) :: lambda_sig_corr_narrow, lambda_sig_corr_broad
    real(xp), intent(in) :: lambda_mu_corr_narrow, lambda_mu_corr_broad
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:), allocatable :: params_lsf
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    integer :: n_beta
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:), allocatable :: residual_noise
    real(xp), dimension(:,:,:), allocatable :: params
    real(xp), dimension(:), allocatable :: b_params
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
    real(xp), dimension(:,:,:), allocatable :: deriv
    real(xp), dimension(:), allocatable :: model

    integer :: start, finish
    real(xp) :: delta

    allocate(deriv(3*n_gauss, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(b_params(n_gauss))
    allocate(params(3*n_gauss, dim_y, dim_x))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    allocate(model(dim_v))

    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    residual = 0._xp    
    params = 0._xp
    model = 0._xp
    
    n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
    do i=1,n_gauss
       b_params(i) = beta((n_beta-n_gauss)+i)
    end do

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          allocate(residual_noise(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v, params_lsf)
          residual(:,i,j) = residual_1D
          if (ABS(std_map(1,i,j)) > 0._xp) then
             residual_noise = residual_1D / std_map(:,i,j)
             f = f + (myfunc_spec(residual_noise))
          end if
          deallocate(residual_1D)
          deallocate(residual_noise)
       end do
    end do

    ! Compute the objective function and the gradient
    do i=1, n_gauss
       !
       conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
       conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
       image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       
       image_amp = params(1+(3*(i-1)),:,:)
       image_mu = params(2+(3*(i-1)),:,:)
       image_sig = params(3+(3*(i-1)),:,:)
       
       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       
       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       
       do l=1, dim_x
          do j=1, dim_y
             !Regularization
             f = f + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
             f = f + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
             f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * ((image_sig(j,l) - b_params(i)) / b_params(i))**2._xp)
             
             ! Regularization for IFS
             if (i .eq. 1) then
                !Corralate sigma fields narrow components
                f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
                !Corralate velocity fields narrow components
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(1-1)),j,l)) - delta_1)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(2-1)),j,l)) - delta_2)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(3-1)),j,l)) - delta_3)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(4-1)),j,l)) - delta_4)**2._xp)

                !Corralate sigma fields broad components
                f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l)) - 1._xp)**2._xp)
                !Corralate velocity fields broad components
                f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(11-1)),j,l) - params(2+(3*(7-1)),j,l)) - delta_1)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(8-1)),j,l)) - delta_2)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(9-1)),j,l)) - delta_3)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(10-1)),j,l)) - delta_4)**2._xp)

                !Corralate amplitude fields narrow and broad components - 1/3 ratio for OIII
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)) - 0.33_xp)**2._xp)
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l)) - 0.33_xp)**2._xp)
                !Corralate amplitude fields narrow and broad components - 0.3 ratio for NII
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)) - 0.3_xp)**2._xp)
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l)) - 0.3_xp)**2._xp)
             elseif (i .ge. 2 .and. i .le. 4) then
                !Corralate sigma fields narrow components
                f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)) - 1._xp)**2._xp)
             elseif (i .ge. 8 .and. i .le. 10) then
                !Corralate sigma fields broad components
                f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l)) - 1._xp)**2._xp)
             end if
             
             ! Gradients ROHSA
             g((n_beta-n_gauss)+i) = g((n_beta-n_gauss)+i) - (lambda_var_sig * (image_sig(j,l) - b_params(i)) * (image_sig(j,l) / b_params(i)**2._xp))
             
             if (ABS(std_map(1,j,l)) > 0._xp) then
                delta = 5._xp * sqrt(params(3+(3*(i-1)),j,l)**2._xp + params_lsf(i)**2._xp)
                ! delta = 5._xp * params(3+(3*(i-1)),j,l)
                start = int(params(2+(3*(i-1)),j,l) - delta)
                finish = int(params(2+(3*(i-1)),j,l) + delta)
                if (start < 1) then 
                   start = 1
                end if
                if (finish > dim_v) then 
                   finish = dim_v
                end if
                do k=start, finish
                ! do k=1, dim_v                          
                   ! deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
                   !      + (dG_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   ! deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
                   !      + (dG_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   ! deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) &
                   !      + (dG_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp))
                   deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
                        + (dGlsf_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
                        + (dGlsf_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) &
                        + (dGlsf_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp))
                end do
             end if

             deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_amp * conv_conv_amp(j,l))
             deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_mu * conv_conv_mu(j,l))
             deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l) + &
                  (lambda_var_sig * (image_sig(j,l) - b_params(i)) / b_params(i)))

             !Gradients IFS
             !Corralate sigma fields narrow and broad components
             if (i .eq. 1) then 
                !Balmer narrow
                deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) - (lambda_sig_corr_narrow * &
                     params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
                     * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                
                
                deriv(3+(3*(1-1)),j,l) = deriv(3+(3*(1-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(5-1)),j,l)&
                     * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
                !Balmer broad
                deriv(3+(3*(11-1)),j,l) = deriv(3+(3*(11-1)),j,l) - (lambda_sig_corr_broad * &
                     params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l)**2._xp & 
                     * (params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))                
                
                deriv(3+(3*(7-1)),j,l) = deriv(3+(3*(7-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(11-1)),j,l)&
                     * (params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))
             elseif (i .ge. 2 .and. i .le. 4) then !Forbidden narrow
                deriv(3+(3*(6-1)),j,l) = deriv(3+(3*(6-1)),j,l) - (lambda_sig_corr_narrow * &
                     params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)**2._xp & 
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))                
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(6-1)),j,l)&
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))
             elseif (i .ge. 8 .and. i .le. 10) then !Forbidden broad
                deriv(3+(3*(12-1)),j,l) = deriv(3+(3*(12-1)),j,l) - (lambda_sig_corr_broad * &
                     params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l)**2._xp & 
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))                
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(12-1)),j,l)&
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))
             end if

             if (i .eq. 1) then             
                ! !Corralate velocity fields narrow and broad components
                !narrow Balmer
                deriv(2+(3*(1-1)),j,l) = deriv(2+(3*(1-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
                     params(2+(3*(1-1)),j,l) - delta_1))
                deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
                     params(2+(3*(1-1)),j,l) - delta_1))

                !narrow Forbidden
                deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(2-1)),j,l) - delta_2))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(2-1)),j,l) - delta_2))

                deriv(2+(3*(3-1)),j,l) = deriv(2+(3*(3-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(3-1)),j,l) - delta_3))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(3-1)),j,l) - delta_3))

                deriv(2+(3*(4-1)),j,l) = deriv(2+(3*(4-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(4-1)),j,l) - delta_4))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(4-1)),j,l) - delta_4))

                !broad Balmer
                deriv(2+(3*(7-1)),j,l) = deriv(2+(3*(7-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(11-1)),j,l) - &
                     params(2+(3*(7-1)),j,l) - delta_1))
                deriv(2+(3*(11-1)),j,l) = deriv(2+(3*(11-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(11-1)),j,l) - &
                     params(2+(3*(7-1)),j,l) - delta_1))

                !broad Forbidden
                deriv(2+(3*(8-1)),j,l) = deriv(2+(3*(8-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(8-1)),j,l) - delta_2))
                deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(8-1)),j,l) - delta_2))

                deriv(2+(3*(9-1)),j,l) = deriv(2+(3*(9-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(9-1)),j,l) - delta_3))
                deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(9-1)),j,l) - delta_3))

                deriv(2+(3*(10-1)),j,l) = deriv(2+(3*(10-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(10-1)),j,l) - delta_4))
                deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
                     params(2+(3*(10-1)),j,l) - delta_4))

                !Corralate amplitude fields narrow and broad components - 1/3 ratio for OIII
                deriv(1+(3*(3-1)),j,l) = deriv(1+(3*(3-1)),j,l) - (lambda_r * &
                     params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)**2._xp & 
                     * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))                
                deriv(1+(3*(2-1)),j,l) = deriv(1+(3*(2-1)),j,l) + (lambda_r / params(1+(3*(3-1)),j,l)&
                     * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))

                deriv(1+(3*(9-1)),j,l) = deriv(1+(3*(9-1)),j,l) - (lambda_r * &
                     params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l)**2._xp & 
                     * (params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l) - 0.33_xp))                
                deriv(1+(3*(8-1)),j,l) = deriv(1+(3*(8-1)),j,l) + (lambda_r / params(1+(3*(9-1)),j,l)&
                     * (params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l) - 0.33_xp))

                !Corralate amplitude fields narrow and broad components - 0.3 ratio for NII
                deriv(1+(3*(6-1)),j,l) = deriv(1+(3*(6-1)),j,l) - (lambda_r * &
                     params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)**2._xp & 
                     * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))                
                deriv(1+(3*(4-1)),j,l) = deriv(1+(3*(4-1)),j,l) + (lambda_r / params(1+(3*(6-1)),j,l)&
                     * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))

                deriv(1+(3*(12-1)),j,l) = deriv(1+(3*(12-1)),j,l) - (lambda_r * &
                     params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l)**2._xp & 
                     * (params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l) - 0.3_xp))                
                deriv(1+(3*(10-1)),j,l) = deriv(1+(3*(10-1)),j,l) + (lambda_r / params(1+(3*(12-1)),j,l)&
                     * (params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l) - 0.3_xp))
             end if

          end do
          !
       end do
    end do        

    call ravel_3D(deriv, g, 3*n_gauss, dim_y, dim_x)

    deallocate(deriv)
    deallocate(residual)
    deallocate(b_params)
    deallocate(params)
    deallocate(conv_amp, conv_mu, conv_sig)
    deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
    deallocate(image_amp, image_mu, image_sig)

  end subroutine f_g_cube_fast_norm

  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_fast_norm_single(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, &
       lambda_mu, lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, lambda_sig_corr_narrow, &
       lambda_sig_corr_broad, lambda_mu_corr_narrow, lambda_mu_corr_broad, lambda_r, params_lsf, &
       delta_1, delta_2, delta_3, delta_4)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    integer, intent(in) :: delta_1, delta_2, delta_3, delta_4
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig, lambda_r
    real(xp), intent(in) :: lambda_sig_corr_narrow, lambda_sig_corr_broad
    real(xp), intent(in) :: lambda_mu_corr_narrow, lambda_mu_corr_broad
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:), allocatable :: params_lsf
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    integer :: n_beta
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:), allocatable :: residual_noise
    real(xp), dimension(:,:,:), allocatable :: params
    real(xp), dimension(:), allocatable :: b_params
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
    real(xp), dimension(:,:,:), allocatable :: deriv
    real(xp), dimension(:), allocatable :: model

    integer :: start, finish
    real(xp) :: delta

    allocate(deriv(3*n_gauss, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(b_params(n_gauss))
    allocate(params(3*n_gauss, dim_y, dim_x))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    allocate(model(dim_v))

    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    residual = 0._xp    
    params = 0._xp
    model = 0._xp
    
    n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
    do i=1,n_gauss
       b_params(i) = beta((n_beta-n_gauss)+i)
    end do

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          allocate(residual_noise(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v, params_lsf)
          residual(:,i,j) = residual_1D
          if (ABS(std_map(1,i,j)) > 0._xp) then
             residual_noise = residual_1D / std_map(:,i,j)
             f = f + (myfunc_spec(residual_noise))
          end if
          deallocate(residual_1D)
          deallocate(residual_noise)
       end do
    end do

    ! Compute the objective function and the gradient
    do i=1, n_gauss
       !
       conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
       conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
       image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       
       image_amp = params(1+(3*(i-1)),:,:)
       image_mu = params(2+(3*(i-1)),:,:)
       image_sig = params(3+(3*(i-1)),:,:)
       
       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       
       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       
       do l=1, dim_x
          do j=1, dim_y
             !Regularization
             f = f + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
             f = f + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
             f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * ((image_sig(j,l) - b_params(i)) / b_params(i))**2._xp)
             
             ! Regularization for IFS
             if (i .eq. 1) then
                !Corralate sigma fields narrow components
                f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
                !Corralate velocity fields narrow components
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(1-1)),j,l)) - delta_1)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(2-1)),j,l)) - delta_2)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(3-1)),j,l)) - delta_3)**2._xp)
                f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(4-1)),j,l)) - delta_4)**2._xp)

                !Corralate amplitude fields narrow components - 1/3 ratio for OIII
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)) - 0.33_xp)**2._xp)
                !Corralate amplitude fields narrow components - 0.3 ratio for NII
                f = f + (0.5_xp * lambda_r * ((params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)) - 0.3_xp)**2._xp)
             elseif (i .ge. 2 .and. i .le. 4) then
                !Corralate sigma fields narrow components
                f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)) - 1._xp)**2._xp)
             end if
             
             ! Gradients ROHSA
             g((n_beta-n_gauss)+i) = g((n_beta-n_gauss)+i) - (lambda_var_sig * (image_sig(j,l) - b_params(i)) * (image_sig(j,l) / b_params(i)**2._xp))
             
             if (ABS(std_map(1,j,l)) > 0._xp) then
                delta = 5._xp * sqrt(params(3+(3*(i-1)),j,l)**2._xp + params_lsf(i)**2._xp)
                ! delta = 5._xp * params(3+(3*(i-1)),j,l)
                start = int(params(2+(3*(i-1)),j,l) - delta)
                finish = int(params(2+(3*(i-1)),j,l) + delta)
                if (start < 0) then 
                   start = 0 
                end if
                if (finish > dim_v) then 
                   finish = dim_v
                end if
                do k=start, finish
                ! do k=1, dim_v                          
                   ! deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
                   !      + (dG_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   ! deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
                   !      + (dG_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   ! deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) &
                   !      + (dG_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
                   !      * (residual(k,j,l)/std_map(k,j,l)**2._xp))
                   deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
                        + (dGlsf_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
                        + (dGlsf_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
                   deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) &
                        + (dGlsf_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
                        * (residual(k,j,l)/std_map(k,j,l)**2._xp))  
                end do
             end if

             deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_amp * conv_conv_amp(j,l))
             deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_mu * conv_conv_mu(j,l))
             deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l) + &
                  (lambda_var_sig * (image_sig(j,l) - b_params(i)) / b_params(i)))

             !Gradients IFS
             !Corralate sigma fields narrow components
             if (i .eq. 1) then 
                !Balmer narrow
                deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) - (lambda_sig_corr_narrow * &
                     params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
                     * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                
                
                deriv(3+(3*(1-1)),j,l) = deriv(3+(3*(1-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(5-1)),j,l)&
                     * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
             elseif (i .ge. 2 .and. i .le. 4) then !Forbidden narrow
                deriv(3+(3*(6-1)),j,l) = deriv(3+(3*(6-1)),j,l) - (lambda_sig_corr_narrow * &
                     params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)**2._xp & 
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))                
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(6-1)),j,l)&
                     * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))
             end if

             if (i .eq. 1) then             
                ! !Corralate velocity fields narrow components
                !narrow Balmer
                deriv(2+(3*(1-1)),j,l) = deriv(2+(3*(1-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
                     params(2+(3*(1-1)),j,l) - delta_1))
                deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
                     params(2+(3*(1-1)),j,l) - delta_1))

                !narrow Forbidden
                deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(2-1)),j,l) - delta_2))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(2-1)),j,l) - delta_2))

                deriv(2+(3*(3-1)),j,l) = deriv(2+(3*(3-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(3-1)),j,l) - delta_3))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(3-1)),j,l) - delta_3))

                deriv(2+(3*(4-1)),j,l) = deriv(2+(3*(4-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(4-1)),j,l) - delta_4))
                deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
                     params(2+(3*(4-1)),j,l) - delta_4))

                !Corralate amplitude fields narrow components - 1/3 ratio for OIII
                deriv(1+(3*(3-1)),j,l) = deriv(1+(3*(3-1)),j,l) - (lambda_r * &
                     params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)**2._xp & 
                     * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))                
                deriv(1+(3*(2-1)),j,l) = deriv(1+(3*(2-1)),j,l) + (lambda_r / params(1+(3*(3-1)),j,l)&
                     * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))

                !Corralate amplitude fields narrow components - 0.3 ratio for NII
                deriv(1+(3*(6-1)),j,l) = deriv(1+(3*(6-1)),j,l) - (lambda_r * &
                     params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)**2._xp & 
                     * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))                
                deriv(1+(3*(4-1)),j,l) = deriv(1+(3*(4-1)),j,l) + (lambda_r / params(1+(3*(6-1)),j,l)&
                     * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))
             end if

          end do
          !
       end do
    end do        

    call ravel_3D(deriv, g, 3*n_gauss, dim_y, dim_x)

    deallocate(deriv)
    deallocate(residual)
    deallocate(b_params)
    deallocate(params)
    deallocate(conv_amp, conv_mu, conv_sig)
    deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
    deallocate(image_amp, image_mu, image_sig)

  end subroutine f_g_cube_fast_norm_single

  ! Objective function to minimize for a spectrum
  pure function  myfunc_spec(residual)
    implicit none
    
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp) :: myfunc_spec
    
    myfunc_spec = 0._xp
    
    myfunc_spec = 0.5_xp * sum(residual**2._xp)    
  end function  myfunc_spec

end module mod_optimize

  ! ! Compute the objective function for a cube and the gradient of the obkective function
  ! subroutine f_g_cube_fast_norm_2(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
  !      lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
  !      lambda_mu_corr_narrow, lambda_mu_corr_broad, lambda_r, params_lsf, delta_1, delta_2, delta_3, delta_4)
  !   implicit none

  !   integer, intent(in) :: n_gauss
  !   integer, intent(in) :: dim_v, dim_y, dim_x
  !   integer, intent(in) :: delta_1, delta_2, delta_3, delta_4
  !   real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig, lambda_r
  !   real(xp), intent(in) :: lambda_sig_corr_narrow, lambda_sig_corr_broad
  !   real(xp), intent(in) :: lambda_mu_corr_narrow, lambda_mu_corr_broad
  !   real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
  !   real(xp), intent(in), dimension(:), allocatable :: beta
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: cube
  !   real(xp), intent(in), dimension(:,:), allocatable :: kernel
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: std_map
  !   real(xp), intent(in), dimension(:), allocatable :: params_lsf
  !   real(xp), intent(inout) :: f
  !   real(xp), intent(inout), dimension(:), allocatable :: g

  !   integer :: i, j, k, l
  !   integer :: n_beta
  !   real(xp), dimension(:,:,:), allocatable :: residual
  !   real(xp), dimension(:), allocatable :: residual_1D
  !   real(xp), dimension(:), allocatable :: residual_noise
  !   real(xp), dimension(:,:,:), allocatable :: params
  !   real(xp), dimension(:), allocatable :: b_params
  !   real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
  !   real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
  !   real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
  !   real(xp), dimension(:,:,:), allocatable :: deriv
  !   real(xp), dimension(:), allocatable :: model

  !   allocate(deriv(3*n_gauss, dim_y, dim_x))
  !   allocate(residual(dim_v, dim_y, dim_x))
  !   allocate(b_params(n_gauss))
  !   allocate(params(3*n_gauss, dim_y, dim_x))
  !   allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
  !   allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
  !   allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
  !   allocate(model(dim_v))

  !   deriv = 0._xp
  !   f = 0._xp
  !   g = 0._xp
  !   residual = 0._xp    
  !   params = 0._xp
  !   model = 0._xp
    
  !   n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

  !   call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
  !   do i=1,n_gauss
  !      b_params(i) = beta((n_beta-n_gauss)+i)
  !   end do

  !   ! Compute the objective function and the gradient
  !   do j=1, dim_x
  !      do i=1, dim_y
  !         allocate(residual_1D(dim_v))
  !         allocate(residual_noise(dim_v))
  !         residual_1D = 0._xp
  !         call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v, params_lsf)
  !         residual(:,i,j) = residual_1D
  !         if (ABS(std_map(1,i,j)) > 0._xp) then
  !            residual_noise = residual_1D / std_map(:,i,j)
  !            f = f + (myfunc_spec(residual_noise))
  !         end if
  !         deallocate(residual_1D)
  !         deallocate(residual_noise)
  !      end do
  !   end do

  !   ! Compute the objective function and the gradient
  !   do i=1, n_gauss
  !      !
  !      conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
  !      conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
  !      image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       
  !      image_amp = params(1+(3*(i-1)),:,:)
  !      image_mu = params(2+(3*(i-1)),:,:)
  !      image_sig = params(3+(3*(i-1)),:,:)
       
  !      call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       
  !      call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       
  !      do l=1, dim_x
  !         do j=1, dim_y
  !            !Regularization
  !            f = f + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
  !            f = f + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
  !            f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * ((image_sig(j,l) - b_params(i)) / b_params(i))**2._xp)             
  !            ! Gradients ROHSA
  !            g((n_beta-n_gauss)+i) = g((n_beta-n_gauss)+i) - (lambda_var_sig * (image_sig(j,l) - b_params(i)) * (image_sig(j,l) / b_params(i)**2._xp))
             
  !            if (ABS(std_map(1,j,l)) > 0._xp) then
  !               do k=1, dim_v                          
  !                  deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
  !                       + (dGlsf_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
  !                       + (dGlsf_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + &
  !                       + (dGlsf_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp))
  !               end do
  !            end if

  !            deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_amp * conv_conv_amp(j,l))
  !            deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_mu * conv_conv_mu(j,l))
  !            deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l) + &
  !                 (lambda_var_sig * (image_sig(j,l) - b_params(i)) / b_params(i)))
  !         end do
  !         !
  !      end do
  !   end do 
    
  !   ! Regularization for IFS
  !   do l=1, dim_x
  !      do j=1, dim_y
  !         !Corralate sigma fields narrow components
  !         f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
  !         !Corralate velocity fields narrow components
  !         f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(1-1)),j,l)) - delta_1)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(2-1)),j,l)) - delta_2)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(3-1)),j,l)) - delta_3)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(6-1)),j,l) - params(2+(3*(4-1)),j,l)) - delta_4)**2._xp)

  !         !Corralate sigma fields broad components
  !         f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l)) - 1._xp)**2._xp)
  !         !Corralate velocity fields broad components
  !         f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(11-1)),j,l) - params(2+(3*(7-1)),j,l)) - delta_1)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(8-1)),j,l)) - delta_2)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(9-1)),j,l)) - delta_3)**2._xp)
  !         f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(12-1)),j,l) - params(2+(3*(10-1)),j,l)) - delta_4)**2._xp)

  !         !Corralate amplitude fields narrow and broad components - 1/3 ratio for OIII
  !         f = f + (0.5_xp * lambda_r * ((params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)) - 0.33_xp)**2._xp)
  !         f = f + (0.5_xp * lambda_r * ((params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l)) - 0.33_xp)**2._xp)
  !         !Corralate amplitude fields narrow and broad components - 0.3 ratio for NII
  !         f = f + (0.5_xp * lambda_r * ((params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)) - 0.3_xp)**2._xp)
  !         f = f + (0.5_xp * lambda_r * ((params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l)) - 0.3_xp)**2._xp)

  !         !Corralate sigma fields narrow and broad components
  !         !Balmer narrow
  !         deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) - (lambda_sig_corr_narrow * &
  !              params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
  !              * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                

  !         deriv(3+(3*(1-1)),j,l) = deriv(3+(3*(1-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(5-1)),j,l)&
  !              * (params(3+(3*(1-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
  !         !Balmer broad
  !         deriv(3+(3*(11-1)),j,l) = deriv(3+(3*(11-1)),j,l) - (lambda_sig_corr_broad * &
  !              params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l)**2._xp & 
  !              * (params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))                

  !         deriv(3+(3*(7-1)),j,l) = deriv(3+(3*(7-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(11-1)),j,l)&
  !              * (params(3+(3*(7-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))

  !         ! !Corralate velocity fields narrow and broad components
  !         !narrow Balmer
  !         deriv(2+(3*(1-1)),j,l) = deriv(2+(3*(1-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
  !              params(2+(3*(1-1)),j,l) - delta_1))
  !         deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(5-1)),j,l) - &
  !              params(2+(3*(1-1)),j,l) - delta_1))

  !         !narrow Forbidden
  !         deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(2-1)),j,l) - delta_2))
  !         deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(2-1)),j,l) - delta_2))

  !         deriv(2+(3*(3-1)),j,l) = deriv(2+(3*(3-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(3-1)),j,l) - delta_3))
  !         deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(3-1)),j,l) - delta_3))

  !         deriv(2+(3*(4-1)),j,l) = deriv(2+(3*(4-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(4-1)),j,l) - delta_4))
  !         deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(6-1)),j,l) - &
  !              params(2+(3*(4-1)),j,l) - delta_4))

  !         !broad Balmer
  !         deriv(2+(3*(7-1)),j,l) = deriv(2+(3*(7-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(11-1)),j,l) - &
  !              params(2+(3*(7-1)),j,l) - delta_1))
  !         deriv(2+(3*(11-1)),j,l) = deriv(2+(3*(11-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(11-1)),j,l) - &
  !              params(2+(3*(7-1)),j,l) - delta_1))

  !         !broad Forbidden
  !         deriv(2+(3*(8-1)),j,l) = deriv(2+(3*(8-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(8-1)),j,l) - delta_2))
  !         deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(8-1)),j,l) - delta_2))

  !         deriv(2+(3*(9-1)),j,l) = deriv(2+(3*(9-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(9-1)),j,l) - delta_3))
  !         deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(9-1)),j,l) - delta_3))

  !         deriv(2+(3*(10-1)),j,l) = deriv(2+(3*(10-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(10-1)),j,l) - delta_4))
  !         deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(12-1)),j,l) - &
  !              params(2+(3*(10-1)),j,l) - delta_4))

  !         !Corralate amplitude fields narrow and broad components - 1/3 ratio for OIII
  !         deriv(1+(3*(3-1)),j,l) = deriv(1+(3*(3-1)),j,l) - (lambda_r * &
  !              params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l)**2._xp & 
  !              * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))                
  !         deriv(1+(3*(2-1)),j,l) = deriv(1+(3*(2-1)),j,l) + (lambda_r / params(1+(3*(3-1)),j,l)&
  !              * (params(1+(3*(2-1)),j,l) / params(1+(3*(3-1)),j,l) - 0.33_xp))

  !         deriv(1+(3*(9-1)),j,l) = deriv(1+(3*(9-1)),j,l) - (lambda_r * &
  !              params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l)**2._xp & 
  !              * (params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l) - 0.33_xp))                
  !         deriv(1+(3*(8-1)),j,l) = deriv(1+(3*(8-1)),j,l) + (lambda_r / params(1+(3*(9-1)),j,l)&
  !              * (params(1+(3*(8-1)),j,l) / params(1+(3*(9-1)),j,l) - 0.33_xp))

  !         !Corralate amplitude fields narrow and broad components - 0.3 ratio for NII
  !         deriv(1+(3*(6-1)),j,l) = deriv(1+(3*(6-1)),j,l) - (lambda_r * &
  !              params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)**2._xp & 
  !              * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))                
  !         deriv(1+(3*(4-1)),j,l) = deriv(1+(3*(4-1)),j,l) + (lambda_r / params(1+(3*(6-1)),j,l)&
  !              * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))

  !         deriv(1+(3*(12-1)),j,l) = deriv(1+(3*(12-1)),j,l) - (lambda_r * &
  !              params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l)**2._xp & 
  !              * (params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l) - 0.3_xp))                
  !         deriv(1+(3*(10-1)),j,l) = deriv(1+(3*(10-1)),j,l) + (lambda_r / params(1+(3*(12-1)),j,l)&
  !              * (params(1+(3*(10-1)),j,l) / params(1+(3*(12-1)),j,l) - 0.3_xp))
  !      end do
  !   end do

  !   do i=2,4
  !      do l=1, dim_x
  !         do j=1, dim_y
  !            !Corralate sigma fields narrow components
  !            f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)) - 1._xp)**2._xp)

  !               !Forbidden narrow
  !               deriv(3+(3*(6-1)),j,l) = deriv(3+(3*(6-1)),j,l) - (lambda_sig_corr_narrow * &
  !                    params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l)**2._xp & 
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))                
  !               deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(6-1)),j,l)&
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))
  !         end do
  !      end do
  !   end do
    
  !   do i=8,10
  !      do l=1, dim_x
  !         do j=1, dim_y
  !            !Corralate sigma fields broad components
  !            f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l)) - 1._xp)**2._xp)

  !            !Forbidden broad
  !            deriv(3+(3*(12-1)),j,l) = deriv(3+(3*(12-1)),j,l) - (lambda_sig_corr_broad * &
  !                 params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l)**2._xp & 
  !                 * (params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))                
  !            deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(12-1)),j,l)&
  !                 * (params(3+(3*(i-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))
  !         end do
  !      end do
  !   end do
    
  !   call ravel_3D(deriv, g, 3*n_gauss, dim_y, dim_x)

  !   deallocate(deriv)
  !   deallocate(residual)
  !   deallocate(b_params)
  !   deallocate(params)
  !   deallocate(conv_amp, conv_mu, conv_sig)
  !   deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
  !   deallocate(image_amp, image_mu, image_sig)

  ! end subroutine f_g_cube_fast_norm_2


  ! ! Compute the objective function for a cube and the gradient of the obkective function
  ! subroutine f_g_cube_fast_norm(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
  !      lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
  !      lambda_mu_corr_narrow, lambda_mu_corr_broad, delta_21, delta_32, lambda_r, delta_27, delta_28, delta_29, &
  !      params_lsf)
  !   implicit none

  !   integer, intent(in) :: n_gauss
  !   integer, intent(in) :: dim_v, dim_y, dim_x
  !   integer, intent(in) :: delta_21, delta_32, delta_27, delta_28, delta_29
  !   real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig, lambda_r
  !   real(xp), intent(in) :: lambda_sig_corr_narrow, lambda_sig_corr_broad
  !   real(xp), intent(in) :: lambda_mu_corr_narrow, lambda_mu_corr_broad
  !   real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
  !   real(xp), intent(in), dimension(:), allocatable :: beta
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: cube
  !   real(xp), intent(in), dimension(:,:), allocatable :: kernel
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: std_map
  !   real(xp), intent(in), dimension(:), allocatable :: params_lsf
  !   real(xp), intent(inout) :: f
  !   real(xp), intent(inout), dimension(:), allocatable :: g

  !   integer :: i, j, k, l
  !   integer :: n_beta
  !   real(xp), dimension(:,:,:), allocatable :: residual
  !   real(xp), dimension(:), allocatable :: residual_1D
  !   real(xp), dimension(:), allocatable :: residual_noise
  !   real(xp), dimension(:,:,:), allocatable :: params
  !   real(xp), dimension(:), allocatable :: b_params
  !   real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
  !   real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
  !   real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
  !   real(xp), dimension(:,:,:), allocatable :: deriv
  !   real(xp), dimension(:), allocatable :: model
  !   ! real(xp) :: gauss

  !   ! real(xp), dimension(:), allocatable :: lsf
  !   ! allocate(lsf(n_gauss))
  !   ! lsf = 2._xp

  !   allocate(deriv(3*n_gauss, dim_y, dim_x))
  !   allocate(residual(dim_v, dim_y, dim_x))
  !   allocate(b_params(n_gauss))
  !   allocate(params(3*n_gauss, dim_y, dim_x))
  !   allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
  !   allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
  !   allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
  !   allocate(model(dim_v))

  !   deriv = 0._xp
  !   f = 0._xp
  !   g = 0._xp
  !   residual = 0._xp    
  !   params = 0._xp
  !   model = 0._xp
  !   ! gauss = 0._xp
    
  !   n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

  !   call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
  !   do i=1,n_gauss
  !      b_params(i) = beta((n_beta-n_gauss)+i)
  !   end do

  !   ! Compute the objective function and the gradient
  !   do j=1, dim_x
  !      do i=1, dim_y
  !         allocate(residual_1D(dim_v))
  !         allocate(residual_noise(dim_v))
  !         residual_1D = 0._xp
  !         call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v, params_lsf)
  !         residual(:,i,j) = residual_1D
  !         if (ABS(std_map(1,i,j)) > 0._xp) then
  !            residual_noise = residual_1D / std_map(:,i,j)
  !            f = f + (myfunc_spec(residual_noise))
  !         end if
  !         deallocate(residual_1D)
  !         deallocate(residual_noise)
  !      end do
  !   end do

  !   ! Compute the objective function and the gradient
  !   do i=1, n_gauss
  !      !
  !      conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
  !      conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
  !      image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       
  !      image_amp = params(1+(3*(i-1)),:,:)
  !      image_mu = params(2+(3*(i-1)),:,:)
  !      image_sig = params(3+(3*(i-1)),:,:)
       
  !      call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       
  !      call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       
  !      do l=1, dim_x
  !         do j=1, dim_y
  !            !Regularization
  !            f = f + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
  !            f = f + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
  !            f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * ((image_sig(j,l) - b_params(i)) / b_params(i))**2._xp)
             
  !            ! Regularization for IFS
  !            !Corralate sigma fields narrow and broad components
  !            if (i .lt. 4) then
  !               f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)) - 1._xp)**2._xp)
  !            elseif (i .ge. 4 .and. i .lt. 7) then
  !               f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
  !            elseif (i .ge. 7 .and. i .lt. 10) then
  !               f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)) - 1._xp)**2._xp)
  !            elseif (i .ge. 10) then
  !               f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
  !            end if

  !            if (i .eq. 1) then !so that it is calculated one (no i index)
  !               !Corralate velocity fields narrow and broad components
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(1-1)),j,l) - params(2+(3*(2-1)),j,l)) + delta_21)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(3-1)),j,l) - params(2+(3*(2-1)),j,l)) - delta_32)**2._xp)

  !               f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(4-1)),j,l) - params(2+(3*(5-1)),j,l)) + delta_21)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(6-1)),j,l) - params(2+(3*(5-1)),j,l)) - delta_32)**2._xp)

  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(2-1)),j,l) - params(2+(3*(7-1)),j,l)) - delta_27)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(2-1)),j,l) - params(2+(3*(8-1)),j,l)) - delta_28)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(2-1)),j,l) - params(2+(3*(9-1)),j,l)) - delta_29)**2._xp)

  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(10-1)),j,l)) - delta_27)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(11-1)),j,l)) - delta_28)**2._xp)
  !               f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(5-1)),j,l) - params(2+(3*(12-1)),j,l)) - delta_29)**2._xp)

  !               !Corralate amplitude fields narrow and broad components - 1/3 ratio
  !               f = f + (0.5_xp * lambda_r * ((params(1,j,l) / params(7,j,l)) - 0.3_xp)**2._xp)
  !               f = f + (0.5_xp * lambda_r * ((params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)) - 0.3_xp)**2._xp)

  !               ! !Corralate sigma OIII lines broad
  !               ! f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(10-1)),j,l) / params(3+(3*(12-1)),j,l)) - 1._xp)**2._xp)
  !               ! !Corralate sigma OIII lines narrow
  !               ! f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(7-1)),j,l) / params(3+(3*(9-1)),j,l)) - 1._xp)**2._xp)     

  !               ! !Corralate sigma NII lines narrow
  !               ! f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(1-1)),j,l) / params(3+(3*(3-1)),j,l)) - 1._xp)**2._xp)
  !               ! !Corralate sigma NII lines broad
  !               ! f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(4-1)),j,l) / params(3+(3*(6-1)),j,l)) - 1._xp)**2._xp)     

  !               ! !Corralate sigma Hab lines narrow
  !               ! f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(2-1)),j,l) / params(3+(3*(8-1)),j,l)) - 1._xp)**2._xp)
  !               ! !Corralate sigma Hab lines broad
  !               ! f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(5-1)),j,l) / params(3+(3*(11-1)),j,l)) - 1._xp)**2._xp)     
  !            end if
             
  !            ! Gradients ROHSA
  !            g((n_beta-n_gauss)+i) = g((n_beta-n_gauss)+i) - (lambda_var_sig * (image_sig(j,l) - b_params(i)) * (image_sig(j,l) / b_params(i)**2._xp))
             
  !            if (ABS(std_map(1,j,l)) > 0._xp) then
  !               do k=1, dim_v                          
  !                  deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
  !                       + (dGlsf_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
  !                       + (dGlsf_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + &
  !                       + (dGlsf_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l), params_lsf(i)) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp))

  !                  ! deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) &
  !                  !      + (dG_da(k, params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
  !                  !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  ! deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) &
  !                  !      + (dG_dmu(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
  !                  !      * (residual(k,j,l)/std_map(k,j,l)**2._xp)) 
                   
  !                  ! deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + &
  !                  !      + (dG_dsig(k, params(1+(3*(i-1)),j,l), params(2+(3*(i-1)),j,l), params(3+(3*(i-1)),j,l)) &
  !                  !      * (residual(k,j,l)/std_map(k,j,l)**2._xp))
  !               end do
  !            end if

  !            deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_amp * conv_conv_amp(j,l))
  !            deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_mu * conv_conv_mu(j,l))
  !            deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l) + &
  !                 (lambda_var_sig * (image_sig(j,l) - b_params(i)) / b_params(i)))

  !            !Gradients IFS
  !            !Corralate sigma fields narrow and broad components
  !            if (i .lt. 4) then
  !               deriv(3+(3*(2-1)),j,l) = deriv(3+(3*(2-1)),j,l) - (lambda_sig_corr_narrow * &
  !                    params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)**2._xp & 
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))                
                
  !               deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(2-1)),j,l)&
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))
  !            elseif (i .ge. 4 .and. i .lt. 7) then
  !               deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) - (lambda_sig_corr_broad * &
  !                    params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                
  !               deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(5-1)),j,l)&
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
  !            elseif (i .ge. 7 .and. i .lt. 10) then
  !               deriv(3+(3*(2-1)),j,l) = deriv(3+(3*(2-1)),j,l) - (lambda_sig_corr_narrow * &
  !                    params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)**2._xp & 
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))                
                
  !               deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(2-1)),j,l)&
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))
  !            elseif (i .ge. 10) then
  !               deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) - (lambda_sig_corr_broad * &
  !                    params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                
  !               deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(5-1)),j,l)&
  !                    * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
  !            end if

  !            if (i .eq. 1) then             
  !               !Corralate velocity fields narrow and broad components
  !               !narrow
  !               deriv(2+(3*(1-1)),j,l) = deriv(2+(3*(1-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(1-1)),j,l) - &
  !                    params(2+(3*(2-1)),j,l) + delta_21))
  !               deriv(2+(3*(3-1)),j,l) = deriv(2+(3*(3-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(3-1)),j,l) - &
  !                    params(2+(3*(2-1)),j,l) - delta_32))
  !               deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) + (lambda_mu_corr_narrow * (2._xp*params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(1-1)),j,l) - params(2+(3*(3-1)),j,l) - delta_21 + delta_32))

  !               !broad
  !               deriv(2+(3*(4-1)),j,l) = deriv(2+(3*(4-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(4-1)),j,l) - &
  !                    params(2+(3*(5-1)),j,l) + delta_21))
  !               deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(6-1)),j,l) - &
  !                    params(2+(3*(5-1)),j,l) - delta_32))
  !               deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_broad * (2._xp*params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(4-1)),j,l) - params(2+(3*(6-1)),j,l) - delta_21 + delta_32))

  !               !narrow
  !               deriv(2+(3*(7-1)),j,l) = deriv(2+(3*(7-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(7-1)),j,l) - delta_27))
  !               deriv(2+(3*(8-1)),j,l) = deriv(2+(3*(8-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(8-1)),j,l) - delta_28))
  !               deriv(2+(3*(9-1)),j,l) = deriv(2+(3*(9-1)),j,l) - (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(9-1)),j,l) - delta_29))

  !               deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(7-1)),j,l) - delta_27))
  !               deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(8-1)),j,l) - delta_28))
  !               deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(2-1)),j,l) - &
  !                    params(2+(3*(9-1)),j,l) - delta_29))

  !               !broad
  !               deriv(2+(3*(10-1)),j,l) = deriv(2+(3*(10-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(10-1)),j,l) - delta_27))
  !               deriv(2+(3*(11-1)),j,l) = deriv(2+(3*(11-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(11-1)),j,l) - delta_28))
  !               deriv(2+(3*(12-1)),j,l) = deriv(2+(3*(12-1)),j,l) - (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(12-1)),j,l) - delta_29))

  !               deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(10-1)),j,l) - delta_27))
  !               deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(11-1)),j,l) - delta_28))
  !               deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(5-1)),j,l) - &
  !                    params(2+(3*(12-1)),j,l) - delta_29))

  !               !Corralate amplitude fields narrow and broad components
  !               deriv(7,j,l) = deriv(7,j,l) - (lambda_r * params(1,j,l) / params(7,j,l)**2._xp & 
  !                    * (params(1,j,l) / params(7,j,l) - 0.3_xp))                
  !               deriv(1,j,l) = deriv(1,j,l) + (lambda_r / params(7,j,l) &
  !                    * (params(1,j,l) / params(7,j,l) - 0.3_xp))

  !               deriv(1+(3*(6-1)),j,l) = deriv(1+(3*(6-1)),j,l) - (lambda_r * &
  !                    params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)**2._xp & 
  !                    * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))                
  !               deriv(1+(3*(4-1)),j,l) = deriv(1+(3*(4-1)),j,l) + (lambda_r / params(1+(3*(6-1)),j,l)&
  !                    * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))

  !               ! !Correlate sigma OIII lines broad
  !               ! deriv(3+(3*(12-1)),j,l) = deriv(3+(3*(12-1)),j,l) - (lambda_sig_corr_broad * &
  !               !      params(3+(3*(10-1)),j,l) / params(3+(3*(12-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(10-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(10-1)),j,l) = deriv(3+(3*(10-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(12-1)),j,l)&
  !               !      * (params(3+(3*(10-1)),j,l) / params(3+(3*(12-1)),j,l) - 1._xp))

  !               ! !Correlate sigma OIII lines narrow
  !               ! deriv(3+(3*(9-1)),j,l) = deriv(3+(3*(9-1)),j,l) - (lambda_sig_corr_narrow * &
  !               !      params(3+(3*(7-1)),j,l) / params(3+(3*(9-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(7-1)),j,l) / params(3+(3*(9-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(7-1)),j,l) = deriv(3+(3*(7-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(9-1)),j,l)&
  !               !      * (params(3+(3*(7-1)),j,l) / params(3+(3*(9-1)),j,l) - 1._xp))

  !               ! !Correlate sigma NII lines narrow
  !               ! deriv(3+(3*(3-1)),j,l) = deriv(3+(3*(3-1)),j,l) - (lambda_sig_corr_narrow * &
  !               !      params(3+(3*(1-1)),j,l) / params(3+(3*(3-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(1-1)),j,l) / params(3+(3*(3-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(1-1)),j,l) = deriv(3+(3*(1-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(3-1)),j,l)&
  !               !      * (params(3+(3*(1-1)),j,l) / params(3+(3*(3-1)),j,l) - 1._xp))

  !               ! !Correlate sigma NII lines broad
  !               ! deriv(3+(3*(6-1)),j,l) = deriv(3+(3*(6-1)),j,l) - (lambda_sig_corr_broad * &
  !               !      params(3+(3*(4-1)),j,l) / params(3+(3*(6-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(4-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(4-1)),j,l) = deriv(3+(3*(4-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(6-1)),j,l)&
  !               !      * (params(3+(3*(4-1)),j,l) / params(3+(3*(6-1)),j,l) - 1._xp))

  !               ! !Correlate sigma Hab lines narrow
  !               ! deriv(3+(3*(8-1)),j,l) = deriv(3+(3*(8-1)),j,l) - (lambda_sig_corr_narrow * &
  !               !      params(3+(3*(2-1)),j,l) / params(3+(3*(8-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(2-1)),j,l) / params(3+(3*(8-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(2-1)),j,l) = deriv(3+(3*(2-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(8-1)),j,l)&
  !               !      * (params(3+(3*(2-1)),j,l) / params(3+(3*(8-1)),j,l) - 1._xp))

  !               ! !Correlate sigma Hab lines broad
  !               ! deriv(3+(3*(11-1)),j,l) = deriv(3+(3*(11-1)),j,l) - (lambda_sig_corr_broad * &
  !               !      params(3+(3*(5-1)),j,l) / params(3+(3*(11-1)),j,l)**2._xp & 
  !               !      * (params(3+(3*(5-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))                
  !               ! deriv(3+(3*(5-1)),j,l) = deriv(3+(3*(5-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(11-1)),j,l)&
  !               !      * (params(3+(3*(5-1)),j,l) / params(3+(3*(11-1)),j,l) - 1._xp))
  !            end if

  !         end do
  !         !
  !      end do
  !   end do        

  !   call ravel_3D(deriv, g, 3*n_gauss, dim_y, dim_x)

  !   deallocate(deriv)
  !   deallocate(residual)
  !   deallocate(b_params)
  !   deallocate(params)
  !   deallocate(conv_amp, conv_mu, conv_sig)
  !   deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
  !   deallocate(image_amp, image_mu, image_sig)

  ! end subroutine f_g_cube_fast_norm

  ! ! Compute the objective function for a cube and the gradient of the obkective function
  ! subroutine f_g_cube_fast_norm(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
  !      lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, lambda_sig_corr_narrow, lambda_sig_corr_broad, &
  !      lambda_mu_corr_narrow, lambda_mu_corr_broad, delta_21, delta_32, lambda_r)
  !   implicit none

  !   integer, intent(in) :: n_gauss
  !   integer, intent(in) :: dim_v, dim_y, dim_x
  !   integer, intent(in) :: delta_21, delta_32
  !   real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig, lambda_r
  !   real(xp), intent(in) :: lambda_sig_corr_narrow, lambda_sig_corr_broad
  !   real(xp), intent(in) :: lambda_mu_corr_narrow, lambda_mu_corr_broad
  !   real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
  !   real(xp), intent(in), dimension(:), allocatable :: beta
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: cube
  !   real(xp), intent(in), dimension(:,:), allocatable :: kernel
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: std_map
  !   real(xp), intent(inout) :: f
  !   real(xp), intent(inout), dimension(:), allocatable :: g

  !   integer :: i, j, k, l
  !   integer :: n_beta
  !   real(xp), dimension(:,:,:), allocatable :: residual
  !   real(xp), dimension(:), allocatable :: residual_1D
  !   real(xp), dimension(:), allocatable :: residual_noise
  !   real(xp), dimension(:,:,:), allocatable :: params
  !   real(xp), dimension(:), allocatable :: b_params
  !   real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
  !   real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
  !   real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
  !   real(xp), dimension(:,:,:), allocatable :: deriv
  !   real(xp), dimension(:), allocatable :: model
  !   real(xp) :: gauss

  !   allocate(deriv(3*n_gauss, dim_y, dim_x))
  !   allocate(residual(dim_v, dim_y, dim_x))
  !   allocate(b_params(n_gauss))
  !   allocate(params(3*n_gauss, dim_y, dim_x))
  !   allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
  !   allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
  !   allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
  !   allocate(model(dim_v))

  !   deriv = 0._xp
  !   f = 0._xp
  !   g = 0._xp
  !   residual = 0._xp    
  !   params = 0._xp
  !   model = 0._xp
  !   gauss = 0._xp
    
  !   n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

  !   call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
  !   do i=1,n_gauss
  !      b_params(i) = beta((n_beta-n_gauss)+i)
  !   end do

  !   ! Compute the objective function and the gradient
  !   do j=1, dim_x
  !      do i=1, dim_y
  !         allocate(residual_1D(dim_v))
  !         allocate(residual_noise(dim_v))
  !         residual_1D = 0._xp
  !         call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v)
  !         residual(:,i,j) = residual_1D
  !         if (ABS(std_map(1,i,j)) > 0._xp) then
  !            residual_noise = residual_1D / std_map(:,i,j)
  !            f = f + (myfunc_spec(residual_noise))
  !         end if
  !         deallocate(residual_1D)
  !         deallocate(residual_noise)
  !      end do
  !   end do

  !   ! Compute the objective function and the gradient
  !   do i=1, n_gauss
  !      !
  !      conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
  !      conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
  !      image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       
  !      image_amp = params(1+(3*(i-1)),:,:)
  !      image_mu = params(2+(3*(i-1)),:,:)
  !      image_sig = params(3+(3*(i-1)),:,:)
       
  !      call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       
  !      call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       
  !      do l=1, dim_x
  !         do j=1, dim_y
  !            !Regularization
  !            f = f + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
  !            f = f + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
  !            f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * ((image_sig(j,l) - b_params(i)) / b_params(i))**2._xp)
             
  !            !Regularization for IFS
  !            !Corralate sigma fields narrow and broad components
  !            if (i .lt. 4) then
  !               f = f + (0.5_xp * lambda_sig_corr_narrow * ((params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)) - 1._xp)**2._xp)
  !            else
  !               f = f + (0.5_xp * lambda_sig_corr_broad * ((params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)) - 1._xp)**2._xp)
  !            end if             

  !            !Corralate velocity fields narrow and broad components
  !            f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(1-1)),j,l) - params(2+(3*(2-1)),j,l)) + delta_21)**2._xp)
  !            f = f + (0.5_xp * lambda_mu_corr_narrow * ((params(2+(3*(3-1)),j,l) - params(2+(3*(2-1)),j,l)) - delta_32)**2._xp)

  !            f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(4-1)),j,l) - params(2+(3*(5-1)),j,l)) + delta_21)**2._xp)
  !            f = f + (0.5_xp * lambda_mu_corr_broad * ((params(2+(3*(6-1)),j,l) - params(2+(3*(5-1)),j,l)) - delta_32)**2._xp)

  !            !Corralate amplitude fields narrow and broad components - 1/3 ratio
  !            f = f + (0.5_xp * lambda_r * ((params(1,j,l) / params(7,j,l)) - 0.3_xp)**2._xp)
  !            f = f + (0.5_xp * lambda_r * ((params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)) - 0.3_xp)**2._xp)
             
  !            !Gradients ROHSA
  !            g((n_beta-n_gauss)+i) = g((n_beta-n_gauss)+i) - (lambda_var_sig * (image_sig(j,l) - b_params(i)) * (image_sig(j,l) / b_params(i)**2._xp))
             
  !            !
  !            if (ABS(std_map(1,j,l)) > 0._xp) then
  !               do k=1, dim_v                          
  !                  deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) &
  !                       / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp) 
                   
  !                  deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params(1+(3*(i-1)),j,l) * &
  !                       ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2._xp) * &
  !                       exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) &
  !                       / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp) 
                   
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params(1+(3*(i-1)),j,l) * &
  !                       ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2._xp / (params(3+(3*(i-1)),j,l)**3._xp) * &
  !                       exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
  !                       * (residual(k,j,l)/std_map(k,j,l)**2._xp)
  !               end do
  !            end if

  !            deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_amp * conv_conv_amp(j,l))
  !            deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_mu * conv_conv_mu(j,l))
  !            deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l) + &
  !                 (lambda_var_sig * (image_sig(j,l) - b_params(i)) / b_params(i)))

  !            !Gradients IFS
  !            !Corralate sigma fields narrow and broad components
  !            if (i .lt. 4) then
  !               if (i .eq. 2) then
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) - (lambda_sig_corr_narrow * &
  !                       params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l)**2._xp & 
  !                       * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))                
  !               else
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_narrow / params(3+(3*(2-1)),j,l)&
  !                       * (params(3+(3*(i-1)),j,l) / params(3+(3*(2-1)),j,l) - 1._xp))
  !               end if
  !            else
  !               if (i .eq. 5) then
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) - (lambda_sig_corr_broad * &
  !                       params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l)**2._xp & 
  !                       * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))                
  !               else
  !                  deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_sig_corr_broad / params(3+(3*(5-1)),j,l)&
  !                       * (params(3+(3*(i-1)),j,l) / params(3+(3*(5-1)),j,l) - 1._xp))
  !               end if
  !            end if

  !            !Corralate velocity fields narrow and broad components
  !            deriv(2+(3*(1-1)),j,l) = deriv(2+(3*(1-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(1-1)),j,l) - &
  !                 params(2+(3*(2-1)),j,l) + delta_21))
  !            deriv(2+(3*(3-1)),j,l) = deriv(2+(3*(3-1)),j,l) + (lambda_mu_corr_narrow * (params(2+(3*(3-1)),j,l) - &
  !                 params(2+(3*(2-1)),j,l) - delta_32))
  !            deriv(2+(3*(2-1)),j,l) = deriv(2+(3*(2-1)),j,l) + (lambda_mu_corr_narrow * (2._xp*params(2+(3*(2-1)),j,l) - &
  !                 params(2+(3*(1-1)),j,l) - params(2+(3*(3-1)),j,l) - delta_21 + delta_32))

  !            deriv(2+(3*(4-1)),j,l) = deriv(2+(3*(4-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(4-1)),j,l) - &
  !                 params(2+(3*(5-1)),j,l) + delta_21))
  !            deriv(2+(3*(6-1)),j,l) = deriv(2+(3*(6-1)),j,l) + (lambda_mu_corr_broad * (params(2+(3*(6-1)),j,l) - &
  !                 params(2+(3*(5-1)),j,l) - delta_32))
  !            deriv(2+(3*(5-1)),j,l) = deriv(2+(3*(5-1)),j,l) + (lambda_mu_corr_broad * (2._xp*params(2+(3*(5-1)),j,l) - &
  !                 params(2+(3*(4-1)),j,l) - params(2+(3*(6-1)),j,l) - delta_21 + delta_32))

  !            !Corralate amplitude fields narrow and broad components
  !            deriv(7,j,l) = deriv(7,j,l) - (lambda_r * params(1,j,l) / params(7,j,l)**2._xp & 
  !                 * (params(1,j,l) / params(7,j,l) - 0.3_xp))                
  !            deriv(1,j,l) = deriv(1,j,l) + (lambda_r / params(7,j,l) &
  !                 * (params(1,j,l) / params(7,j,l) - 0.3_xp))

  !            deriv(1+(3*(6-1)),j,l) = deriv(1+(3*(6-1)),j,l) - (lambda_r * &
  !                 params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l)**2._xp & 
  !                 * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))                
  !            deriv(1+(3*(4-1)),j,l) = deriv(1+(3*(4-1)),j,l) + (lambda_r / params(1+(3*(6-1)),j,l)&
  !                 * (params(1+(3*(4-1)),j,l) / params(1+(3*(6-1)),j,l) - 0.3_xp))

  !         end do
  !         !
  !      end do
  !   end do        
    
  !   call ravel_3D(deriv, g, 3*n_gauss, dim_y, dim_x)

  !   deallocate(deriv)
  !   deallocate(residual)
  !   deallocate(b_params)
  !   deallocate(params)
  !   deallocate(conv_amp, conv_mu, conv_sig)
  !   deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
  !   deallocate(image_amp, image_mu, image_sig)

  ! end subroutine f_g_cube_fast_norm
