!< Currently the code is being written with the following equation in mind
!< \tau_\theta \frac{\partial \theta}{\partial t} = 
!<          \nabla \cdot \left ( 2h(\eta) \nabla \theta + g(\eta) \frac{\nabla \theta}{\left | \nabla \theta \right |} \right ) 
!< This is the equation 40 mentioned in https://doi.org/10.1080/14786435.2012.713135
!< The orientation is considered as a scalar quantity, 
!< in reality it is quaternion. 

module grid_orientation_spectral

!public :: &
 real(pReal), dimension(3,3), private :: D_ref
 
contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_init
 integer :: i, j, k, cell
 
 cell = 0
 D_ref = 0.0_pReal
 mobility_ref = 0.0_pReal
 
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
   D_ref = D_ref + grid_orientation_calculate_diffusivity(1,cell)
   mobility_ref = mobility_ref + grid_orientation_calculate_mobility(1,cell)
 enddo; enddo; enddo
 
end subroutine grid_thermal_spectral_init
!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral crystallinity residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_orientation_formResidual(in,x_scal,f_scal,ierr)

 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
   in
 PetscScalar, dimension( &
   XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, dimension( &
   X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscErrorCode :: ierr
 integer :: i, j, k, cell

  PetscErrorCode :: ierr
  integer :: i, j, k, cell
!  real(pReal)   :: Tdot, dTdot_dT

! evaluate polarization field
  scalarField_real = 0.0_pReal
  scalarField_real(1:grid(1),1:grid(2),1:grid3) = orientation_current 
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of  orientation field
  call utilities_FFTvectorBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    ! vectorField_real(1:3,i,j,k) = matmul(grid_orientation_calculate_diffusivity(1,cell) - D_ref, &
                                               ! vectorField_real(1:3,i,j,k))
    vectorField_real(1:3,i,j,k) = (grid_orientation_calculate_diffusivity(1,cell) - D_ref)* &
                                               vectorField_real(1:3,i,j,k)     !orientation conduction in 
                                                                               !<'2 h(eta) + \frac{g(eta)}{|\nabla \theta|}' in eqn 40
                                                                               !< D_ref = 2 h(eta) + \frac{g(eta)}{|\nabla \theta|} initial time step
                                                                               !< assuming that diffusivity is a scalar
  enddo; enddo; enddo
  call utilities_FFTvectorForward
  call utilities_fourierVectorDivergence                                                            !< calculate orientation divergence in fourier field
!assuming that orientation is a scalar quantity
  call utilities_FFTscalarBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
!    call thermal_conduction_getSourceAndItsTangent(Tdot, dTdot_dT, temperature_current(i,j,k), 1, cell) 
! no Tdot like terms in orientation field solvers
    scalarField_real(i,j,k) = params%timeinc*scalarField_real(i,j,k) + &
                              mobility(1,cell)*(temperature_lastInc(i,j,k)  - &
                                                                          temperature_current(i,j,k)) + &
                              mobility_ref*temperature_current(i,j,k)
  enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of damage field with green operator
  call utilities_FFTscalarForward
  call utilities_fourierGreenConvolution(D_ref, mobility_ref, params%timeinc)
  call utilities_FFTscalarBackward
 
!--------------------------------------------------------------------------------------------------
! constructing residual
  f_scal = orientation_current - scalarField_real(1:grid(1),1:grid(2),1:grid3)

end subroutine spectral_orientation_formResidual

!--------------------------------------------------------------------------------------------------
!> @brief calculates diffusivity
!--------------------------------------------------------------------------------------------------
function grid_orientation_calculate_diffusivity(ip,el)
!< currently the structure is quite similar to thermal_conduction_getConductivity33
 integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el 

 real :: grid_orientation_calculate_diffusivity
 
 integer :: i, j, k, cell
 write(6,*) 'calculating orientation diffusivity'; flush(6)

 cell=0
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
   crystallinity_c(i,j,k) = crystallinity(material_homogenizationAt(cell))% &
                            p(crystallinityMapping(material_homogenizationAt(cell))%p(1,cell))
 enddo; enddo;enddo

  scalarField_real = 0.0_pReal
  scalarField_real(1:grid(1),1:grid(2),1:grid3) = orientation_current 
  call utilities_FFTscalarForward
  call utilities_fourierScalarGradient                                                              !< calculate gradient of  orientation field
  call utilities_FFTvectorBackward
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    oriNorm(i,j,k) = norm2(vectorField_real(1:3,i,j,k))
  enddo; enddo;enddo
  
 grid_orientation_calculate_diffusivity(ip,el) = 2*(d*crystallinity_c**2 + e) + &
               (a1*crystallinity_c + a2*crystallinity_c**2 + a3*crystallinity_c**3) &
                *math_normByItself(oriNorm)
                !< assuming the diffusivity to be a scalar currently
                !< diffusivity = 2h(eta) + g(eta)/(|\nabla theta|)

end function grid_orientation_calculate_diffusivity
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates mobility
!--------------------------------------------------------------------------------------------------
function grid_orientation_calculate_mobility(ip,el)

 integer :: i, j, k, cell
 write(6,*) 'calculating orientation mobility'; flush(6)
 
 cell=0
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
   temp(i,j,k)           =  temperature(material_homogenizationAt(cell))% &
                            p(thermalMapping(material_homogenizationAt(cell))%p(1,cell))
   crystallinity_c(i,j,k) = crystallinity(material_homogenizationAt(cell))% &
                            p(crystallinityMapping(material_homogenizationAt(cell))%p(1,cell))
 enddo; enddo;enddo
 scalarField_real = 0.0_pReal
 scalarField_real(1:grid(1),1:grid(2),1:grid3) = crystallinity_c 
 call utilities_FFTscalarForward()
 call utilities_fourierScalarGradient()                                                             !< calculate gradient of damage field
 call utilities_FFTvectorBackward()
 grad_crystallinity(1:3,1:grid(1),1:grid(2),1:grid3) = vectorField_real(1:3,1:grid(1),1:grid(2),1:grid3)
                                  
 grid_orientation_calculate_mobility(ip,el) =  beta_ori0 * exp(Q/(kB*temp)) * omega**2 *crystallinity_c**2 &  ! philmag 3638 for arrhenius term
             *(1.0_pReal + beta_r*(1.0_pReal- exp(-beta_s*crystallinity_c*norm2(math_normByItself(grad_crystallinity))))) !< ToDo: correct order of norm2/math_normByItself?


end function grid_orientation_calculate_mobility