module grid_orientation_spectral

!public :: &
 
contains

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
  call utilities_fourierScalarGradient                                                              !< calculate gradient of damage field
  call utilities_FFTvectorBackward
  cell = 0
  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
    cell = cell + 1
    vectorField_real(1:3,i,j,k) = matmul(orientation_conduction_getConductivity33(1,cell) - D_ref, &
                                               vectorField_real(1:3,i,j,k))    !orientation conduction in 'h(eta)' in eqn 40
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
                              tau(1,cell)*(temperature_lastInc(i,j,k)  - &
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
 
