!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for nonlocal orientation
!--------------------------------------------------------------------------------------------------
module grid_orientation_spectral 
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScdmda
  use PETScsnes

  use prec
  use spectral_utilities
  use mesh_grid
  use orientation_nonlocal
  use numerics
  use RXconstants
  use crystallite
  use homogenization
  use material
  use FEsolving
  use math


  implicit none
  private

!--------------------------------------------------------------------------------------------------
! derived types
 type(tSolutionParams), private :: params

!--------------------------------------------------------------------------------------------------
! PETSc data
 SNES, private :: orientation_snes
 Vec,  private :: solution_vec
 PetscInt, private :: xstart, xend, ystart, yend, zstart, zend
 real(pReal), private, dimension(:,:,:,:), allocatable :: &
   orientation_current, &                                                                           !< field of current orientation
   orientation_lastInc, &                                                                           !< field of previous orientation
   orientation_stagInc, &                                                                           !< field of staggered orientation
   grad_crystallinity, &
   diffusivity, &
   oriSource
   
 real(pReal), private, dimension(:,:,:), allocatable :: &
   oriNorm, &
   crystallinity_c, &
   temp, &
   mobility

!--------------------------------------------------------------------------------------------------
! reference diffusion tensor, mobility etc. 
 integer,               private :: totalIter = 0                                         !< total iteration in current increment
 real(pReal), dimension(4),   private :: D_ref
 real(pReal), private                 :: mobility_ref
 
 character (len=*), parameter, public :: &
   spectral_orientation_label = 'spectralorientation'
 real(pReal), dimension(:,:),allocatable :: orientation_array                                      !< dummy array to store the orientations
   

 public :: &
   grid_orientation_spectral_init, &
   grid_orientation_spectral_solution, &
   grid_orientation_spectral_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
!--------------------------------------------------------------------------------------------------
subroutine grid_orientation_spectral_init

 PetscInt, dimension(0:worldsize-1) :: localK
 integer :: i, j, k, cell, &
            e, &                                                                                   !< counter in element loop
            c                                                                                      !< counter in integration point component loop
 DM :: orientation_grid
 PetscErrorCode :: ierr
 PetscObject    :: dummy
 PetscScalar, pointer :: x_scal(:,:,:,:)     


 write(6,'(/,a)') ' <<<+-  spectral_orientation init  -+>>>'

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,'-orientation_snes_type newtonls -orientation_snes_mf &
                               &-orientation_snes_ksp_ew -orientation_ksp_type fgmres',ierr)
 CHKERRQ(ierr)
 call PETScOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
 CHKERRQ(ierr)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
 call SNESCreate(PETSC_COMM_WORLD,orientation_snes,ierr); CHKERRQ(ierr)
 call SNESSetOptionsPrefix(orientation_snes,'orientation_',ierr);CHKERRQ(ierr) 
  localK            = 0
  localK(worldrank) = grid3
  call MPI_Allreduce(MPI_IN_PLACE,localK,worldsize,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call DMDACreate3D(PETSC_COMM_WORLD, &
        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                     !< cut off stencil at boundary
        DMDA_STENCIL_BOX, &                                                                         !< Moore (26) neighborhood around central point
        grid(1),grid(2),grid(3), &                                                                  !< global grid
        1, 1, worldsize, &
        4, 0, &                                                                                     !< #dof (orientation phase field), ghost boundary width (domain overlap)
        [grid(1)],[grid(2)],localK, &                                                                   !< local grid
        orientation_grid,ierr)                                                                           !< handle, error
 CHKERRQ(ierr)
 call SNESSetDM(orientation_snes,orientation_grid,ierr); CHKERRQ(ierr)                                        !< connect snes to da
  call DMsetFromOptions(orientation_grid,ierr); CHKERRQ(ierr)
  call DMsetUp(orientation_grid,ierr); CHKERRQ(ierr) 
 call DMCreateGlobalVector(orientation_grid,solution_vec,ierr); CHKERRQ(ierr)                                !< global solution vector (grid x 1, i.e. every def grad tensor)
 call DMDASNESSetFunctionLocal(orientation_grid,INSERT_VALUES,spectral_orientation_formResidual,dummy,ierr)   !< residual vector of same shape as solution vector
 CHKERRQ(ierr) 
  call SNESSetFromOptions(orientation_snes,ierr); CHKERRQ(ierr)                                          ! pull it all together with additional CLI arguments

!--------------------------------------------------------------------------------------------------
! init fields             
 call DMDAGetCorners(orientation_grid,xstart,ystart,zstart,xend,yend,zend,ierr)
 CHKERRQ(ierr)
 xend = xstart + xend - 1
 yend = ystart + yend - 1
 zend = zstart + zend - 1 
 call VecSet(solution_vec,1.0,ierr); CHKERRQ(ierr)                                      !< VS: dont we need upper and lower bounds as done in grid_crystallinity_spectral?
 allocate(crystallinity_c(grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(mobility(grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(temp(grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(diffusivity(4,grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(grad_crystallinity(3,grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(oriSource(4,grid(1),grid(2),grid3), source=0.0_pReal)
 allocate(oriNorm(grid(1),grid(2),grid3), source=0.0_pReal)

 allocate(orientation_array(size(crystallite_orientation,3),4))
 allocate(orientation_current(4,grid(1),grid(2),grid3), source=1.0_pReal)
 allocate(orientation_lastInc(4,grid(1),grid(2),grid3), source=1.0_pReal)
 allocate(orientation_stagInc(4,grid(1),grid(2),grid3), source=1.0_pReal)

 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
     do c = 1,homogenization_Ngrains(material_homogenizationAt(e))
       orientation_array(e,:) = crystallite_orientation(c,i,e)%asQuaternion()
     enddo
   enddo
 enddo

 orientation_current = reshape(orientation_array,[4,grid(1),grid(2),grid3])
 orientation_lastInc = reshape(orientation_array,[4,grid(1),grid(2),grid3])
 orientation_stagInc = reshape(orientation_array,[4,grid(1),grid(2),grid3])
 
 call DMDAVecGetArrayF90(orientation_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)                              !< get the data out of PETSc to work with
 x_scal(0:3,xstart:xend,ystart:yend,zstart:zend) = orientation_current
 call DMDAVecRestoreArrayF90(orientation_grid,solution_vec,x_scal,ierr); CHKERRQ(ierr)

                                  
 quaternionField_fourier = cmplx(0.0,0.0)
 quaternionField_real(1:4,1:grid(1),1:grid(2),1:grid3) = reshape(orientation_array,[4,grid(1),grid(2),grid3])
 call utilities_FFTquaternionForward()
 call utilities_fourierQuaternionGradient()
 call utilities_FFTquatGradBackward()
 
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   oriNorm(i,j,k)=norm2([norm2(quatGradField_real(1,:,i,j,k)),norm2(quatGradField_real(2,:,i,j,k)),&
            norm2(quatGradField_real(3,:,i,j,k)),norm2(quatGradField_real(4,:,i,j,k))])
 enddo; enddo; enddo
!--------------------------------------------------------------------------------------------------
! orientation reference diffusion update
 call grid_orientation_spectral_calculateMobilityAndDiffusivity()
 D_ref = sum(sum(sum(diffusivity,2),2),2) *wgt
 call MPI_Allreduce(MPI_IN_PLACE,D_ref,4,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 mobility_ref = sum(sum(sum(mobility,1),1),1) *wgt
 call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 write(6,*) 'max|min|avg mob orientation init', maxval(mobility), minval(mobility),mobility_ref

end subroutine grid_orientation_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the spectral orientation scheme with internal iterations
!--------------------------------------------------------------------------------------------------
type(tSolutionState) function grid_orientation_spectral_solution(timeinc,timeinc_old)

!--------------------------------------------------------------------------------------------------
! input data for solution
 real(pReal), intent(in) :: &
   timeinc, &                                                                                       !< increment in time for current solution
   timeinc_old                                                                                      !< increment in time of last increment
 integer :: i, j, k, cell
 PetscReal ::  minorientation, maxorientation, stagNorm, solnNorm
  real(pReal), dimension(3,3) :: maxRate, rot
!--------------------------------------------------------------------------------------------------
! PETSc Data
 PetscErrorCode :: ierr   
 SNESConvergedReason :: reason

 grid_orientation_spectral_solution%converged =.false.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data 
 params%timeinc = timeinc
 params%timeincOld = timeinc_old
 call grid_orientation_spectral_calculateMobilityAndDiffusivity()
 
 call SNESSolve(orientation_snes,PETSC_NULL_VEC,solution_vec,ierr); CHKERRQ(ierr)
 call SNESGetConvergedReason(orientation_snes,reason,ierr); CHKERRQ(ierr)

 if (reason < 1) then
   grid_orientation_spectral_solution%converged = .false.
   grid_orientation_spectral_solution%iterationsNeeded = itmax
 else
   grid_orientation_spectral_solution%converged = .true.
   grid_orientation_spectral_solution%iterationsNeeded = totalIter
 endif
 stagNorm = maxval(abs(orientation_current - orientation_stagInc))
 solnNorm = maxval(abs(orientation_current))
 call MPI_Allreduce(MPI_IN_PLACE,stagNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,solnNorm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 orientation_stagInc = orientation_current
 grid_orientation_spectral_solution%stagConverged =     stagNorm < err_damage_tolAbs &
                                         .or. stagNorm < err_damage_tolRel*solnNorm

 quaternionField_fourier = cmplx(0.0,0.0)
 quaternionField_real(1:4,1:grid(1),1:grid(2),1:grid3) = reshape(orientation_array,[4,grid(1),grid(2),grid3])
 call utilities_FFTquaternionForward()
 call utilities_fourierQuaternionGradient()
 call utilities_FFTquatGradBackward()

 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   oriNorm(i,j,k)=norm2([norm2(quatGradField_real(1,:,i,j,k)),norm2(quatGradField_real(2,:,i,j,k)),&             ! 
                         norm2(quatGradField_real(3,:,i,j,k)),norm2(quatGradField_real(4,:,i,j,k))])
 enddo; enddo; enddo
 !grid_orientation_spectral_solution%converged =.true.
!--------------------------------------------------------------------------------------------------
! updating orientation state 
 cell = 0                                                                                      !< material point = 0
 maxRate = 0.0_pReal
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   cell = cell + 1
!   rot = math_mul33x33(transpose(math_qToR(orientation_lastInc(1:4,i,j,k))),&
!                                 math_qToR(orientation_current(1:4,i,j,k)))/params%timeinc
!   rot(1,1) = 0.0_pReal
!   rot(2,2) = 0.0_pReal
!   rot(3,3) = 0.0_pReal
!   crystallinity_Li(:,:,cell) = rot
!   if (norm2(rot) > norm2(maxRate)) &
!     maxRate = rot
   call orientation_nonlocal_putNonLocalorientation([orientation_current(1:4,i,j,k),oriNorm(i,j,k)],1,cell)
 enddo; enddo; enddo
   write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'max Rate:',&
                transpose(maxRate)
 call VecMin(solution_vec,PETSC_NULL_INTEGER,minorientation,ierr); CHKERRQ(ierr)
 call VecMax(solution_vec,PETSC_NULL_INTEGER,maxorientation,ierr); CHKERRQ(ierr)

 if (worldrank == 0) then 
   if (grid_orientation_spectral_solution%converged) &
     write(6,'(/,a)') ' ... nonlocal orientation converged .....................................'
   write(6,'(/,a,f8.6,2x,f8.6,2x,f8.6,/)',advance='no') ' Minimum|Maximum|Delta orientation      = ',&
                                                       minorientation, maxorientation, stagNorm
   write(6,'(/,a)') ' ==========================================================================='
   flush(6) 
 endif 

end function grid_orientation_spectral_solution


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral orientation residual vector
!--------------------------------------------------------------------------------------------------
subroutine spectral_orientation_formResidual(in,x_scal,f_scal,dummy,ierr)

 DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
   in
 PetscScalar, dimension( &
   4,XG_RANGE,YG_RANGE,ZG_RANGE) :: &
   x_scal
 PetscScalar, dimension( &
   4,X_RANGE,Y_RANGE,Z_RANGE) :: &
   f_scal
 PetscObject :: dummy
 PetscErrorCode :: ierr
 integer :: i, j, k, o, cell
 real(pReal)   :: phiDot, dPhiDot_dPhi, H=1.0_pReal, T = 273.0_pReal
 integer, dimension(3) :: other

 orientation_current = x_scal
 
 write(6,*) 'calculating orientation residual'; flush(6)
!--------------------------------------------------------------------------------------------------
! evaluate polarization field
 quaternionField_real = 0.0_pReal
 quaternionField_real(1:4,1:grid(1),1:grid(2),1:grid3) = orientation_current 
 call utilities_FFTquaternionForward()
 call utilities_fourierQuaternionGradient()                                                             !< calculate gradient of orientation field
 call utilities_FFTquatGradBackward()
 
 ! 10.1209/epl/i2005-10081-7 eq(2)
 do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
   oriNorm(i,j,k)=norm2([norm2(quatGradField_real(1,:,i,j,k)),norm2(quatGradField_real(2,:,i,j,k)),&
                         norm2(quatGradField_real(3,:,i,j,k)),norm2(quatGradField_real(4,:,i,j,k))])
 enddo; enddo; enddo
 ! D_ref = sum(sum(sum(diffusivity,2),2),2) *wgt
 ! call MPI_Allreduce(MPI_IN_PLACE,D_ref,4,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 ! mobility_ref = sum(sum(sum(mobility,1),1),1) *wgt
 ! call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 

 ! do o = 1, 4
   ! scalarField_real = 0.0_pReal
   ! do i =1, 3
     ! vectorField_real(i,1:grid(1),1:grid(2),1:grid3) = &
          ! diffusivity(o,1:grid(1),1:grid(2),1:grid3)*quatGradField_real(o,i,1:grid(1),1:grid(2),1:grid3)
   ! enddo
   ! call utilities_FFTvectorForward()
   ! call utilities_fourierVectorDivergence()
   ! call utilities_FFTscalarBackward()
   ! oriSource(o,1:grid(1),1:grid(2),1:grid3) = &
      ! orientation_current(o,1:grid(1),1:grid(2),1:grid3)*scalarField_real(1:grid(1),1:grid(2),1:grid3)
 ! enddo
 do o = 1, 4
   other = mod([o,o+1,o+2],4) + 1
   vectorField_real = 0.0_pReal
   vectorField_real(1:3,1:grid(1),1:grid(2),1:grid3) = quatGradField_real(o,1:3,1:grid(1),1:grid(2),1:grid3)
   do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
     vectorField_real(1:3,i,j,k) = matmul(math_I3*(diffusivity(o,i,j,k) - D_ref(o)), &
                                                vectorField_real(1:3,i,j,k))
   enddo; enddo; enddo
   call utilities_FFTvectorForward()
   call utilities_fourierVectorDivergence()                                                           !< calculate orientation divergence in fourier field
   call utilities_FFTscalarBackward()

   do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
     scalarField_real(i,j,k) = params%timeinc*scalarField_real(i,j,k) + &
                               ! params%timeinc*orientation_current(o,i,j,k)* &
                                ! (  oriSource(other(1),i,j,k) &
                                 ! + oriSource(other(2),i,j,k) &
                                 ! + oriSource(other(3),i,j,k))+ &
                               mobility(i,j,k)*orientation_lastInc(o,i,j,k) - &
                               mobility(i,j,k)*orientation_current(o,i,j,k) + &
                               mobility_ref*orientation_current(o,i,j,k)
   enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! convolution of orientation field with green operator
   call utilities_FFTscalarForward()
   call utilities_fourierGreenConvolution(D_ref(o)*math_I3, mobility_ref, params%timeinc)
   call utilities_FFTscalarBackward()
   quaternionField_real(o,1:grid(1),1:grid(2),1:grid3) = scalarField_real(1:grid(1),1:grid(2),1:grid3)
 enddo
!--------------------------------------------------------------------------------------------------
! constructing residual
 f_scal = quaternionField_real(1:4,1:grid(1),1:grid(2),1:grid3) - orientation_current


end subroutine spectral_orientation_formResidual

!--------------------------------------------------------------------------------------------------
!> @brief spectral orientation forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_orientation_spectral_forward(cutBack)

 logical,     intent(in)                     :: cutBack
 PetscErrorCode                              :: ierr
 integer                               :: i, j, k, cell
 DM :: dm_local
 PetscScalar,  dimension(:,:,:,:), pointer     :: x_scal

 if (cutBack) then 
   orientation_current = orientation_lastInc
   orientation_stagInc = orientation_lastInc
!--------------------------------------------------------------------------------------------------
! reverting orientation field state 
   cell = 0
   call SNESGetDM(orientation_snes,dm_local,ierr); CHKERRQ(ierr)
   call DMDAVecGetArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)                            !< get the data out of PETSc to work with
   x_scal(0:3,xstart:xend,ystart:yend,zstart:zend) = orientation_current
   call DMDAVecRestoreArrayF90(dm_local,solution_vec,x_scal,ierr); CHKERRQ(ierr)
   do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
     cell = cell + 1                                                                           
     call orientation_nonlocal_putNonLocalorientation([orientation_current(:,i,j,k),oriNorm(i,j,k)],1,cell)  !< VS: works only for one IP
   enddo; enddo; enddo
 else
!--------------------------------------------------------------------------------------------------
! update rate and forward last inc
 orientation_lastInc = orientation_current

 endif  

!--------------------------------------------------------------------------------------------------
! orientation reference diffusion update
 call grid_orientation_spectral_calculateMobilityAndDiffusivity()
 D_ref = sum(sum(sum(diffusivity,2),2),2) *wgt
 call MPI_Allreduce(MPI_IN_PLACE,D_ref,4,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 mobility_ref = sum(sum(sum(mobility,1),1),1) *wgt
 call MPI_Allreduce(MPI_IN_PLACE,mobility_ref,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 write(6,*) 'max|min|avg mob orientation forward', maxval(mobility), minval(mobility),mobility_ref

 end subroutine grid_orientation_spectral_forward
 
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates mobility
!--------------------------------------------------------------------------------------------------
subroutine grid_orientation_spectral_calculateMobilityAndDiffusivity()

 integer :: proc
 integer :: i, j, k, cell
 write(6,*) 'calculating orientation mobility and updating crystallinity'; flush(6)

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
                                  
 mobility =  beta_ori0 * exp(Q/(kB*temp)) * omega**2 *crystallinity_c**2 &  ! philmag 3638 for arrhenius term
             *(1.0_pReal + beta_r*(1.0_pReal- exp(-beta_s*crystallinity_c*norm2(math_normByItself(grad_crystallinity))))) !< ToDo: correct order of norm2/math_normByItself?
 
 do i=1, 4
   diffusivity(i,:,:,:) = (&
                            2* (d*crystallinity_c**2 + e)         &                                 ! 2 h(eta)
                          + (a1*crystallinity_c + a2*crystallinity_c**2 + a3*crystallinity_c**3) &
                          *math_normByItself(oriNorm) &
                          )
 enddo
end subroutine grid_orientation_spectral_calculateMobilityAndDiffusivity

end module grid_orientation_spectral
