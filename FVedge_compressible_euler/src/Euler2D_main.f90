PROGRAM Euler2D_main
  USE constants
  USE data_types
  USE mesh_data
  USE dataIO, only : vtkIO_legacy
  USE initCalc
  USE initMesh
  USE mesh_check
  USE test_prob
  USE solver
!!$  USE meshOps
!!$  USE initGlobals
!!$  USE variables
!!$  USE timeCLF
!!$  USE fluxes
!!$  USE artVis
!!$  USE limitFCT
!!$  USE BCs
!!$  USE converge
!!$  USE solvers
!!$  USE derivatives

  IMPLICIT NONE 

  INTEGER,PARAMETER :: mesh_type = 8
  INTEGER :: i,j,k                  !Used for Indexing
  INTEGER :: ndim               !number of spatial dimentions
  CHARACTER (LEN=30):: zfemVfld,zfemSfld              ! Output Files: zfem Format
  CHARACTER (LEN=30):: densOut,presOut,enerOut,velOut ! Output Files: zfem Format
  CHARACTER (LEN=80):: eulerOut ! Output File: vtk Format

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Read Input Data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !define values
  nunks     = 4
  max_steps = 1000
  gamma     = 1.40d0
  pgmin     = 0.0d0
  Cr        = 0.950d0
!!$  Cr        = 0.10d0
  Clap      = 2.0d0


  limiter_type = 'none'
  limiter_type = 'vanalbada'
!!$  flux_type    = 'roe'
!!$  flux_type    = 'rhll'
  flux_type    = 'hllc'

  if(mesh_type == 8)then
     prob_type    = 'shock_diff'
  elseif(mesh_type < 8)then
     prob_type    = 'airfoil'
  elseif(mesh_type == 9)then
     prob_type    = 'kelvin_helmholtz'
  elseif(mesh_type == 10)then
     prob_type    = 'shock_tube_2D'
  else
     write(*,*), 'Error: invalid mesh type'
     stop
  endif

  !set problem specific values
  if(trim(prob_type) == 'shock_diff')then
     tf        = 0.180d0
     Minf      = 0.0d0
     step_type = 'global'
  elseif(trim(prob_type) == 'airfoil')then
     tf        = 1000
     alpha     = 4.0d0*pi/180.0d0 
     Minf      = 0.40d0
     Minf      = 0.80d0
     step_type = 'local'
  elseif(trim(prob_type) == 'kelvin_helmholtz')then
     tf        = 1.0d0
     Minf      = 0.0d0
     step_type = 'global'
  elseif(trim(prob_type) == 'shock_tube_2D')then
     tf        = 0.10d0
     alpha     = 0.0d0 
     Minf      = 0.0d0
     step_type = 'global'
  else
     write(*,*), 'Error: problem type not recognized'
     stop
  endif

  if (mesh_type == 8)then
     WRITE(eulerOut,101)
  elseif(mesh_type < 8)then
     if(Minf < 0.50d0)then
        WRITE(eulerOut,102)
     else
        WRITE(eulerOut,103)
     endif
  elseif(mesh_type == 9)then
     WRITE(eulerOut,104)
  elseif(mesh_type == 10)then
     WRITE(eulerOut,105)
  endif
101 FORMAT('data/diffraction_shock_2D.vtk')
102 FORMAT('data/naca0012_04.vtk')
103 FORMAT('data/naca0012_08.vtk')
104 FORMAT('data/kelvin_helmholtz.vtk')
105 FORMAT('data/shock_tube_2D.vtk')

  !read in grid data
  CALL setMesh(mesh_type,ndim)

  !create mesh
  CALL initStruc(ndim)

  !test mesh quantities
  CALL mesh_tests

  !initialize problem
  if(trim(prob_type) == 'shock_diff')then
     CALL init_shock_diffraction
  elseif(trim(prob_type) == 'airfoil')then
     CALL init_airfoil
  elseif(trim(prob_type) == 'kelvin_helmholtz')then
     CALL init_kelvin_helmholtz
  elseif(trim(prob_type) == 'shock_tube_2D')then
     CALL init_shock_tube_2D
  else
     write(*,*), 'Error: problem type not recognized'
     stop
  endif

  CALL FVedge_solver

  CALL vtkIO_legacy(eulerOut)

END PROGRAM Euler2D_main
