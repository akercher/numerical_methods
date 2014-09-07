!!$ Program: Euler2D_initCalc.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module to initialize calculation.

MODULE initCalc

  INTERFACE setMesh
     MODULE PROCEDURE setMesh_tri
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SET MESH (TRIANGELE) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setMesh_tri(mesh_type,ndim)!,nx,ny,nz)
    USE constants
    USE data_types
    USE mesh_data
    USE dataIO

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: mesh_type
    INTEGER,INTENT(INOUT) :: ndim!,nx,ny,nz
    INTEGER :: i,j,k
    INTEGER :: cnt
    INTEGER :: bufferI          !buffer for integer values
    DOUBLE PRECISION :: TA          !Total Area of Computational Grid
    DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: bndCoor !Bounding Box Coordinates
    DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: initVal !initial values
    CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
    CHARACTER (LEN=50):: fid_in  ! Input File
    CHARACTER (LEN=50):: bid_in  ! Boundary condition File
    CHARACTER (LEN=30):: pts,ele,bpS,bpL ! Input Files
    LOGICAL :: fndELE

    ! Format file names of input data
    IF(mesh_type .EQ. 1)THEN
       WRITE(fid_in,101)
       WRITE(bid_in,100)
    ELSEIF(mesh_type .EQ. 2)THEN   
       WRITE(fid_in,102)
       WRITE(bid_in,100)
    ELSEIF(mesh_type .EQ. 3)THEN   
       WRITE(fid_in,103)
       WRITE(bid_in,100)
    ELSEIF(mesh_type .EQ. 4)THEN   
       WRITE(fid_in,104)
    ELSEIF(mesh_type .EQ. 5)THEN   
       WRITE(fid_in,105)
    ELSEIF(mesh_type .EQ. 6)THEN   
       WRITE(fid_in,106)
    ELSEIF(mesh_type .EQ. 7)THEN   
       WRITE(fid_in,107)
    ELSEIF(mesh_type .EQ. 8)THEN   
       WRITE(fid_in,108)
       WRITE(bid_in,1081)
    ELSEIF(mesh_type .EQ. 9)THEN   
       WRITE(fid_in,109)
       WRITE(bid_in,1091)
    ELSEIF(mesh_type .EQ. 10)THEN   
       WRITE(fid_in,1010)
       WRITE(bid_in,10101)
    ENDIF

100 FORMAT('naca0012/naca0012.bcmap')
101 FORMAT('naca0012/naca0012.mesh.coarse')
102 FORMAT('naca0012/naca0012.mesh.medium')
103 FORMAT('naca0012/naca0012.mesh.fine')
104 FORMAT('naca0012/simple.mesh_3x3')
105 FORMAT('naca0012/simple.mesh_3x3_wall')
106 FORMAT('naca0012/simple.mesh_4x4')
107 FORMAT('naca0012/simple.mesh_4x4_wall')
108 FORMAT('project.grid')
1081 FORMAT('project.bcmap')
109 FORMAT('mesh_generator/periodic.grid')
!!$1091 FORMAT('../mesh_generator/project.bcmap')
1091 FORMAT('mesh_generator/periodic.bcmap')
1010 FORMAT('mesh_generator/sod2D.grid')
10101 FORMAT('mesh_generator/sod2D.bcmap')

    CALL dataInput(fid_in,bid_in)

  END SUBROUTINE setMesh_tri



END MODULE initCalc
