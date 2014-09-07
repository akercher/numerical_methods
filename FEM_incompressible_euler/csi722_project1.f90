!!$ Program: csi722_project1.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 12/15/2011
!!$ Description: Main program for CSI 722 project 1. 

PROGRAM csi722_project1
USE facedef
USE handleData
USE createStruc
USE elementOps
USE solveOps

  IMPLICIT NONE 

  INTEGER,PARAMETER :: test = 2,mesh_type = 3,kf = 50000
  INTEGER :: i,j,k                  !Used for Indexing
  INTEGER :: cnt
  INTEGER :: ndim               !number of spatial dimentions
  INTEGER :: fdim               !number of dimentions for field
  INTEGER :: npoin              !number of points
  INTEGER :: nnode              !nodes per element
  INTEGER :: nelem              !number of elements
  INTEGER :: nboun              !number of boundary points
  INTEGER :: nface              !number of faces
  INTEGER :: nunks              !number of unknowns
  INTEGER :: noutput            !corresponds to project number
  INTEGER :: bufferI          !buffer for integer values
  DOUBLE PRECISION,PARAMETER :: ca = 2.0d0          !Artificial speed of sound
  DOUBLE PRECISION :: xmin,xmax,ymin,ymax,zmin,zmax !Min/Max Coordinates
  DOUBLE PRECISION :: TA          !Total Area of Computational Grid
  DOUBLE PRECISION :: vx,vy,v2,eninf
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: lnode       !Conectivity Array
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: bpts !Boundary Node Array, point and type(0:wall,4:farfield)
  INTEGER,ALLOCATABLE,DIMENSION(:) :: belem !Boundary Element Array
  INTEGER,ALLOCATABLE,DIMENSION(:) :: cnt_array
  DOUBLE PRECISION,DIMENSION(2,2) :: bndCoor        !Bounding Box Coordinates
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: dt   !Element timestep array
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: dtp   !node timestep array
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: EA   !Element Area Array
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Pr  !pressure
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: xyz  !point coordinates unstructured
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: xyz3D  !point coordinates unstructured
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: vfld
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: vxyz !velocity field
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Fxyz !fluxes at points
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Nxyz !shape function derivatives
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: uxyz,ulast,duxyz,uele  !unknowns
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: initVal !initial values
  TYPE(FACE),ALLOCATABLE,DIMENSION(:) :: bface !Boundary Face
  TYPE(NODE),ALLOCATABLE,DIMENSION(:) :: bnode !Boundary Node 
  CHARACTER (LEN=50):: fid_in                 ! Input File
  CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
  CHARACTER (LEN=30):: zfemVfld,zfemSfld              ! Output Files: zfem Format
  CHARACTER (LEN=30):: densOut,presOut,enerOut,velOut ! Output Files: zfem Format
  CHARACTER (LEN=30):: eulerOut ! Output File: vtk Format
  CHARACTER (LEN=30):: pts,ele,bpS,bpL                 ! Input Files
  LOGICAL :: WallFLAG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Read Input Data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(mesh_type .EQ. 4)THEN 
     ! Format file names of input data
     WRITE(pts,11)
11 FORMAT('data/pts')
     WRITE(ele,21)
!!$21 FORMAT('data/ele')
21 FORMAT('data/ele.wall')
     WRITE(bpS,31)
!!$31 FORMAT('data/bp.side')
31 FORMAT('data/bp.wall')

     ! Open files containing input data
     OPEN(10,FILE = pts)
     OPEN(20,FILE = ele)
     OPEN(30,FILE = bpS)
     
     ! READ Parameters Related Grid
     ndim = 2
     nnode = 3
     nunks = 3
     
     READ(10,*),npoin
     READ(20,*),nelem
     READ(30,*),nboun
     
     !allocations relating to grid
     ALLOCATE(xyz(ndim,npoin),vxyz(ndim,npoin),lnode(nnode,nelem),bpts(2,nboun),initVal(4,npoin))
     ALLOCATE(Fxyz(ndim*nnode,nelem),Nxyz(ndim*nnode,nelem))
     ALLOCATE(uxyz(nunks,npoin),ulast(nunks,npoin),uele(nunks,nelem),duxyz(nunks,npoin))
     ALLOCATE(dt(nelem),EA(nelem),Pr(npoin))

     lnode(:,:) = 0
     bpts(:,:) = 0
     xyz(:,:) = 0.0d0
     EA(:) = 0.0d0

     DO i=1,npoin
        READ(10,*),xyz(:,i)
     ENDDO

     DO i=1,nelem
        READ(20,*),lnode(:,i)      ! nodes of element i
     ENDDO

     DO i=1,nboun
        READ(30,*),bpts(1,i),bpts(2,i)      ! nodes of bp.side
     ENDDO

     nface = nboun
     ALLOCATE(bface(nface))
     ALLOCATE(bnode(nboun))
     
     bface(1)%eleid = 1
     bface(1)%btype = 4
     bface(1)%fnode(1) = 1
     bface(1)%fnode(2) = 2
     
     bface(2)%eleid = 2
     bface(2)%btype = 4
     bface(2)%fnode(1) = 5
     bface(2)%fnode(2) = 1
     
     bface(3)%eleid = 3
     bface(3)%btype = 4
     bface(3)%fnode(1) = 2
     bface(3)%fnode(2) = 3

     bface(4)%eleid = 5
     bface(4)%btype = 4
     bface(4)%fnode(1) = 3
     bface(4)%fnode(2) = 4

     bface(5)%eleid = 5
     bface(5)%btype = 4
     bface(5)%fnode(1) = 4
     bface(5)%fnode(2) = 8

     bface(6)%eleid = 8
     bface(6)%btype = 4
     bface(6)%fnode(1) = 9
     bface(6)%fnode(2) = 5

     bface(7)%eleid = 10  
     bface(7)%btype = 4
     bface(7)%fnode(1) = 8
     bface(7)%fnode(2) = 12

     bface(8)%eleid = 13
     bface(8)%btype = 4
     bface(8)%fnode(1) = 14
     bface(8)%fnode(2) = 13
     
     bface(9)%eleid = 13  
     bface(9)%btype = 4
     bface(9)%fnode(1) = 13
     bface(9)%fnode(2) = 9
     
     bface(10)%eleid = 15
     bface(10)%btype = 4
     bface(10)%fnode(1) = 15
     bface(10)%fnode(2) = 14
     
     bface(11)%eleid = 16
     bface(11)%btype = 4
     bface(11)%fnode(1) = 12
     bface(11)%fnode(2) = 16

     bface(12)%eleid = 17
     bface(12)%btype = 4
     bface(12)%fnode(1) = 16
     bface(12)%fnode(2) = 15
     
     bface(13)%eleid = 4
     bface(13)%btype = 0
     bface(13)%fnode(1) = 7
     bface(13)%fnode(2) = 6

     bface(14)%eleid = 11
     bface(14)%btype = 0
     bface(14)%fnode(1) = 11
     bface(14)%fnode(2) = 7

     bface(15)%eleid = 9
     bface(15)%btype = 0
     bface(15)%fnode(1) = 6
     bface(15)%fnode(2) = 11

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  ELSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  ! Format file names of input data
  IF(mesh_type .EQ. 1)THEN
     WRITE(fid_in,101)
  ELSEIF(mesh_type .EQ. 2)THEN   
     WRITE(fid_in,102)
  ELSEIF(mesh_type .EQ. 3)THEN   
     WRITE(fid_in,103)
  ENDIF
101 FORMAT('naca0012/naca0012.mesh.coarse')
102 FORMAT('naca0012/naca0012.mesh.medium')
103 FORMAT('naca0012/naca0012.mesh.fine')

  ! Open file containing input data
  OPEN(1001,FILE = fid_in)

  !read in problem defining data:ndim,nnodes...
  READ(1001,*) bufferI
  READ(1001,*) bufferC
  READ(1001,*) bufferC
  READ(1001,*) bufferC
  READ(1001,*) bufferC
  READ(1001,*) ndim,nnode
  READ(1001,*) bufferC
  READ(1001,*) nelem,npoin,nboun

  CLOSE(1001)

  !define number of unknowns
  nunks = 3

  !allocations relating to grid
  ALLOCATE(xyz(ndim,npoin),vxyz(ndim,npoin),lnode(nnode,nelem),bpts(2,nboun),initVal(4,npoin))
  ALLOCATE(Fxyz(ndim*nnode,nelem),Nxyz(ndim*nnode,nelem))
  ALLOCATE(uxyz(nunks,npoin),ulast(nunks,npoin),uele(nunks,nelem),duxyz(nunks,npoin))
  ALLOCATE(dt(nelem),dtp(npoin),EA(nelem),Pr(npoin))

  !read in rest of data
  CALL dataInput(ndim,npoin,nnode,nelem,nboun,lnode,bpts,xyz,initVal,fid_in)

  dt(:) = 0.0d0
  dtp(:) = 0.0d0
  EA(:) = 0.0d0
  Fxyz(:,:) = 0.0d0
  Nxyz(:,:) = 0.0d0

  DEALLOCATE(initVal)

  nface = nboun
  ALLOCATE(bface(nface))
  ALLOCATE(bnode(nboun))

  !create boundary face data structure
  CALL createFaces(ndim,npoin,nnode,nelem,nboun,nface,lnode,bpts,xyz,bface)
ENDIF

  DO i=1,nboun
     bnode(i)%bpt = bpts(1,i)
     bnode(i)%btype = bpts(2,i)
  ENDDO

  !set initial unknowns
  uxyz(:,:) = 0.0d0             !unknowns at nodes
  uele(:,:) = 0.0d0             !unknowns for element
  ulast(:,:) = 0.0d0            !previous value of unknowns at nodes
  duxyz(:,:) = 0.0d0            !change in unknowns at nodes
  dt(:) = 0.0d0

  uxyz(1,:) = Prinf
  uxyz(2,:) = vxinf
  uxyz(3,:) = vyinf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Execute Calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TA = 0.0d0                    ! Initialize Total Area of Computational Grid

  ! Call function to calculate area of each element
  EA = ElementArea(npoin,nnode,nelem,lnode,xyz)
  TA = SUM(EA)

  ! Call function to compute bounding box coordinates
  Bndcoor = BndBox(npoin,xyz)

  xmin  = bndCoor(1,1)
  xmax  = bndCoor(1,2)
  ymin  = bndCoor(2,1)
  ymax  = bndCoor(2,2)

  !calculate shape function derivatives for each element
  CALL shapeDerivatives(ndim,npoin,nelem,nnode,lnode,EA,xyz,Nxyz)

  !calculate face normals
  CALL faceNormal(ndim,npoin,nnode,nelem,nface,lnode,xyz,Nxyz,bface)

  !calculate boundary normals
  CALL bnodeNormal(ndim,nface,nboun,bface,bnode)

!!$GOTO 12345

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! converge to steady state !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=1,kf

     !set previous value to current value
     ulast(:,:) = uxyz(:,:)

     !calculate timestep for each element
     CALL timestepFEM(ndim,npoin,nelem,nnode,nunks,ca,lnode,xyz,Nxyz,uxyz,dt)

     !calculate element unknowns at halfstep
     CALL eulerHalfStep(ndim,npoin,nelem,nnode,nunks,ca,lnode,dt,Nxyz,uxyz,uele)  

     !solve system of equations
     CALL solveEuler(ndim,npoin,nelem,nnode,nface,nunks,ca,lnode,bface,dt,EA,xyz,Nxyz,uxyz,uele,duxyz)

     uxyz(:,:) = ulast(:,:) + duxyz(:,:)

     !apply boundary conditions
     CALL eulerBCs(ndim,npoin,nelem,nnode,nboun,nunks,ca,lnode,bnode,xyz,uxyz)

     !sum the squares of the difference in unkowns
     CALL sumSq(k,nunks,npoin,duxyz)
     
  ENDDO

12345 CONTINUE

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Project 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write results to file in zfem format
  noutput = 2

  !field dimension
  fdim = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(xyz3D(fdim,npoin),vfld(fdim,npoin))
  xyz3D(:,:) = 0.0d0
  vfld(:,:) = 0.0d0

  xyz3D(1,:) = xyz(1,:)
  xyz3D(2,:) = xyz(2,:)

  vfld(1,:) = uxyz(2,:)
  vfld(2,:) = uxyz(3,:)

  WRITE(eulerOut,"(A22)"),'data/incompress_05.vtk'

  CALL vtkIO(npoin,nnode,nelem,lnode,xyz3D,uxyz(1,:),vfld,eulerOut)

  WRITE(velOut,"(A18)"),'data/velocity.zfem'
!!$  WRITE(densOut,"(A17)"),'data/density.zfem'
  WRITE(presOut,"(A18)"),'data/pressure.zfem'
!!$  WRITE(enerOut,"(A16)"),'data/energy.zfem'

  CALL zfemIO(fdim,npoin,nnode,nelem,lnode,xyz3D,vfld,velOut)
  !field dimension
  fdim = 1
  CALL zfemIO(fdim,npoin,nnode,nelem,lnode,xyz3D,uxyz(1,:),presOut)


END PROGRAM csi722_project1
