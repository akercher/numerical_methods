!!$ Program: csi721_project2.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 04/01/2011
!!$ Description: Main program for CSI 721 project 2. 

PROGRAM csi721_project2
USE facedef
USE handleData
USE createStruc
USE elementOps
USE solveOps

  IMPLICIT NONE 

!!$  INTEGER,PARAMETER :: nnode = 3   !number of nodes per element
  INTEGER,PARAMETER :: test = 2,mesh_type = 3,kf = 10000
  INTEGER :: i,j,k                  !Used for Indexing
  INTEGER :: cnt
  INTEGER :: nval
  INTEGER :: ndim               !number of spatial dimentions
  INTEGER :: fdim               !number of dimentions for field
!!$  INTEGER :: npoin              !number of points
!!$  INTEGER :: nnode              !nodes per element
!!$  INTEGER :: nelem              !number of elements
!!$  INTEGER :: nboun              !number of boundary points
  INTEGER :: nface              !number of faces
!!$  INTEGER :: nunks              !number of unknowns
  INTEGER :: noutput            !corresponds to project number
  INTEGER :: bufferI          !buffer for integer values
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
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: elvol   !Element Area Array
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: rho  !density.
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Pr  !pressure
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: en  !energy
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: xyz  !point coordinates unstructured
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: xyz3D  !point coordinates unstructured
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: vfld
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: vxyz !velocity field
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
  ReaD(1001,*) bufferI
  ReaD(1001,*) bufferC
  ReaD(1001,*) bufferC
  ReaD(1001,*) bufferC
  ReaD(1001,*) bufferC
  ReaD(1001,*) ndim,nnode
  ReaD(1001,*) bufferC
  ReaD(1001,*) nelem,npoin,nboun

  CLOSE(1001)

  !define number of unknowns
  nunks = ndim + 2

  !allocations relating to grid
  ALLOCATE(xyz(ndim,npoin),vxyz(ndim,npoin),lnode(nnode,nelem),bpts(2,nboun),initVal(4,npoin))
  ALLOCATE(Nxyz(ndim*nnode,nelem))
  ALLOCATE(uxyz(nunks,npoin),ulast(nunks,npoin),uele(nunks,nelem),duxyz(nunks,npoin))
  ALLOCATE(dt(nelem),dtp(npoin),elvol(nelem),rho(npoin),Pr(npoin),en(npoin))

  !read in rest of data
  CALL dataInput(ndim,lnode,bpts,xyz,initVal,fid_in)

  dt(:) = 0.0d0
  dtp(:) = 0.0d0
  elvol(:) = 0.0d0
  Pr(:) = 0.0d0
  Nxyz(:,:) = 0.0d0

  v2 = vxinf*vxinf + vyinf*vyinf
  eninf = (Prinf/((gamma - 1.0d0)*rhoinf)) + v2/2.0d0

  DO i=1,npoin
     rho(i) = rhoinf
     vxyz(1,i) = vxinf
     vxyz(2,i) = vyinf
     en(i) = eninf
  ENDDO
  DEALLOCATE(initVal)

  nface = nboun
  ALLOCATE(bface(nface))
  ALLOCATE(bnode(nboun))

  !create boundary face data structure
  CALL createFaces(ndim,nface,lnode,bpts,xyz,bface)

!!$  ENDIF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i=1,nboun
     bnode(i)%bpt = bpts(1,i)
     bnode(i)%btype = bpts(2,i)
  ENDDO

  !set initial unknowns
  uxyz(:,:) = 0.0d0             !unknowns at nodes
  uele(:,:) = 0.0d0             !unknowns for element
  ulast(:,:) = 0.0d0            !previous value of unknowns at nodes
  duxyz(:,:) = 0.0d0            !change in unknowns at nodes
  DO i=1,npoin
     uxyz(1,i) = rho(i)
     uxyz(2,i) = rho(i)*vxyz(1,i)
     uxyz(3,i) = rho(i)*vxyz(2,i)
     uxyz(4,i) = rho(i)*en(i)
  ENDDO

  !set previous value to current value
  ulast(:,:) = uxyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Execute Calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TA = 0.0d0                    ! Initialize Total Area of Computational Grid

  ! Call function to calculate area of each element
  elvol = ElementArea(lnode,xyz)
  TA = SUM(elvol)

  ! Call function to compute bounding box coordinates
  Bndcoor = BndBox(xyz)

  xmin  = bndCoor(1,1)
  xmax  = bndCoor(1,2)
  ymin  = bndCoor(2,1)
  ymax  = bndCoor(2,2)

  !calculate shape function derivatives for each element
  CALL shapeDerivatives(ndim,lnode,elvol,xyz,Nxyz)

  !calculate face normals
  CALL faceNormal(ndim,nface,lnode,xyz,Nxyz,bface)

  !calculate boundary normals
  CALL bnodeNormal(ndim,nface,bface,bnode)

!!$GOTO 12345

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! converge to steady state !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  uele(:,:) = 0.0d0
  ulast(:,:) = uxyz(:,:)

  dt(:) = 0.0d0
  duxyz(:,:) = 0.0d0   
  DO k=1,kf
     ulast(:,:) = uxyz(:,:)

     !calculate timestep for each element
     CALL timestepFEM(ndim,lnode,elvol,xyz,Nxyz,uxyz,dt)

     !calculate element unknowns at halfstep
     CALL eulerHalfStep(ndim,lnode,dt,Nxyz,uxyz,uele)  

     !solve system of equations
     CALL solveEuler(ndim,nface,lnode,bface,dt,elvol,xyz,Nxyz,uxyz,uele,duxyz)

     uxyz(:,:) = ulast(:,:) + duxyz(:,:)

     !apply boundary conditions
     CALL eulerBCs(ndim,lnode,bnode,xyz,ulast,uxyz)

     !sum the squares of the difference in unkowns
     CALL sumSq(k,duxyz)
     
  ENDDO

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Project 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write results to file in zfem format
  noutput = 2
  DO i=1,npoin
     rho(i) = uxyz(1,i)
     vxyz(1,i) = uxyz(2,i)/rho(i)
     vxyz(2,i) = uxyz(3,i)/rho(i)
     en(i) = uxyz(4,i)/rho(i)
  ENDDO

  !field dimension
  fdim = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(xyz3D(fdim,npoin),vfld(fdim,npoin))
  xyz3D(:,:) = 0.0d0
  vfld(:,:) = 0.0d0

  xyz3D(1,:) = xyz(1,:)
  xyz3D(2,:) = xyz(2,:)

  vfld(1,:) = vxyz(1,:)
  vfld(2,:) = vxyz(2,:)

  IF(Minf .LT. 0.50d0)THEN
     WRITE(eulerOut,"(A12,I1,A5)"),'data/euler_0',4,'a.vtk'
!!$     WRITE(eulerOut,"(A12,I1,A5)"),'data/euler_0',4,'b.vtk'
  ELSEIF(Minf .GT. 0.50d0)THEN
     WRITE(eulerOut,"(A12,I1,A5)"),'data/euler_0',8,'a.vtk'
!!$     WRITE(eulerOut,"(A12,I1,A5)"),'data/euler_0',8,'b.vtk'
  ENDIF

  Pr(:) = 0.0d0
  !calculate pressure
  nval = npoin
  CALL pressureCalc(ndim,nval,uxyz,Pr)    

  CALL vtkIO(lnode,xyz3D,Pr,vfld,eulerOut)

  WRITE(velOut,"(A18)"),'data/velocity.zfem'
  WRITE(densOut,"(A17)"),'data/density.zfem'
  WRITE(presOut,"(A18)"),'data/pressure.zfem'
  WRITE(enerOut,"(A16)"),'data/energy.zfem'

  CALL zfemIO(fdim,lnode,xyz3D,vfld,velOut)
  !field dimension
  fdim = 1
  CALL zfemIO(fdim,lnode,xyz3D,rho,densOut)
  CALL zfemIO(fdim,lnode,xyz3D,Pr,presOut)
  CALL zfemIO(fdim,lnode,xyz3D,en,enerOut)

12345 CONTINUE

END PROGRAM csi721_project2
