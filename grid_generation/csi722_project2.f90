!!$ Program: csi722_project2.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 10/22/2011
!!$ Description: Main program for CSI 722 project 2, advancing front grid generation. 

PROGRAM csi722_project2
USE facedef
USE advFrontOps

  IMPLICIT NONE 

  INTEGER,PARAMETER :: maxsteps = 500 !aviod infinite loop
  INTEGER :: i,j,k,cnt,cntsteps          !Used for Indexing
  INTEGER :: iface,ipoin,iact          !Used for Indexing
  INTEGER :: imin          !Index of face with minimum length
  INTEGER :: nface_new          !number of active fronts
  INTEGER :: nact          !number of active fronts
  INTEGER :: npoin              !number of points
  INTEGER :: nnode              !nodes per element
  INTEGER :: nelem              !number of elements
  INTEGER :: nboun              !number of boundary points
  INTEGER :: nface              !number of faces
  INTEGER :: nfron              !number of active fronts
  INTEGER :: bufferI            !buffer for integer values
  INTEGER :: rout               !do two rays intersect (0:no,1:yes)
  INTEGER :: dir                !direction of normal, 1:inward; -1:outward
  DOUBLE PRECISION :: TA          !Total Area of Computational Grid
  DOUBLE PRECISION :: EA          !Element Area
  DOUBLE PRECISION :: rsep        !magnitude of separation vector
  DOUBLE PRECISION :: rsep_new        !magnitude of separation vector
  DOUBLE PRECISION :: lmin
  DOUBLE PRECISION :: dh
  INTEGER,ALLOCATABLE,DIMENSION(:) :: bpts !Boundary Node Array
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: lnode       !Conectivity Array
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: a,b,c,d !nodes of rays
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: Xba,Xca,Xda,Xdc,Xac,Xbc
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: pmid !midpoint
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: pnew !new point to add
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Xnew !point coordinates of potential new face
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: lpoin  !point coordinates unstructured
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Nlpoin !shape function derivatives
  TYPE(FACE),ALLOCATABLE,DIMENSION(:) :: faces !element faces
  TYPE(FACE),ALLOCATABLE,DIMENSION(:) :: front !advancing front
  CHARACTER (LEN=50):: fid_in                 !input file
  CHARACTER (LEN=50):: face_out                 !faces output file
  CHARACTER (LEN=50):: pts_out                 !points output file
  CHARACTER (LEN=50):: ele_out                 !elements output file
  CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
  LOGICAL :: AddPt,AddFa

  !format output file names
  WRITE(face_out ,"(A10)"),'data/faces'
  WRITE(pts_out ,"(A8)"),'data/pts'
  WRITE(ele_out ,"(A8)"),'data/ele'

  !open files
  OPEN(1001, FILE = face_out)
  OPEN(1002, FILE = pts_out)
  OPEN(1003, FILE = ele_out)

  !define parameters
  nnode = 3
  nboun = 6

  ALLOCATE(bpts(nboun),lpoin(ndim,npoin_max),lnode(nnode,nelem_max),faces(nface_max))
  ALLOCATE(a(ndim),b(ndim),c(ndim),d(ndim))
  ALLOCATE(Xba(ndim),Xca(ndim),Xda(ndim),Xdc(ndim),Xac(ndim),Xbc(ndim),Xnew(ndim,nnofa))
  ALLOCATE(pmid(ndim),pnew(ndim))

  nelem = 0
  lnode(:,:) = 0

  !set direction to inward normal
  dir = 1
  CALL getInitGrid(ndim,nnode,dir,npoin,nface,lpoin,faces)
  nact = nface

  !write current faces to file
  WRITE(1001,*) nface
  DO i=1,nface
     WRITE(1001,*) faces(i)%state,faces(i)%fnode(1) - 1,faces(i)%fnode(2) - 1
  ENDDO

  !uniform grid
  IF(gtype .EQ. 1)THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     WRITE(*,*) 'Generating uniform grid'

     !set height and area of element
     he = 0.250d0
     he2 = he*he

     cntsteps = 0
     DO WHILE ((nact .GT. 0) .AND. (cntsteps .LT. maxsteps))

        !loop over active faces
        nface_new = nface
        DO iface=1,nface_new
           
           !make sure face is still active
           IF(faces(iface)%state .EQ. 1)THEN
              
              !set logical value to add a point and face to TRUE
              AddPt = .TRUE.
              
              !get point coordinates of face
              a(:) = lpoin(:,faces(iface)%fnode(1))
              c(:) = lpoin(:,faces(iface)%fnode(2))
              
              !get midpoint
              b(:) = (c(:) + a(:))/2.0d0
              !define separation vector
              Xba(:) = b(:) - a(:)
              
              !define midpoint
              pmid(:) = a(:) + Xba(:)
              
              !define new point to add
              pnew(:) = pmid(:) + he*faces(iface)%nhat(:)
              
              !check if new point is within certain distance of an existing point
              DO j=1,iface-1
                 DO k=1,nnofa
                    a(:) = lpoin(:,faces(j)%fnode(k))
                    Xba(:) = pnew(:) - a(:)
                    rsep = Xba(1)*Xba(1) + Xba(2)*Xba(2)
                    IF(rsep .LE. he2)THEN
                       AddPt = .FALSE.
                       pnew(:) = a(:)
                       ipoin = faces(j)%fnode(k)
                    ENDIF
                 ENDDO
              ENDDO
              
              DO j=iface+1,nface
                 DO k=1,nnofa
                    a(:) = lpoin(:,faces(j)%fnode(k))
                    Xba(:) = pnew(:) - a(:)
                    rsep = Xba(1)*Xba(1) + Xba(2)*Xba(2)
                    IF(rsep .LE. he2)THEN
                       AddPt = .FALSE.
                       pnew(:) = a(:)
                       ipoin = faces(j)%fnode(k)
                    ENDIF
                 ENDDO
              ENDDO
              
              !check area of new element and make sure it is uniform
              a(:) = lpoin(:,faces(iface)%fnode(1))
              c(:) = lpoin(:,faces(iface)%fnode(2))
              
              Xba(:) = pnew(:) - a(:)
              Xca(:) = c(:) - a(:)
              
              EA = 0.50d0*(Xba(1)*Xca(2) - Xba(2)*Xca(1))
              EA = DABS(EA)
              IF(EA .NE. he2)THEN
                 !if area does not match uniform frid then skip face
                 GOTO 1234
              ENDIF

              !check to see if the new element cross any face
              CALL intersectTest(ndim,npoin,nface,iface,pnew,lpoin,faces,rout)
              
              IF(AddPt .EQV. .TRUE.)THEN
                 npoin = npoin + 1
                 lpoin(:,npoin) = pnew(:)
                 ipoin = npoin
              ENDIF
              
              nelem = nelem + 1
              CALL updateFaces(ndim,nnode,iface,ipoin,nelem,nface,nact,faces,lnode)
              
              !calculate inward normal
              DO j=1,nface
                 Xba(:) = lpoin(:,faces(j)%fnode(2)) - lpoin(:,faces(j)%fnode(1))
                 CALL unitNormal(ndim,dir,Xba,faces(j)%flen,faces(j)%nhat(:))
                 !adjust length
                 faces(j)%flen = faces(j)%flen*(1.0d0 + DBLE(j)*leps)       
              ENDDO
              
!!$              CALL printFaces(ndim,npoin,nface,lpoin(:,1:npoin),faces(1:nface))
!!$              CALL printElements(ndim,nnode,npoin,nelem,lpoin(:,1:npoin),lnode(:,1:nelem))
              
              !write current faces to file
              WRITE(1001,*) nface
              DO i=1,nface
                 WRITE(1001,*) faces(i)%state,faces(i)%fnode(1) - 1,faces(i)%fnode(2) - 1
              ENDDO

              cntsteps = cntsteps + 1

1234          CONTINUE
           ENDIF
        ENDDO
     ENDDO

  !nonunifrom grid
  ELSE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     WRITE(*,*) 'Generating nonuniform grid'

     cntsteps = 0
     iface = 0     
     DO WHILE ((nact .GT. 0) .AND. (cntsteps .LT. maxsteps))

        !find smallest face
        nface_new = nface
        !set lmin and imin to large numbers
        lmin = 1.0E+16
        imin  = nface + 1
        DO iface=1,nface_new
           !make sure face is still active
           IF(faces(iface)%state .EQ. 1)THEN
              IF(faces(iface)%flen .LT. lmin)THEN
                 lmin = faces(iface)%flen
                 imin = faces(iface)%fid
              ENDIF
           ENDIF
        ENDDO
        !set iface equal to imin
        iface = imin

        !set logical value to add a point and face to TRUE
        AddPt = .TRUE.
              
        !get point coordinates of face
        a(:) = lpoin(:,faces(iface)%fnode(1))
        c(:) = lpoin(:,faces(iface)%fnode(2))
        
        !get midpoint
        b(:) = (c(:) + a(:))/2.0d0
        !define separation vector
        Xba(:) = b(:) - a(:)

        rout = 1
        he = 0.260d0
        DO WHILE((rout .EQ. 1) .AND. (he .GT. 0.0d0))

           !define height
!!$           dh = dh + 1.0d0
!!$           he = (tan60*DSQRT(Xba(1)*Xba(1) + Xba(2)*Xba(2)))/dh
           he = he + 0.01;

           he2 = he*he

           !define midpoint
           pmid(:) = a(:) + Xba(:)
           
           !define new point to add
           pnew(:) = pmid(:) + he*faces(iface)%nhat(:)

           !find closest point to pnew
           rsep = 1.0E+16
           DO j=1,iface-1
              DO k=1,nnofa
                 a(:) = lpoin(:,faces(j)%fnode(k))
                 Xca(:) = pnew(:) - a(:)
                 rsep_new = Xca(1)*Xca(1) + Xca(2)*Xca(2)
                 IF(rsep_new .LE. rsep)THEN
                    !check to see if new face intersects any faces
                    CALL intersectTest(ndim,npoin,nface,iface,a,lpoin,faces,rout)
                    IF(rout .EQ. 0)THEN
                       ipoin = faces(j)%fnode(k)
                       rsep = rsep_new
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO

           DO j=iface+1,nface
              DO k=1,nnofa
                 a(:) = lpoin(:,faces(j)%fnode(k))
                 Xca(:) = pnew(:) - a(:)
                 rsep_new = Xca(1)*Xca(1) + Xca(2)*Xca(2)
                 IF(rsep_new .LE. rsep)THEN
                    !check to see if new face intersects any faces
                    CALL intersectTest(ndim,npoin,nface,iface,a,lpoin,faces,rout)
                    IF(rout .EQ. 0)THEN
                       ipoin = faces(j)%fnode(k)
                       rsep = rsep_new
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO

           !check if new point is within certain distance of an existing point
           IF(rsep .LE. he2)THEN
              AddPt = .FALSE.
              pnew(:) = lpoin(:,ipoin)
           ENDIF

           !check to see if the new element cross any face
           CALL intersectTest(ndim,npoin,nface,iface,pnew,lpoin,faces,rout)

        ENDDO

        IF(AddPt .EQV. .TRUE.)THEN
           npoin = npoin + 1
           lpoin(:,npoin) = pnew(:)
           ipoin = npoin
        ENDIF
              
        nelem = nelem + 1
        CALL updateFaces(ndim,nnode,iface,ipoin,nelem,nface,nact,faces,lnode)
        
        !calculate inward normal
        DO j=1,nface
           Xba(:) = lpoin(:,faces(j)%fnode(2)) - lpoin(:,faces(j)%fnode(1))
           CALL unitNormal(ndim,dir,Xba,faces(j)%flen,faces(j)%nhat(:))
           !adjust length
           faces(j)%flen = faces(j)%flen*(1.0d0 + DBLE(j)*leps)       
        ENDDO
              
!!$              CALL printFaces(ndim,npoin,nface,lpoin(:,1:npoin),faces(1:nface))
!!$              CALL printElements(ndim,nnode,npoin,nelem,lpoin(:,1:npoin),lnode(:,1:nelem))
              
        !write current faces to file
        WRITE(1001,*) nface
        DO i=1,nface
           WRITE(1001,*) faces(i)%state,faces(i)%fnode(1) - 1,faces(i)%fnode(2) - 1
        ENDDO
        
        cntsteps = cntsteps + 1
        
     ENDDO
     
  ENDIF

  WRITE(*,*) 'nelem =',nelem
  WRITE(*,*) 'npoin =',npoin

  !write points to file
  WRITE(1002,*) npoin
  DO i=1,npoin
     WRITE(1002,*) lpoin(1,i),lpoin(2,i)
  ENDDO
  
END PROGRAM csi722_project2
