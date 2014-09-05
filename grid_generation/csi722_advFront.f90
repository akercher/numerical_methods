!!$ Program: csi722_advFront.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 10/22/2011
!!$ Description: Module includes subroutines for advancing front.

MODULE advFrontOps
USE facedef

  INTERFACE unitNormal
     MODULE PROCEDURE unitNormal_2D
  END INTERFACE

  INTERFACE updateFaces
     MODULE PROCEDURE updateFaces_2D
  END INTERFACE

  INTERFACE intersectTest
     MODULE PROCEDURE intersectTest_2D
  END INTERFACE

  INTERFACE getInitGrid
     MODULE PROCEDURE getInitGrid_2D
  END INTERFACE

  INTERFACE printFaces
     MODULE PROCEDURE printFaces_2D
  END INTERFACE

  INTERFACE printElements
     MODULE PROCEDURE printElements_2D
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! unitNormal_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE unitNormal_2D(ndim,dir,Xba,len,nout)
    !Calculate unit normal of face in 2D. 
    !dir: direction of normal, 1:inward; -1:outward
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,dir
    INTEGER :: i,j,k
    DOUBLE PRECISION :: X_inv   !1 over magnitude of Xin multiplied by direction
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim) :: Xba !separation vector
    DOUBLE PRECISION,INTENT(OUT) :: len  !length of separation vector
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(ndim) :: nout !unit normal

    nout(:) = 0.0d0

    len = DSQRT(Xba(1)*Xba(1) + Xba(2)*Xba(2))
    X_inv = DBLE(dir)/len

    nout(1) = -Xba(2)*X_inv
    nout(2) = Xba(1)*X_inv

    DO i=1,ndim
       IF(DABS(nout(i)) .LE. tol)THEN
          nout(i) = 0.0d0
       ENDIF
    ENDDO

  END SUBROUTINE unitNormal_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! updateFaces_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE updateFaces_2D(ndim,nnode,iface,ipoin,nelem,nface,nact,faces,lnode)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nnode,iface,ipoin,nelem
    INTEGER,INTENT(INOUT) :: nface,nact
    INTEGER :: dnface,kface
    INTEGER :: i,j,k,cnt
    INTEGER :: iact,inofa,jface,jnode,poin1,poin2
    TYPE(FACE),INTENT(INOUT),DIMENSION(nface_max) :: faces      !current faces
    INTEGER,INTENT(INOUT),DIMENSION(nnode,nelem_max) :: lnode      !conectivity
    INTEGER,DIMENSION(nnofa,2) :: elfa      !other faces of element
    
    !check if nface is greater then nface_max
    IF(nface .GE. (nface_max - 2))THEN
       WRITE(*,*) 'nface greater than nface_max'
       !CALL subroutine to reduce face array to only active faces
    ENDIF

    !create element
    lnode(1,nelem) = faces(iface)%fnode(1)
    lnode(2,nelem) = faces(iface)%fnode(2)
    lnode(3,nelem) = ipoin

    !make old face inactive and set element number
    faces(iface)%eleid = nelem       
    faces(iface)%state = 0
    
    !update number of active fronts
    nact = nact - 1

    !define other faces of element
    elfa(1,1) = faces(iface)%fnode(1)
    elfa(2,1) = ipoin

    elfa(1,2) = ipoin
    elfa(2,2) = faces(iface)%fnode(2)

    dnface = 0
    DO i=1,2

       !loop over active front
       jface = 0
       cnt = 0

       DO WHILE ((jface .LT. nface) .AND. (cnt .LT. 2))
          jface = jface + 1
          cnt = 0
          poin1 = faces(jface)%fnode(1)
          poin2 = faces(jface)%fnode(2)
          !count number of equal points

          DO inofa=1,nnofa
             IF((elfa(inofa,i) .EQ. poin1) .OR. (elfa(inofa,i) .EQ. poin2))THEN
                cnt = cnt + 1
             ENDIF
          ENDDO
       ENDDO

       IF(cnt .EQ. 2)THEN
          !make face inactive
          faces(jface)%state = 0
          !update number of active fronts
          nact = nact - 1
       ELSE
          !create face
          !update increment for nface
          dnface = dnface + 1
          kface = nface + dnface
          faces(kface)%fid = kface          
          faces(kface)%state = 1
          faces(kface)%eleid = nelem
          faces(kface)%fnode(1) = elfa(1,i)
          faces(kface)%fnode(2) = elfa(2,i)
          !update number of active fronts
          nact = nact + 1
       ENDIF

    ENDDO    

    nface = nface + dnface

  END SUBROUTINE updateFaces_2D
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! intersectTest_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE intersectTest_2D(ndim,npoin,nface,iface,pnew,lpoin,faces,rout)
    !Do two rays intersect? 
    !   0:false, do not intersect 
    !   1:true, intersect
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nface,iface
    INTEGER :: i,j,k,inode,kface
    DOUBLE PRECISION :: n1,n2
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,d !nodes of rays
    DOUBLE PRECISION,DIMENSION(ndim) :: Xba,Xca,Xda,Xdc,Xac,Xbc
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim) :: pnew !new point coordinate
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: lpoin !point coordinates
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: faces !face data
    INTEGER,INTENT(OUT) :: rout

    rout = 0
    b(:) = pnew(:)
    DO inode = 1,nnofa
       !get point coordinates
       a(:) = lpoin(:,faces(iface)%fnode(inode))

       kface = 1

       DO WHILE((rout .EQ. 0) .AND. (kface .LE. nface))

          c(:) = lpoin(:,faces(kface)%fnode(1))
          d(:) = lpoin(:,faces(kface)%fnode(2))
          
          Xba(:) = b(:) - a(:)       !1st ray
          Xca(:) = c(:) - a(:)
          Xda(:) = d(:) - a(:)
          
          Xdc(:) = d(:) - c(:)       !2nd ray
          Xac(:) = a(:) - c(:)
          Xbc(:) = b(:) - c(:)
          
          n1 = 0.0d0
          n2 = 0.0d0
          
          !calculate cross products and dot products
          n1 = (Xba(1)*Xca(2) - Xba(2)*Xca(1))*(Xba(1)*Xda(2) - Xba(2)*Xda(1))
          n2 = (Xdc(1)*Xac(2) - Xdc(2)*Xac(1))*(Xdc(1)*Xbc(2) - Xdc(2)*Xbc(1))
          
          !check to see if they intersect
          IF((n1 .LT. 0.0d0) .AND. (n2 .LT. 0.0d0))THEN
             rout = 1           !they intersect
          ENDIF
          kface = kface + 1
       ENDDO
    ENDDO

  END SUBROUTINE intersectTest_2D   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! getInitGrid_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE getInitGrid_2D(ndim,nnode,dir,npoin,nface,lpoin,faces)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nnode,dir
    INTEGER,INTENT(OUT) :: npoin,nface
    INTEGER :: i,j,k,cnt
    INTEGER :: ipoin,jpoin,kpoin,iface,jface,kface
    INTEGER :: nx,ny
    DOUBLE PRECISION :: dx,dy,xpoin,ypoin,len
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(ndim,npoin_max) :: lpoin !points
    TYPE(FACE),INTENT(OUT),DIMENSION(nface_max) :: faces      !faces
    DOUBLE PRECISION,DIMENSION(ndim) :: Xba !separartion vector

    !initialize
    Xba(:) = 0.0d0
    lpoin(:,:) = 0.0d0
    DO i=1,nface_max
       faces(i)%fid = 0
       faces(i)%eleid = 0
       faces(i)%state = 0
    ENDDO

    !define initial parmeters
    nface = 36
    npoin = 36

    !define points
    dx = 0.50d0
    dy = 0.50d0
    !define larger/outer rectangle
    nx = 8
    ny = 7
    xpoin = 0.0d0
    ypoin = 0.0d0

    kpoin = 0
    DO ipoin = 1,nx-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       xpoin = xpoin + dx       
    ENDDO
    DO ipoin = 1,ny-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       ypoin = ypoin + dy
    ENDDO
    DO ipoin = 1,nx-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       xpoin = xpoin - dx
    ENDDO
    DO ipoin = 1,ny-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       ypoin = ypoin - dy
    ENDDO

    !define faces
    kface = kpoin
    DO iface=1,kface-1
       faces(iface)%fid = iface
       faces(iface)%state = 1
       faces(iface)%fnode(1) = iface       
       faces(iface)%fnode(2) = iface + 1       
    ENDDO

    faces(kface)%fid = kface
    faces(kface)%state = 1
    faces(kface)%fnode(1) = kface       
    faces(kface)%fnode(2) = 1

    !define smaller/inner rectangle
    nx = 4
    ny = 3
    xpoin = 1.0d0
    ypoin = 1.0d0

    DO ipoin = 1,ny-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       ypoin = ypoin + dy
       faces(kpoin)%fid = kpoin
    ENDDO
    DO ipoin = 1,nx-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       xpoin = xpoin + dx
    ENDDO
    DO ipoin = 1,ny-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       ypoin = ypoin - dy
    ENDDO
    DO ipoin = 1,nx-1
       kpoin = kpoin + 1
       lpoin(1,kpoin) = xpoin
       lpoin(2,kpoin) = ypoin
       xpoin = xpoin - dx
    ENDDO

    !define faces
    kface = kface + 1
    DO iface=kface,kpoin-1
       faces(iface)%fid = iface
       faces(iface)%state = 1
       faces(iface)%fnode(1) = iface       
       faces(iface)%fnode(2) = iface + 1       
    ENDDO

    faces(kpoin)%fid = kpoin
    faces(kpoin)%state = 1
    faces(kpoin)%fnode(1) = kpoin       
    faces(kpoin)%fnode(2) = kface
    

!!$       lpoin(1,1) = 0.0d0
!!$       lpoin(2,1) = 0.0d0
!!$    
!!$    lpoin(1,2) = 0.50d0
!!$    lpoin(2,2) = 0.0d0
!!$    
!!$    lpoin(1,3) = 1.0d0
!!$    lpoin(2,3) = 0.0d0
!!$    
!!$    lpoin(1,4) = 1.0d0
!!$    lpoin(2,4) = 0.50d0
!!$    
!!$    lpoin(1,5) = 0.50d0
!!$    lpoin(2,5) = 0.50d0
!!$    
!!$    lpoin(1,6) = 0.0d0
!!$    lpoin(2,6) = 0.50d0

        
    !define initial advancing faces
!!$    faces(1)%fid = 1
!!$    faces(1)%eleid = 0
!!$    faces(1)%state = 1
!!$    faces(1)%fnode(1) = 1
!!$    faces(1)%fnode(2) = 2
!!$    
!!$    faces(2)%fid = 2
!!$    faces(2)%eleid = 0
!!$    faces(2)%state = 1
!!$    faces(2)%fnode(1) = 2
!!$    faces(2)%fnode(2) = 3
!!$    
!!$    faces(3)%fid = 3
!!$    faces(3)%eleid = 0
!!$    faces(3)%state = 1
!!$    faces(3)%fnode(1) = 3
!!$    faces(3)%fnode(2) = 4
!!$    
!!$    faces(4)%fid = 4
!!$    faces(4)%eleid = 0
!!$    faces(4)%state = 1
!!$    faces(4)%fnode(1) = 4
!!$    faces(4)%fnode(2) = 5
!!$    
!!$    faces(5)%fid = 5
!!$    faces(5)%eleid = 0
!!$    faces(5)%state = 1
!!$    faces(5)%fnode(1) = 5
!!$    faces(5)%fnode(2) = 6
!!$    
!!$    faces(6)%fid = 6
!!$    faces(6)%eleid = 0
!!$    faces(6)%state = 1
!!$    faces(6)%fnode(1) = 6
!!$    faces(6)%fnode(2) = 1
    
    !calculate inward normal and face length
    DO iface=1,nface
       Xba(:) = lpoin(:,faces(iface)%fnode(2)) - lpoin(:,faces(iface)%fnode(1))
       CALL unitNormal(ndim,dir,Xba,faces(iface)%flen,faces(iface)%nhat(:))
       !adjust length
       faces(iface)%flen = faces(iface)%flen*(1.0d0 + DBLE(iface)*leps)       
    ENDDO
    
  END SUBROUTINE getInitGrid_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! printFaces_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE printFaces_2D(ndim,npoin,nface,lpoin,faces)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nface
    INTEGER :: i,j,k
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: lpoin
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: faces      !current faces

    WRITE(*,*) 'Faces:'
    DO i=1,nface
       WRITE(*,*) i,'fid:',faces(i)%fid,'state:',faces(i)%state
       WRITE(*,*) ' ',faces(i)%fnode(1),lpoin(1,faces(i)%fnode(1)),lpoin(2,faces(i)%fnode(1))
       WRITE(*,*) ' ',faces(i)%fnode(2),lpoin(1,faces(i)%fnode(2)),lpoin(2,faces(i)%fnode(2))
    ENDDO

  END SUBROUTINE printFaces_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! printElements_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE printElements_2D(ndim,nnode,npoin,nelem,lpoin,lnode)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nnode,npoin,nelem
    INTEGER :: i,j,k
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: lpoin
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    
    WRITE(*,*) 'Elements:'
    DO i=1,nelem
       WRITE(*,*) i
       WRITE(*,*) ' ',lnode(1,i),lpoin(1,lnode(1,i)),lpoin(2,lnode(1,i))
       WRITE(*,*) ' ',lnode(2,i),lpoin(1,lnode(2,i)),lpoin(2,lnode(2,i))
       WRITE(*,*) ' ',lnode(3,i),lpoin(1,lnode(3,i)),lpoin(2,lnode(3,i))
    ENDDO

  END SUBROUTINE printElements_2D

END MODULE advFrontOps
