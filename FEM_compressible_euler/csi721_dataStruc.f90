!!$ Program: csi721_dataStruc.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 05/01/2011
!!$ Description: Module to create grid data structures

MODULE createStruc
USE facedef

  INTERFACE createFaces
     MODULE PROCEDURE createFaces_tri
  END INTERFACE

  INTERFACE faceNormal
     MODULE PROCEDURE faceNormal_tri
  END INTERFACE

  INTERFACE bnodeNormal
     MODULE PROCEDURE bnodeNormal_tri
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FACES (TRIANGELE) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE createFaces_tri(ndim,nface,lnode,bpts,xyz,bface)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nface
    INTEGER :: i,j,k
    INTEGER :: cnt
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    INTEGER,INTENT(IN),DIMENSION(2,nboun) :: bpts
    INTEGER,ALLOCATABLE,DIMENSION(:) :: cnt_array
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz
    TYPE(FACE),INTENT(OUT),DIMENSION(nface) :: bface
    LOGICAL :: fndELE

    !initialize face data
    DO i=1,nface
       bface(i)%eleid = -1
       bface(i)%fnode(:) = 0
    ENDDO

    ALLOCATE(cnt_array(nboun))
    cnt_array(:) = 0

    cnt = 1

    !loop over boundary points
    DO i=1,nboun

       DO WHILE(cnt_array(i) .LT. 2)

          fndELE = .FALSE.
          j = 0

          DO WHILE((fndELE .EQV. .FALSE.) .AND. (j .LT. nelem))
             j = j + 1
             !check to see if element face already accounted for
             DO k=1,nface
                IF(bface(k)%eleid .EQ. j)THEN
                   GOTO 1000
                ENDIF
             ENDDO

             !check 1st component of connectivity
             IF(bpts(1,i) .EQ. lnode(1,j))THEN
                bface(cnt)%eleid = j
                bface(cnt)%btype = bpts(2,i)
                bface(cnt)%fnode(1) = bpts(1,i)
                
                cnt_array(i) = cnt_array(i) + 1

                !find other boundary point
                DO k=1,i-1
                   IF(bpts(1,k) .EQ. lnode(2,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(3,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO
                DO k=i+1,nboun
                   IF(bpts(1,k) .EQ. lnode(2,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(3,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO
                
             !check 2nd component of connectivity
             ELSEIF(bpts(1,i) .EQ. lnode(2,j))THEN
                bface(cnt)%eleid = j
                bface(cnt)%btype = bpts(2,i)
                bface(cnt)%fnode(1) = bpts(1,i)
                
                cnt_array(i) = cnt_array(i) + 1
                
                !find other boundary point
                DO k=1,i-1
                   IF(bpts(1,k) .EQ. lnode(1,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(3,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO
                DO k=i+1,nboun
                   IF(bpts(1,k) .EQ. lnode(1,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(3,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO
                
             !check 3rd component of connectivity             
             ELSEIF(bpts(1,i) .EQ. lnode(3,j))THEN
                bface(cnt)%eleid = j
                bface(cnt)%btype = bpts(2,i)
                bface(cnt)%fnode(1) = bpts(1,i)
                
                cnt_array(i) = cnt_array(i) + 1

                !find other boundary point
                DO k=1,i-1
                   IF(bpts(1,k) .EQ. lnode(1,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(2,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO
                DO k=i+1,nboun
                   IF(bpts(1,k) .EQ. lnode(1,j))THEN
                      bface(cnt)%fnode(2) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ELSEIF(bpts(1,k) .EQ. lnode(2,j))THEN
                      bface(cnt)%fnode(2) = bface(cnt)%fnode(1)
                      bface(cnt)%fnode(1) = bpts(1,k)
                      fndELE = .TRUE.
                      cnt_array(k) = cnt_array(k) + 1
                      cnt = cnt + 1
                   ENDIF
                ENDDO

             ENDIF

1000         CONTINUE

          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(cnt_array)

  END SUBROUTINE createFaces_tri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FACE NORMALS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE faceNormal_tri(ndim,nface,lnode,xyz,Nxyz,bface)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nface
    INTEGER :: i,j,k
    INTEGER :: fn1,fn2
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy,Nmag
    DOUBLE PRECISION :: ndot,xn,yn,xt,yt    !dot product of normal
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba
    DOUBLE PRECISION,DIMENSION(ndim) :: fhat !normal for proscribed for (4 degrees) 
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    TYPE(FACE),INTENT(INOUT),DIMENSION(nface) :: bface

    fhat(:) = 0.0d0
    fhat(1) = DCOS(alpha)
    fhat(2) = DSIN(alpha)

    DO i=1,nface
       !get element
       k = bface(i)%eleid
       
       fn1 = bface(i)%fnode(1)
       fn2 = bface(i)%fnode(2)

       !compute unit tangent vector
       a(:) = xyz(:,fn1)
       b(:) = xyz(:,fn2)

       Xba = b(:) - a(:)
       Xba(:) = Xba(:)/DSQRT(Xba(1)*Xba(1) + Xba(2)*Xba(2))

       bface(i)%that(:) = Xba(:)

       !shape derivatives
       Nax = Nxyz(1,k)          !node a of element i, x-direction
       Nay = Nxyz(2,k)          !node a of element i, y-direction

       Nbx = Nxyz(3,k)          !node b of element i, x-direction
       Nby = Nxyz(4,k)          !node b of element i, y-direction

       Ncx = Nxyz(5,k)          !node c of element i, x-direction
       Ncy = Nxyz(6,k)          !node c of element i, y-direction 
       
       Nmag = 0.0d0

       IF((lnode(1,k) .NE. fn1) .AND. (lnode(1,k).NE. fn2))THEN
          Nmag = DSQRT(Nax*Nax + Nay*Nay)
          bface(i)%nhat(1) = -Nax/Nmag
          bface(i)%nhat(2) = -Nay/Nmag
       ELSEIF((lnode(2,k) .NE. fn1) .AND. (lnode(2,k).NE. fn2))THEN
          Nmag = DSQRT(Nbx*Nbx + Nby*Nby)
          bface(i)%nhat(1) = -Nbx/Nmag
          bface(i)%nhat(2) = -Nby/Nmag
       ELSEIF((lnode(3,k) .NE. fn1) .AND. (lnode(3,k).NE. fn2))THEN
          Nmag = DSQRT(Ncx*Ncx + Ncy*Ncy)
          bface(i)%nhat(1) = -Ncx/Nmag
          bface(i)%nhat(2) = -Ncy/Nmag
       ENDIF

       !check boundary: inflow or outflow
       ndot = fhat(1)*bface(i)%nhat(1) + fhat(2)*bface(i)%nhat(2)
       IF(ndot .LE. 0.0d0)THEN
          bface(i)%flow = 0
       ELSE
          bface(i)%flow = 1
       ENDIF

       bface(i)%nhat(1) = bface(i)%nhat(1)
       bface(i)%nhat(2) = bface(i)%nhat(2)
       bface(i)%that(1) = -bface(i)%nhat(2)
       bface(i)%that(2) = bface(i)%nhat(1)

       !compute (INWARD) normal and tangential components of velocity at infinity
       bface(i)%vninf = -(vxinf*bface(i)%nhat(1) + vyinf*bface(i)%nhat(2))
       bface(i)%vtinf = -(vxinf*bface(i)%that(1) + vyinf*bface(i)%that(2))

    ENDDO

  END SUBROUTINE faceNormal_tri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BOUNDARY NODE NORMALS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bnodeNormal_tri(ndim,nface,bface,bnode)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nface
    INTEGER :: i,j,k,cnt
    INTEGER :: bn,fn1,fn2,n1,n2
    DOUBLE PRECISION :: ndot,xn,yn,xt,yt    !dot product of normal
    DOUBLE PRECISION :: nmag
    DOUBLE PRECISION,DIMENSION(ndim) :: fhat !normal for prescribed for (4 degrees) 
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: bface
    TYPE(NODE),INTENT(INOUT),DIMENSION(nboun) :: bnode

    fhat(:) = 0.0d0
    fhat(1) = DCOS(alpha)
    fhat(2) = DSIN(alpha)

    nmag = 0.0d0
    DO i=1,nboun
       bn = bnode(i)%bpt
       cnt = 0
       j = 1
       DO WHILE(cnt .LT. 2)
          fn1 = bface(j)%fnode(1)
          fn2 = bface(j)%fnode(2)
          IF((bn .EQ. fn1) .OR. (bn .EQ. fn2))THEN
             cnt = cnt + 1
             bnode(i)%fid(cnt) = j
          ENDIF
          j = j+1
       ENDDO
       n1 = bnode(i)%fid(1)
       n2 = bnode(i)%fid(2)       

       !find tip of wing
       IF((bface(n1)%btype .EQ. 0) .AND. (bface(n2)%btype .EQ. 0))THEN
          bnode(i)%flow = -1
          ndot = bface(n1)%nhat(1)*bface(n2)%nhat(1) + bface(n1)%nhat(2)*bface(n2)%nhat(2)
          IF(ndot .LT. 0.0d0)THEN
             bnode(i)%flow = -2
          ENDIF
       ENDIF

       bnode(i)%nhat(1) = (bface(n1)%nhat(1) + bface(n2)%nhat(1))/2.0d0
       bnode(i)%nhat(2) = (bface(n1)%nhat(2) + bface(n2)%nhat(2))/2.0d0
       !make sure vector is unit length
       nmag  = DSQRT(bnode(i)%nhat(1)*bnode(i)%nhat(1) + bnode(i)%nhat(2)*bnode(i)%nhat(2))
       bnode(i)%nhat(1) = bnode(i)%nhat(1)/nmag
       bnode(i)%nhat(2) = bnode(i)%nhat(2)/nmag

       !compute if boundary is in flow or outflow
       IF((bface(n1)%btype .EQ. 4) .AND. (bface(n2)%btype .EQ. 4))THEN
          ndot = fhat(1)*bnode(i)%nhat(1) + fhat(2)*bnode(i)%nhat(2)
          IF(ndot .LE. 0.0d0)THEN
             bnode(i)%flow = 0
          ELSE
             bnode(i)%flow = 1
          ENDIF
       ENDIF

       bnode(i)%nhat(1) = bnode(i)%nhat(1)
       bnode(i)%nhat(2) = bnode(i)%nhat(2)
       bnode(i)%that(1) = -bnode(i)%nhat(2)
       bnode(i)%that(2) = bnode(i)%nhat(1)

       !compute normal and tangential components of velocity at infinity
       bnode(i)%vninf = (vxinf*bnode(i)%nhat(1) + vyinf*bnode(i)%nhat(2))
       bnode(i)%vtinf = (vxinf*bnode(i)%that(1) + vyinf*bnode(i)%that(2))

    ENDDO

  END SUBROUTINE bnodeNormal_tri


END MODULE createStruc
