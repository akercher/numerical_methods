!!$ Program: csi721_elementOps.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 04/01/2011
!!$ Description: Module includes subroutines to gather element information.

MODULE elementOps
USE facedef

  INTERFACE BndBox
     MODULE PROCEDURE BndBox_2D, BndBox_3D
  END INTERFACE

  INTERFACE AvgWeightedNorm
     MODULE PROCEDURE AvgWeightedNorm_2D,AvgWeightedNorm_3D
  END INTERFACE   

  INTERFACE ElementArea
     MODULE PROCEDURE ElementArea_2D,ElementArea_3D
  END INTERFACE

  INTERFACE shapeDerivatives
     MODULE PROCEDURE shapeDerivatives_linear
  END INTERFACE   

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BndBox_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION BndBox_2D(xyz) result(Bndcoor)
    IMPLICIT NONE

    DOUBLE PRECISION :: xmax,xmin,ymax,ymin,zmax,zmin
    DOUBLE PRECISION,INTENT(IN),DIMENSION(2,npoin) :: xyz
    DOUBLE PRECISION,DIMENSION(2,2) :: Bndcoor
 
    Bndcoor(1,1) = MINVAL(xyz(1,:))    ! Min x-coor 
    Bndcoor(1,2) = MAXVAL(xyz(1,:))    ! Max x-coor
    Bndcoor(2,1) = MINVAL(xyz(2,:))    ! Min y-coor
    Bndcoor(2,2) = MAXVAL(xyz(2,:))    ! Max y-coor

  END FUNCTION BndBox_2D   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BndBox_3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION BndBox_3D(ndim,xyz) result(Bndcoor)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    DOUBLE PRECISION :: xmax,xmin,ymax,ymin,zmax,zmin
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz
    DOUBLE PRECISION,DIMENSION(ndim,2) :: Bndcoor
 
    Bndcoor(1,1) = MINVAL(xyz(1,:))    ! Min x-coor 
    Bndcoor(1,2) = MAXVAL(xyz(1,:))    ! Max x-coor
    Bndcoor(2,1) = MINVAL(xyz(2,:))    ! Min y-coor
    Bndcoor(2,2) = MAXVAL(xyz(2,:))    ! Max y-coor
    Bndcoor(3,1) = MINVAL(xyz(3,:))    ! Min z-coor
    Bndcoor(3,2) = MAXVAL(xyz(3,:))    ! Max z-coor

  END FUNCTION BndBox_3D   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AvgWeightedNorm_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION AvgWeightedNorm_2D(Xba,Xca) RESULT(AWN)
    IMPLICIT NONE

    INTEGER :: i,j,k
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: Xba,Xca
    DOUBLE PRECISION :: AWN ! Average Weighted Normal
    
    ! Calulate Cross Product Divided by 2
    AWN = 0.50d0*(Xba(1)*Xca(2) - Xba(2)*Xca(1))

  END FUNCTION AvgWeightedNorm_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AvgWeightedNorm_3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION AvgWeightedNorm_3D(ndim,Xba,Xca) RESULT(AWN)
    IMPLICIT NONE

    INTEGER :: i,j,k
    INTEGER,INTENT(IN) :: ndim
    DOUBLE PRECISION,DIMENSION(ndim),INTENT(IN) :: Xba,Xca
    DOUBLE PRECISION,DIMENSION(ndim) :: AWN ! Average Weighted Normal
    
    ! Calulate Cross Product Divided by 2
    AWN(1) = 0.50d0*(Xba(2)*Xca(3) - Xba(3)*Xca(2))
    AWN(2) = 0.50d0*(Xba(1)*Xca(3) - Xba(3)*Xca(1))
    AWN(3) = 0.50d0*(Xba(1)*Xca(2) - Xba(2)*Xca(1))

  END FUNCTION AvgWeightedNorm_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ElementArea_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ElementArea_2D(lnode,xyz) result(elvol)
    IMPLICIT NONE

    INTEGER :: i,j,k
    INTEGER :: n1,n2,n3
    DOUBLE PRECISION :: AWN ! Average Weighted Normal
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(2,npoin) :: xyz
    INTEGER,DIMENSION(3) :: n_array ! Array for temp. storing nodes of element
    DOUBLE PRECISION,DIMENSION(2) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nelem) :: elvol ! Element Area
 
    DO i=1,nelem
       
       n1 = lnode(1,i)
       n2 = lnode(2,i)
       n3 = lnode(3,i)

       ! Nodes of Current Element
       a(:) = xyz(:,n1)
       b(:) = xyz(:,n2)
       c(:) = xyz(:,n3)

       Xba(:) = b(:) - a(:)
       Xca(:) = c(:) - a(:)
     
       ! Calulate Avergae Weight Norm
       AWN = AvgWeightedNorm(Xba,Xca)

       ! Calculate Element Area
       elvol(i) = DABS(AWN)

    ENDDO

  END FUNCTION ElementArea_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ElementArea_3D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION ElementArea_3D(ndim,lnode,xyz) result(elvol)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: i,j,k
    INTEGER :: n1,n2,n3
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz
    INTEGER,DIMENSION(3) :: n_array ! Array for temp. storing nodes of element
    DOUBLE PRECISION,DIMENSION(3) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(3) :: AWN ! Average Weighted Normal
    DOUBLE PRECISION,DIMENSION(nelem) :: elvol ! Element Area
 
    DO i=1,nelem
       
       n1 = lnode(1,i)
       n2 = lnode(2,i)
       n3 = lnode(3,i)

       ! Nodes of Current Element
       a(:) = xyz(:,n1)
       b(:) = xyz(:,n2)
       c(:) = xyz(:,n3)

       Xba(:) = b(:) - a(:)
       Xca(:) = c(:) - a(:)
     
       ! Calulate Avergae Weight Norm
       AWN(:) = AvgWeightedNorm(Xba,Xca)

       ! Calculate Element Area
       elvol(i) = (AWN(1)**2.0d0 + AWN(2)**2.0d0 + AWN(3)**2.0d0)**(1.0d0/2.0d0)

    ENDDO

  END FUNCTION ElementArea_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! shapeDerivatives_linear !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE shapeDerivatives_linear(ndim,lnode,elvol,xyz,Nxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: i,j,k,kstart
    INTEGER :: n1,n2,n3
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz ! Point Coordinates
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(ndim*nnode,nelem) :: Nxyz !element shape derivatives
    DOUBLE PRECISION,DIMENSION(2) :: a,b,c,Xba,Xca

    DO i=1,nelem
       n1 = lnode(1,i)
       n2 = lnode(2,i)
       n3 = lnode(3,i)

       ! Nodes of Current Element
       a(:) = xyz(:,n1)
       b(:) = xyz(:,n2)
       c(:) = xyz(:,n3)

       Xba(:) = b(:) - a(:)
       Xca(:) = c(:) - a(:)
       
       ! Calculate Shape Derivative for Element
       Nxyz(1,i) = (-Xca(2) + Xba(2))/(2.0d0*elvol(i)) !node 1 of element i, x-direction
       Nxyz(2,i) = (Xca(1) - Xba(1))/(2.0d0*elvol(i))  !node 1 of element i, y-direction

       Nxyz(3,i) = Xca(2)/(2.0d0*elvol(i))  !node 2 of element i, x-direction
       Nxyz(4,i) = -Xca(1)/(2.0d0*elvol(i)) !node 2 of element i, y-direction

       Nxyz(5,i) = -Xba(2)/(2.0d0*elvol(i)) !node 3 of element i, x-direction
       Nxyz(6,i) = Xba(1)/(2.0d0*elvol(i))  !node 3 of element i, y-direction   

    ENDDO
    
  END SUBROUTINE shapeDerivatives_linear


END MODULE elementOps

