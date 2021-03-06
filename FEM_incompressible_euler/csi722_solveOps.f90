!!$ Program: csi722_solveOps.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 12/15/2011
!!$ Description: Module includes operations related to CSI 722 project 1.

MODULE solveOps
USE facedef
USE elementOps

  INTERFACE pressureCalc
     MODULE PROCEDURE pressureCalc_compress
  END INTERFACE

  INTERFACE calcGrad
     MODULE PROCEDURE calcGrad_linear
  END INTERFACE   

  INTERFACE solveLaplace
     MODULE PROCEDURE solveLaplace_linear
  END INTERFACE   

  INTERFACE eulerHalfStep
     MODULE PROCEDURE eulerHalfStep_compress,eulerHalfStep_incompress
  END INTERFACE   

  INTERFACE solveEuler
     MODULE PROCEDURE solveEuler_compress,solveEuler_incompress
  END INTERFACE   

  INTERFACE lapidusArtVis
     MODULE PROCEDURE lapidusArtVis_eulerCompress,lapidusArtVis_eulerIncompress
  END INTERFACE

  INTERFACE eulerBCs
     MODULE PROCEDURE eulerBCs_farCompress,eulerBCs_farIncompress
  END INTERFACE   

  INTERFACE timestepFEM
     MODULE PROCEDURE timestepFEM_eulerCompress,timestepFEM_eulerIncompress
  END INTERFACE      

  INTERFACE sumSq
     MODULE PROCEDURE sumSq_euler
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! pressureCalc_compress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE pressureCalc_compress(ndim,nval,nunks,un,Pr)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nval,nunks
    INTEGER :: i,j,k
    DOUBLE PRECISION :: v2,rho,vx,vy,vz,en
!!$    DOUBLE PRECISION,DIMENSION(npoin) :: rho  !density
!!$    DOUBLE PRECISION,DIMENSION(npoin) :: en  !energy
!!$    DOUBLE PRECISION,DIMENSION(ndim,npoin) :: vxyz !velocity vector
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,nval) :: un !unknowns
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nval) :: Pr !pressure
 
    !initialize pressure
    Pr(:) = 0.0d0

    DO i=1,nval

       rho = un(1,i)
       vx = un(2,i)/rho
       vy = un(3,i)/rho
       en = un(4,i)/rho

       v2 = vx*vx + vy*vy

       Pr(i) = (gamma - 1.0d0)*rho*(en - v2/2.0d0)

    ENDDO
    
  END SUBROUTINE pressureCalc_compress   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GRADIENT (LINEAR) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE calcGrad_linear(npoin,nnode,nelem,lnode,EA,xyz,fxyz,gsol)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: npoin,nnode,nelem
    INTEGER :: i,j,k
    INTEGER :: n1,n2,n3
    DOUBLE PRECISION :: u,v
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(npoin) :: fxyz  ! Scaler Field
    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: xyz ! Point Coordinates
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(3,npoin) :: gsol ! Gradient Vector
    DOUBLE PRECISION,DIMENSION(1:3) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(3,npoin) :: Nxyz ! Shape Derivatives
    DOUBLE PRECISION,DIMENSION(npoin) :: ML ! Mass Lumped Matrix
    DOUBLE PRECISION,DIMENSION(3,npoin) :: rxyz ! Right Hand Side
 
    ! Initialize Mass Lumped Matrix and Right Hand Side
    ML(:) = 0.0d0
    rxyz(:,:) = 0.0d0

    DO i=1,nelem
       n1 = lnode(1,i)
       n2 = lnode(2,i)
       n3 = lnode(3,i)

       ! Nodes of Current Element
       a(:) = xyz(:,n1)
       b(:) = xyz(:,n2)
       c(:) = xyz(:,n3)

       FORALL(j=1:3) Xba(j) = b(j) - a(j)
       FORALL(j=1:3) Xca(j) = c(j) - a(j)
       
       ! Calculate Shape Derivative for Element
       Nxyz(1,n1) = (-Xca(2) + Xba(2))/(2.0d0*EA(i))
       Nxyz(1,n2) = Xca(2)/(2.0d0*EA(i))
       Nxyz(1,n3) = -Xba(2)/(2.0d0*EA(i))

       Nxyz(2,n1) = (Xca(1) - Xba(1))/(2.0d0*EA(i))
       Nxyz(2,n2) = -Xca(1)/(2.0d0*EA(i))
       Nxyz(2,n3) = Xba(1)/(2.0d0*EA(i))     
     
       u=(EA(i)/3.0d0)*(Nxyz(1,n1)*fxyz(n1) + Nxyz(1,n2)*fxyz(n2) + Nxyz(1,n3)*fxyz(n3))
       v=(EA(i)/3.0d0)*(Nxyz(2,n1)*fxyz(n1) + Nxyz(2,n2)*fxyz(n2) + Nxyz(2,n3)*fxyz(n3))

       ! Update rhs for x and y coordinates
       rxyz(1,n1) = rxyz(1,n1) + u
       rxyz(1,n2) = rxyz(1,n2) + u
       rxyz(1,n3) = rxyz(1,n3) + u

       rxyz(2,n1) = rxyz(2,n1) + v
       rxyz(2,n2) = rxyz(2,n2) + v
       rxyz(2,n3) = rxyz(2,n3) + v

       ! Update Lumped Mass Matrix
       ML(n1) = ML(n1) + (EA(i)/3.0d0)
       ML(n2) = ML(n2) + (EA(i)/3.0d0)
       ML(n3) = ML(n3) + (EA(i)/3.0d0)

    ENDDO

    ! Solve for gradient field
    DO i=1,npoin
       gsol(1,i) = rxyz(1,i)/ML(i)
       gsol(2,i) = rxyz(2,i)/ML(i)
    ENDDO
    
  END SUBROUTINE calcGrad_linear   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LAPLACE (LINEAR) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE solveLaplace_linear(npoin,nelem,nnode,nbpS,nbpL,lnode,bnodeS,bnodeL,EA,xyz,fsol)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: npoin,nelem,nnode,nbpS,nbpL
    INTEGER :: i,j,k,kstart
    INTEGER :: n1,n2,n3,nj,nk,INFO
    DOUBLE PRECISION :: u,v
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    INTEGER,INTENT(IN),DIMENSION(nbpS) :: bnodeS
    INTEGER,INTENT(IN),DIMENSION(nbpL) :: bnodeL
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: xyz ! Point Coordinates
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(npoin) :: fsol ! Solution to Laplace
    INTEGER,DIMENSION(3) :: n_array ! Array for temp. storing nodes of element
    INTEGER,DIMENSION(npoin) :: IPIV
    DOUBLE PRECISION,DIMENSION(3) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(npoin,1) :: rxyz ! Right Hand Side
    DOUBLE PRECISION,DIMENSION(3,npoin) :: Nxyz ! Shape Derivatives
    DOUBLE PRECISION,DIMENSION(npoin,npoin) :: KL ! Matrix for Laplace Eq.

    rxyz(:,1) = 0.0d0
    KL(:,:) = 0.0d0

    DO i=1,nelem
       n1 = lnode(1,i)
       n2 = lnode(2,i)
       n3 = lnode(3,i)

       ! The elements of this array are n1,n2,n3 defined above
       n_array(:) = lnode(:,i)

       ! Nodes of Current Element
       a(:) = xyz(:,n1)
       b(:) = xyz(:,n2)
       c(:) = xyz(:,n3)

       Xba(:) = b(:) - a(:)
       Xca(:) = c(:) - a(:)
       
       ! Calculate Shape Derivative for Element
       Nxyz(1,n1) = (-Xca(2) + Xba(2))/(2.0d0*EA(i))
       Nxyz(1,n2) = Xca(2)/(2.0d0*EA(i))
       Nxyz(1,n3) = -Xba(2)/(2.0d0*EA(i))

       Nxyz(2,n1) = (Xca(1) - Xba(1))/(2.0d0*EA(i))
       Nxyz(2,n2) = -Xca(1)/(2.0d0*EA(i))
       Nxyz(2,n3) = Xba(1)/(2.0d0*EA(i))     

       ! Update Global Matrix
       kstart = 1
       DO j=1,nnode
          DO k=kstart,nnode
             nj = n_array(j)
             nk = n_array(k)
             KL(nj,nk) = KL(nj,nk) + EA(i)*(Nxyz(1,nj)*Nxyz(1,nk) + Nxyz(2,nj)*Nxyz(2,nk))
             KL(nk,nj) = KL(nj,nk)
          ENDDO
          kstart = kstart + 1
       ENDDO
    ENDDO

    ! Apply Side Boundary Conditions
    DO i=1,nbpS
       n1 = bnodeS(i)
       
       ! Side Boundaries: phi = 0
       rxyz(n1,1) = 0.0d0
       KL(n1,:) = 0.0d0
       KL(n1,n1) = 1.0d0

    ENDDO

    ! Apply Lower Boundary Conditions
    DO i=1,nbpL
       n1 = bnodeL(i)
       
       ! Lower Boundary: db/dx = -0.8x, integral over boundary is subtracted
       rxyz(n1,1) = rxyz(n1,1) + 0.80d0*xyz(1,n1)
       KL(n1,:) = 0.0d0
       KL(n1,n1) = 1.0d0

    ENDDO
    
    ! Solve System Using Lapack Library
    CALL DGESV(npoin,1,KL,npoin,IPIV,rxyz,npoin,INFO)
    fsol(:) = rxyz(:,1)
    
  END SUBROUTINE solveLaplace_linear


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER HALF-STEP (Compress) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerHalfStep_compress(ndim,npoin,nelem,nnode,nunks,lnode,dt,Nxyz,uxyz,uele)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    DOUBLE PRECISION :: vx,vy,vxavg,vyavg,enavg
    DOUBLE PRECISION :: Fax,Fay,Fbx,Fby,Fcx,Fcy
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: dti     !dt(i)/2
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: dt !element timestep
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nunks,nelem) :: uele !element unknowns
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    DOUBLE PRECISION,DIMENSION(ndim,npoin) :: vxyz !velocity vector
    DOUBLE PRECISION,DIMENSION(npoin) :: rho             !density
    DOUBLE PRECISION,DIMENSION(npoin) :: en              !energy
    DOUBLE PRECISION,DIMENSION(npoin) :: Pr              !pressure
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nunks) :: uavg,Fele

    Fax = 0.0d0
    Fay = 0.0d0
    Fbx = 0.0d0
    Fby = 0.0d0
    Fcx = 0.0d0
    Fcy = 0.0d0

    uele(:,:) = 0.0d0

    rho(:) = 0.0d0
    vxyz(:,:) = 0.0d0
    en(:) = 0.0d0
    Pr(:) = 0.0d0
    
    !calculate pressure
    CALL pressureCalc(ndim,npoin,nunks,uxyz,Pr)    

    !calculate quantities based off unknowns
    DO i=1,npoin
       rho(i) = uxyz(1,i)
       vxyz(1,i) = uxyz(2,i)/rho(i)
       vxyz(2,i) = uxyz(3,i)/rho(i)
       en(i) = uxyz(4,i)/rho(i)
    ENDDO

    !calculate flux at nodes
    DO i=1,nelem

       dti = dt(i)/2.0d0

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !calcualte average value of uxyz
       uavg(:) = 0.0d0
       uavg(:) = uavg(:) + uxyz(:,na)
       uavg(:) = uavg(:) + uxyz(:,nb)
       uavg(:) = uavg(:) + uxyz(:,nc)
       uavg(:) = uavg(:)/DBLE(nnode)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction 

       !calculate fluxes
       Fele(:) = 0.0d0

       !Flux for rho
       Fax = rho(na)*vxyz(1,na)
       Fay = rho(na)*vxyz(2,na)

       Fbx = rho(nb)*vxyz(1,nb)
       Fby = rho(nb)*vxyz(2,nb)

       Fcx = rho(nc)*vxyz(1,nc)
       Fcy = rho(nc)*vxyz(2,nc)

       Fele(1) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(1) = dti*Fele(1)

       !Flux for rho*vx
       Fax = rho(na)*vxyz(1,na)*vxyz(1,na) + Pr(na)
       Fay = rho(na)*vxyz(1,na)*vxyz(2,na)

       Fbx = rho(nb)*vxyz(1,nb)*vxyz(1,nb) + Pr(nb)
       Fby = rho(nb)*vxyz(1,nb)*vxyz(2,nb)

       Fcx = rho(nc)*vxyz(1,nc)*vxyz(1,nc) + Pr(nc)
       Fcy = rho(nc)*vxyz(1,nc)*vxyz(2,nc)

       Fele(2) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(2) = dti*Fele(2)

       !Flux for rho*vy
       Fax = rho(na)*vxyz(1,na)*vxyz(2,na) 
       Fay = rho(na)*vxyz(2,na)*vxyz(2,na) + Pr(na)

       Fbx = rho(nb)*vxyz(1,nb)*vxyz(2,nb) 
       Fby = rho(nb)*vxyz(2,nb)*vxyz(2,nb) + Pr(nb)

       Fcx = rho(nc)*vxyz(1,nc)*vxyz(2,nc) 
       Fcy = rho(nc)*vxyz(2,nc)*vxyz(2,nc) + Pr(nc)

       Fele(3) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(3) = dti*Fele(3)

       !Flux for rho*en
       Fax = vxyz(1,na)*(rho(na)*en(na) + Pr(na))
       Fay = vxyz(2,na)*(rho(na)*en(na) + Pr(na))

       Fbx = vxyz(1,nb)*(rho(nb)*en(nb) + Pr(nb))
       Fby = vxyz(2,nb)*(rho(nb)*en(nb) + Pr(nb))

       Fcx = vxyz(1,nc)*(rho(nc)*en(nc) + Pr(nc))
       Fcy = vxyz(2,nc)*(rho(nc)*en(nc) + Pr(nc))

       Fele(4) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(4) = dti*Fele(4)

       !calculate unknowns for elements
       uele(:,i) = uavg(:) - Fele(:)

    ENDDO

  END SUBROUTINE eulerHalfStep_compress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! solveEuler_compress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE solveEuler_compress(ndim,npoin,nelem,nnode,nface,nunks,lnode,bface,dt,EA,xyz,Nxyz,uxyz,uele,duxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nface,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION :: dx,dy
    DOUBLE PRECISION :: Fx,Fy
    DOUBLE PRECISION :: F1,F2,F3,F4
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: KLL,Lx,Ly
    DOUBLE PRECISION :: vx,vy,vz !velocity vector
    DOUBLE PRECISION :: rho      !density
    DOUBLE PRECISION :: en       !energy
    DOUBLE PRECISION,PARAMETER :: beta1 = 0.10d0
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: bface    !boundary faces
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: dt !element timestep array
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !point coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,nelem) :: uele !element unknowns
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nunks,npoin) :: duxyz !unknowns at point
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    INTEGER,DIMENSION(npoin) :: IPIV
    DOUBLE PRECISION,DIMENSION(npoin) :: dtp !timestep array
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nelem) :: Pr              !pressure
    DOUBLE PRECISION,DIMENSION(npoin) :: ML !lumped mass matrix
    DOUBLE PRECISION,DIMENSION(npoin,npoin) :: MC !consistent mass matrix
    DOUBLE PRECISION,DIMENSION(nunks) :: Fbnd !boundary fluxes
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs !Right Hand Side for rho,rho*vx,rho*vy,rho*e
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs_Lap  !Lapidus Artificial vis.
!!$    DOUBLE PRECISION,DIMENSION(npoin,1) :: rxyz,rtmp

    Pr(:) = 0.0d0
    !calculate pressure
    CALL pressureCalc(ndim,nelem,nunks,uele,Pr)    

    Fx = 0.0d0
    Fy = 0.0d0

    dtp(:) = 1E+16
    duxyz(:,:) = 0.0d0
    rhs(:,:) = 0.0d0
    rhs_Lap(:,:) = 0.0d0
    ML(:) = 0.0d0

    !calculate flux of element
    DO i=1,nelem

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       rho = uele(1,i)
       vx = uele(2,i)/rho
       vy = uele(3,i)/rho
       en = uele(4,i)/rho

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction 

       !scale shape derivatives by element area
       Nax = EA(i)*Nax
       Nay = EA(i)*Nay

       Nbx = EA(i)*Nbx
       Nby = EA(i)*Nby

       Ncx = EA(i)*Ncx
       Ncy = EA(i)*Ncy

       !Flux for rho
       Fx = rho*vx
       Fy = rho*vy

       rhs(1,na) = rhs(1,na) + Nax*Fx + Nay*Fy
       rhs(1,nb) = rhs(1,nb) + Nbx*Fx + Nby*Fy
       rhs(1,nc) = rhs(1,nc) + Ncx*Fx + Ncy*Fy

       !Flux for rho*vx
       Fx = rho*vx*vx + Pr(i)
       Fy = rho*vx*vy

       rhs(2,na) = rhs(2,na) + Nax*Fx + Nay*Fy
       rhs(2,nb) = rhs(2,nb) + Nbx*Fx + Nby*Fy
       rhs(2,nc) = rhs(2,nc) + Ncx*Fx + Ncy*Fy

       !Flux for rho*vy
       Fx = rho*vx*vy 
       Fy = rho*vy*vy + Pr(i)

       rhs(3,na) = rhs(3,na) + Nax*Fx + Nay*Fy
       rhs(3,nb) = rhs(3,nb) + Nbx*Fx + Nby*Fy
       rhs(3,nc) = rhs(3,nc) + Ncx*Fx + Ncy*Fy

       !Flux for rho*en
       Fx = vx*(rho*en + Pr(i))
       Fy = vy*(rho*en + Pr(i))

       rhs(4,na) = rhs(4,na) + Nax*Fx + Nay*Fy
       rhs(4,nb) = rhs(4,nb) + Nbx*Fx + Nby*Fy
       rhs(4,nc) = rhs(4,nc) + Ncx*Fx + Ncy*Fy

       !Update Lumped Mass Matrix
       ML(na) = ML(na) + (EA(i)/3.0d0)
       ML(nb) = ML(nb) + (EA(i)/3.0d0)
       ML(nc) = ML(nc) + (EA(i)/3.0d0)

       dtp(na) = MIN(dtp(na),dt(i))
       dtp(nb) = MIN(dtp(nb),dt(i))
       dtp(nc) = MIN(dtp(nc),dt(i))
    ENDDO

    a(:) = 0.0d0
    b(:) = 0.0d0
    Fbnd(:) = 0.0d0

    DO i=1,npoin
       dtp(i) = beta1*dtp(i)
    ENDDO

    !adjust right-hand side for boundary integral
    DO i=1,nface

       k = bface(i)%eleid
       na = bface(i)%fnode(1)
       nb = bface(i)%fnode(2)

       !calculate variables
       rho = uele(1,k)
       vx = uele(2,k)/rho
       vy = uele(3,k)/rho
       en = uele(4,k)/rho

       a(:) = xyz(:,na)
       b(:) = xyz(:,nb)

       !calculate unit normal
       dx = b(1) - a(1)
       dy = b(2) - a(2)
       
       !calculate flux at nodes
       Fx = rho*vx
       Fy = rho*vy
       
       !evaluate boundary flux (counterclockwise)
       Fbnd(1) = (Fx*dy - Fy*dx)/2.0d0
       
       !Flux for rho*vx
       Fx = rho*vx*vx + Pr(k)
       Fy = rho*vx*vy

       !evaluate boundary flux (counterclockwise)
       Fbnd(2) = (Fx*dy - Fy*dx)/2.0d0

       !Flux for rho*vy
       Fx = rho*vx*vy 
       Fy = rho*vy*vy + Pr(k)

       !evaluate boundary flux (counterclockwise)
       Fbnd(3) = (Fx*dy - Fy*dx)/2.0d0
       
       !Flux for rho*en
       Fx = vx*(rho*en + Pr(k))
       Fy = vy*(rho*en + Pr(k))

       !evaluate boundary flux (counterclockwise)
       Fbnd(4) = (Fx*dy - Fy*dx)/2.0d0
       
       !check if boundary is wall
       IF(bface(i)%btype .EQ. 0)THEN
          Fbnd(1) = 0.0d0
          Fbnd(2) = dy*Pr(k)/2.0d0
          Fbnd(3) = -dx*Pr(k)/2.0d0
          Fbnd(4) = 0.0d0
       ENDIF

       !adjust rhs
       rhs(1,na) = rhs(1,na) - Fbnd(1)
       rhs(2,na) = rhs(2,na) - Fbnd(2)
       rhs(3,na) = rhs(3,na) - Fbnd(3)
       rhs(4,na) = rhs(4,na) - Fbnd(4)
       rhs(1,nb) = rhs(1,nb) - Fbnd(1)
       rhs(2,nb) = rhs(2,nb) - Fbnd(2)
       rhs(3,nb) = rhs(3,nb) - Fbnd(3)
       rhs(4,nb) = rhs(4,nb) - Fbnd(4)

    ENDDO

    !calulate Lapidus artificial viscosity
    CALL lapidusArtVis(ndim,npoin,nelem,nnode,nunks,lnode,EA,xyz,Nxyz,uxyz,rhs_Lap)

    rhs(1,:) = rhs(1,:) + rhs_Lap(1,:)
    rhs(2,:) = rhs(2,:) + rhs_Lap(2,:)
    rhs(3,:) = rhs(3,:) + rhs_Lap(3,:)
    rhs(4,:) = rhs(4,:) + rhs_Lap(4,:)

    !Solve for unknowns
!!$    duxyz(:,:) = 0.0d0
    DO i=1,npoin
       duxyz(1,i) = dtp(i)*rhs(1,i)/ML(i)
       duxyz(2,i) = dtp(i)*rhs(2,i)/ML(i)
       duxyz(3,i) = dtp(i)*rhs(3,i)/ML(i)
       duxyz(4,i) = dtp(i)*rhs(4,i)/ML(i)
    ENDDO

  END SUBROUTINE solveEuler_compress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER HALF-STEP (Incompress) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerHalfStep_incompress(ndim,npoin,nelem,nnode,nunks,ca,lnode,dt,Nxyz,uxyz,uele)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    DOUBLE PRECISION,INTENT(IN) :: ca     !artifical speed of sound
    DOUBLE PRECISION :: ca2     !artifical speed of sound squared
    DOUBLE PRECISION :: vx,vy,vxavg,vyavg,enavg
    DOUBLE PRECISION :: Fax,Fay,Fbx,Fby,Fcx,Fcy
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: dti     !dt(i)/2
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: dt !element timestep
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nunks,nelem) :: uele !element unknowns
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nunks) :: uavg,Fele

    Fax = 0.0d0
    Fay = 0.0d0
    Fbx = 0.0d0
    Fby = 0.0d0
    Fcx = 0.0d0
    Fcy = 0.0d0

    uele(:,:) = 0.0d0

    ca2 = ca*ca

    !calculate flux at nodes
    DO i=1,nelem

       dti = dt(i)/2.0d0

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !calcualte average value of uxyz
       uavg(:) = 0.0d0
       uavg(:) = uavg(:) + uxyz(:,na)
       uavg(:) = uavg(:) + uxyz(:,nb)
       uavg(:) = uavg(:) + uxyz(:,nc)
       uavg(:) = uavg(:)/DBLE(nnode)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction 

       !calculate fluxes
       Fele(:) = 0.0d0

       !Flux for Pr
!!$       Fax = ca2*uxyz(2,na)
!!$       Fay = ca2*uxyz(3,na)
!!$
!!$       Fbx = ca2*uxyz(2,nb)
!!$       Fby = ca2*uxyz(3,nc)
!!$
!!$       Fcx = ca2*uxyz(2,nc)
!!$       Fcy = ca2*uxyz(3,nc)

       Fax = uxyz(2,na)
       Fay = uxyz(3,na)

       Fbx = uxyz(2,nb)
       Fby = uxyz(3,nc)

       Fcx = uxyz(2,nc)
       Fcy = uxyz(3,nc)

       Fele(1) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(1) = dti*Fele(1)

       !Flux for vx
!!$       Fax = uxyz(2,na)*uxyz(2,na) + uxyz(1,na)
!!$       Fay = uxyz(2,na)*uxyz(3,na)
!!$
!!$       Fbx = uxyz(2,nb)*uxyz(2,nb) + uxyz(1,nb)
!!$       Fby = uxyz(2,nb)*uxyz(3,nb)
!!$
!!$       Fcx = uxyz(2,nc)*uxyz(2,nc) + uxyz(1,nc)
!!$       Fcy = uxyz(2,nc)*uxyz(3,nc)
!!$
       Fax = uxyz(2,na)*uxyz(2,na) + ca2*uxyz(1,na)
       Fay = uxyz(2,na)*uxyz(3,na)

       Fbx = uxyz(2,nb)*uxyz(2,nb) + ca2*uxyz(1,nb)
       Fby = uxyz(2,nb)*uxyz(3,nb)

       Fcx = uxyz(2,nc)*uxyz(2,nc) + ca2*uxyz(1,nc)
       Fcy = uxyz(2,nc)*uxyz(3,nc)

       Fele(2) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(2) = dti*Fele(2)

       !Flux for vy
!!$       Fax = uxyz(2,na)*uxyz(3,na) 
!!$       Fay = uxyz(3,na)*uxyz(3,na) + uxyz(1,na)
!!$
!!$       Fbx = uxyz(2,nb)*uxyz(3,nb) 
!!$       Fby = uxyz(3,nc)*uxyz(3,nc) + uxyz(1,nc)
!!$
!!$       Fcx = uxyz(2,nc)*uxyz(3,nc) 
!!$       Fcy = uxyz(3,nc)*uxyz(3,nc) + uxyz(1,nc)
!!$
       Fax = uxyz(2,na)*uxyz(3,na) 
       Fay = uxyz(3,na)*uxyz(3,na) + ca2*uxyz(1,na)

       Fbx = uxyz(2,nb)*uxyz(3,nb) 
       Fby = uxyz(3,nc)*uxyz(3,nc) + ca2*uxyz(1,nc)

       Fcx = uxyz(2,nc)*uxyz(3,nc) 
       Fcy = uxyz(3,nc)*uxyz(3,nc) + ca2*uxyz(1,nc)

       Fele(3) = Nax*Fax + Nay*Fay + Nbx*Fbx + Nby*Fby + Ncx*Fcx + Ncy*Fcy
       Fele(3) = dti*Fele(3)

       !calculate unknowns for elements
       uele(:,i) = uavg(:) - Fele(:)

    ENDDO

  END SUBROUTINE eulerHalfStep_incompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! solveEuler_incompress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE solveEuler_incompress(ndim,npoin,nelem,nnode,nface,nunks,ca,lnode,bface,dt,EA,xyz,Nxyz,uxyz,uele,duxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nface,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION,INTENT(IN) :: ca     !artifical speed of sound
    DOUBLE PRECISION :: ca2     !artifical speed of sound squared
    DOUBLE PRECISION :: dx,dy
    DOUBLE PRECISION :: Fx,Fy
    DOUBLE PRECISION :: F1,F2,F3,F4
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: KLL,Lx,Ly
    DOUBLE PRECISION :: vx,vy,vz !velocity vector
    DOUBLE PRECISION,PARAMETER :: beta = 1.0d0
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: bface    !boundary faces
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: dt !element timestep array
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !point coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,nelem) :: uele !element unknowns
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nunks,npoin) :: duxyz !unknowns at point
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    INTEGER,DIMENSION(npoin) :: IPIV
    DOUBLE PRECISION,DIMENSION(npoin) :: dtp !timestep array
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(npoin) :: ML !lumped mass matrix
    DOUBLE PRECISION,DIMENSION(nunks) :: Fbnd !boundary fluxes
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs !Right Hand Side for rho,rho*vx,rho*vy,rho*e
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs_Lap  !Lapidus Artificial vis.

    Fx = 0.0d0
    Fy = 0.0d0

    dtp(:) = 1E+16
    duxyz(:,:) = 0.0d0
    rhs(:,:) = 0.0d0
    rhs_Lap(:,:) = 0.0d0
    ML(:) = 0.0d0

    ca2 = ca*ca

    !calculate flux of element
    DO i=1,nelem

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction 

       !scale shape derivatives by element area
       Nax = EA(i)*Nax
       Nay = EA(i)*Nay

       Nbx = EA(i)*Nbx
       Nby = EA(i)*Nby

       Ncx = EA(i)*Ncx
       Ncy = EA(i)*Ncy

       !Flux for Pr
!!$       Fx = ca2*uele(2,i)
!!$       Fy = ca2*uele(3,i)
       Fx = uele(2,i)
       Fy = uele(3,i)


       rhs(1,na) = rhs(1,na) + Nax*Fx + Nay*Fy
       rhs(1,nb) = rhs(1,nb) + Nbx*Fx + Nby*Fy
       rhs(1,nc) = rhs(1,nc) + Ncx*Fx + Ncy*Fy

       !Flux for vx
!!$       Fx = uele(2,i)*uele(2,i) + uele(1,i)
!!$       Fy = uele(2,i)*uele(3,i)
       Fx = uele(2,i)*uele(2,i) + ca2*uele(1,i)
       Fy = uele(2,i)*uele(3,i)

       rhs(2,na) = rhs(2,na) + Nax*Fx + Nay*Fy
       rhs(2,nb) = rhs(2,nb) + Nbx*Fx + Nby*Fy
       rhs(2,nc) = rhs(2,nc) + Ncx*Fx + Ncy*Fy

       !Flux for vy
!!$       Fx = uele(2,i)*uele(3,i) 
!!$       Fy = uele(3,i)*uele(3,i) + uele(1,i)
       Fx = uele(2,i)*uele(3,i) 
       Fy = uele(3,i)*uele(3,i) + ca2*uele(1,i)

       rhs(3,na) = rhs(3,na) + Nax*Fx + Nay*Fy
       rhs(3,nb) = rhs(3,nb) + Nbx*Fx + Nby*Fy
       rhs(3,nc) = rhs(3,nc) + Ncx*Fx + Ncy*Fy

       !Update Lumped Mass Matrix
       ML(na) = ML(na) + (EA(i)/3.0d0)
       ML(nb) = ML(nb) + (EA(i)/3.0d0)
       ML(nc) = ML(nc) + (EA(i)/3.0d0)

       !determine minimum timestep
       dtp(na) = MIN(dtp(na),dt(i))
       dtp(nb) = MIN(dtp(nb),dt(i))
       dtp(nc) = MIN(dtp(nc),dt(i))
    ENDDO

    a(:) = 0.0d0
    b(:) = 0.0d0
    Fbnd(:) = 0.0d0

    DO i=1,npoin
       dtp(i) = beta*dtp(i)
    ENDDO

    !adjust right-hand side for boundary integral
    DO i=1,nface

       k = bface(i)%eleid
       na = bface(i)%fnode(1)
       nb = bface(i)%fnode(2)

       a(:) = xyz(:,na)
       b(:) = xyz(:,nb)

       !calculate unit normal
       dx = b(1) - a(1)
       dy = b(2) - a(2)
       
       !calculate flux at nodes
       !flux for Pr
!!$       Fx = ca2*uele(2,k)
!!$       Fy = ca2*uele(3,k)
       Fx = uele(2,k)
       Fy = uele(3,k)
       
       !evaluate boundary flux (counterclockwise)
       Fbnd(1) = (Fx*dy - Fy*dx)/2.0d0
       
       !Flux for vx
!!$       Fx = uele(2,k)*uele(2,k) + uele(1,k)
!!$       Fy = uele(2,k)*uele(3,k)
       Fx = uele(2,k)*uele(2,k) + ca2*uele(1,k)
       Fy = uele(2,k)*uele(3,k)

       !evaluate boundary flux (counterclockwise)
       Fbnd(2) = (Fx*dy - Fy*dx)/2.0d0

       !Flux for vy
!!$       Fx = uele(2,k)*uele(3,k) 
!!$       Fy = uele(3,k)*uele(3,k) + uele(1,k)
       Fx = uele(2,k)*uele(3,k) 
       Fy = uele(3,k)*uele(3,k) + ca2*uele(1,k)

       !evaluate boundary flux (counterclockwise)
       Fbnd(3) = (Fx*dy - Fy*dx)/2.0d0
       
       !check if boundary is wall
       IF(bface(i)%btype .EQ. 0)THEN
          Fbnd(1) = 0.0d0
!!$          Fbnd(2) = dy*uele(1,k)/2.0d0
!!$          Fbnd(3) = -dx*uele(1,k)/2.0d0
          Fbnd(2) = dy*ca2*uele(1,k)/2.0d0
          Fbnd(3) = -dx*ca2*uele(1,k)/2.0d0
       ENDIF

       !adjust rhs
       rhs(1,na) = rhs(1,na) - Fbnd(1)
       rhs(2,na) = rhs(2,na) - Fbnd(2)
       rhs(3,na) = rhs(3,na) - Fbnd(3)

       rhs(1,nb) = rhs(1,nb) - Fbnd(1)
       rhs(2,nb) = rhs(2,nb) - Fbnd(2)
       rhs(3,nb) = rhs(3,nb) - Fbnd(3)

    ENDDO

    !calulate Lapidus artificial viscosity
    CALL lapidusArtVis(ndim,npoin,nelem,nnode,nunks,ca,lnode,EA,xyz,Nxyz,uxyz,rhs_Lap)

    rhs(1,:) = rhs(1,:) + rhs_Lap(1,:)
    rhs(2,:) = rhs(2,:) + rhs_Lap(2,:)
    rhs(3,:) = rhs(3,:) + rhs_Lap(3,:)

    !Solve for unknowns
!!$    duxyz(:,:) = 0.0d0
    DO i=1,npoin
       duxyz(1,i) = dtp(i)*rhs(1,i)/ML(i)
       duxyz(2,i) = dtp(i)*rhs(2,i)/ML(i)
       duxyz(3,i) = dtp(i)*rhs(3,i)/ML(i)
    ENDDO

  END SUBROUTINE solveEuler_incompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! lapidusArtVis_eulerCompress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE lapidusArtVis_eulerCompress(ndim,npoin,nelem,nnode,nunks,lnode,EA,xyz,Nxyz,uxyz,rhs_Lap)
    !Artificial Viscosity: Lapidus
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION,PARAMETER :: c1 = 4.0d0
    DOUBLE PRECISION,PARAMETER :: tol = 1E-16
    DOUBLE PRECISION :: dx,dy
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: NaL,NbL,NcL
    DOUBLE PRECISION :: vLa,vLb,vLc         
    DOUBLE PRECISION :: ua,ub,uc,duL         
    DOUBLE PRECISION :: Lx,Ly 
    DOUBLE PRECISION :: KLL
    DOUBLE PRECISION :: he2 !height of element squared
    DOUBLE PRECISION :: MLinv !inverse lumped mass matrix
    DOUBLE PRECISION :: magdVa,magdVb,magdVc       !magnitude of grad. velocity
    DOUBLE PRECISION :: dmagVx,dmagVy
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !point coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nunks,npoin) :: rhs_Lap  !Lapidus Artificial vis.
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    INTEGER,DIMENSION(npoin) :: IPIV
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nnode) :: magdV !magnitude of gradient of velocity
    DOUBLE PRECISION,DIMENSION(nunks) :: VLap
    DOUBLE PRECISION,DIMENSION(nelem) :: Pr              !pressure
    DOUBLE PRECISION,DIMENSION(npoin) :: rho             !density
    DOUBLE PRECISION,DIMENSION(npoin) :: vx 
    DOUBLE PRECISION,DIMENSION(npoin) :: vy 
    DOUBLE PRECISION,DIMENSION(npoin) :: magV
    DOUBLE PRECISION,DIMENSION(npoin) :: en              !energy

    rho(:) = 0.0d0
    vx(:) = 0.0d0
    vy(:) = 0.0d0
    en(:) = 0.0d0
    DO i=1,npoin
       rho(i) = uxyz(1,i)
       vx(i) = uxyz(2,i)/rho(i)
       vy(i) = uxyz(3,i)/rho(i)
       en(i) = uxyz(4,i)/rho(i)
       magV(i) = DSQRT(vx(i)*vx(i) + vy(i)*vy(i))
    ENDDO

    rhs_Lap(:,:) = 0.0d0

    !calculate Lapidus artificial viscosity
    DO i=1,nelem

       he2 = EA(i)

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction

       dmagVx = Nax*magV(na) + Nbx*magV(nb) + Ncx*magV(nc)
       dmagVy = Nay*magV(na) + Nby*magV(nb) + Ncy*magV(nc)

       IF((dmagVx .GT. tol) .OR. (dmagVy .GT. tol))THEN       
          Lx = dmagVx/DSQRT(dmagVx*dmagVx + dmagVy*dmagVy)
          Ly = dmagVy/DSQRT(dmagVx*dmagVx + dmagVy*dmagVy)
       ELSE
          Lx = 0.0d0
          Ly = 0.0d0
       ENDIF

       !transform shape function derivatives to be w.r.t. Lx,Ly
       NaL = Nax*Lx + Nay*Ly
       NbL = Nbx*Lx + Nby*Ly
       NcL = Ncx*Lx + Ncy*Ly

       KLL = 0.0d0
       !some over all nodes
       KLL =  NaL*(vx(na)*Lx + vy(na)*Ly) + NbL*(vx(nb)*Lx + vy(nb)*Ly) + NcL*(vx(nc)*Lx + vy(nc)*Ly)
       KLL = c1*he2*KLL

       !calc. for each unknown
       VLap(:) = NaL*uxyz(:,na) + NbL*uxyz(:,nb) + NcL*uxyz(:,nc)
       
       !Multiply by ABS(KLL)
       VLap(:) = EA(i)*ABS(KLL)*VLap(:) 
       
       rhs_Lap(:,na) = rhs_Lap(:,na) - NaL*VLap(:)
       rhs_Lap(:,nb) = rhs_Lap(:,nb) - NbL*Vlap(:)
       rhs_Lap(:,nc) = rhs_Lap(:,nc) - NcL*Vlap(:)
!!$       IF(i .EQ. 0)THEN
!!$          WRITE(*,*)' for element: ',i
!!$          WRITE(*,*)' with nodes: ',na,nb,nc
!!$          WRITE(*,*)' unkno(na)= ',(uxyz(j,na),j=1,nunks)
!!$          WRITE(*,*)' unkno(nb)= ',(uxyz(j,nb),j=1,nunks)
!!$          WRITE(*,*)' unkno(nc)= ',(uxyz(j,nc),j=1,nunks)
!!$          WRITE(*,*)' Lx,Ly= ',Lx,Ly
!!$          WRITE(*,*)' Nax,Nay= ',Nax,Nay
!!$          WRITE(*,*)' Nbx,Nby= ',Nbx,Nby
!!$          WRITE(*,*)' Ncx,Ncy= ',Ncx,Ncy
!!$          WRITE(*,*)' NaL,NbL,NcL= ',NaL,NbL,NcL
!!$          WRITE(*,*)' KLL= ',KLL
!!$          WRITE(*,*)' EA,Vlap= ',EA(i),(Vlap(j),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_a= ',(rhs_Lap(j,na),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_b= ',(rhs_Lap(j,nb),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_c= ',(rhs_Lap(j,nc),j=1,nunks)
!!$       ENDIF
    ENDDO

  END SUBROUTINE lapidusArtVis_eulerCompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! lapidusArtVis_eulerIncompress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE lapidusArtVis_eulerIncompress(ndim,npoin,nelem,nnode,nunks,ca,lnode,EA,xyz,Nxyz,uxyz,rhs_Lap)
    !Artificial Viscosity: Lapidus
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION,INTENT(IN) :: ca     !artifical speed of sound
    DOUBLE PRECISION :: ca2     !artifical speed of sound squared
    DOUBLE PRECISION,PARAMETER :: c1 = 1.0d0
    DOUBLE PRECISION,PARAMETER :: tol = 1E-16
    DOUBLE PRECISION :: dx,dy
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: NaL,NbL,NcL
    DOUBLE PRECISION :: vLa,vLb,vLc         
    DOUBLE PRECISION :: ua,ub,uc,duL         
    DOUBLE PRECISION :: Lx,Ly 
    DOUBLE PRECISION :: KLL
    DOUBLE PRECISION :: he2 !height of element squared
    DOUBLE PRECISION :: MLinv !inverse lumped mass matrix
    DOUBLE PRECISION :: magdVa,magdVb,magdVc       !magnitude of grad. velocity
    DOUBLE PRECISION :: dmagVx,dmagVy
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA ! Area of Element
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !point coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nunks,npoin) :: rhs_Lap  !Lapidus Artificial vis.
    INTEGER,DIMENSION(nnode) :: n_array ! array for temp. storing nodes of element
    INTEGER,DIMENSION(npoin) :: IPIV
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca
    DOUBLE PRECISION,DIMENSION(nnode) :: magdV !magnitude of gradient of velocity
    DOUBLE PRECISION,DIMENSION(nunks) :: VLap
!!$    DOUBLE PRECISION,DIMENSION(nelem) :: Pr              !pressure
    DOUBLE PRECISION,DIMENSION(npoin) :: vx 
    DOUBLE PRECISION,DIMENSION(npoin) :: vy 
    DOUBLE PRECISION,DIMENSION(npoin) :: magV

    DO i=1,npoin
       magV(i) = DSQRT(uxyz(2,i)*uxyz(2,i) + uxyz(3,i)*uxyz(3,i))
    ENDDO

    rhs_Lap(:,:) = 0.0d0

    !calculate Lapidus artificial viscosity
    DO i=1,nelem

       he2 = EA(i)

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction

       dmagVx = Nax*magV(na) + Nbx*magV(nb) + Ncx*magV(nc)
       dmagVy = Nay*magV(na) + Nby*magV(nb) + Ncy*magV(nc)

       IF((dmagVx .GT. tol) .OR. (dmagVy .GT. tol))THEN       
          Lx = dmagVx/DSQRT(dmagVx*dmagVx + dmagVy*dmagVy)
          Ly = dmagVy/DSQRT(dmagVx*dmagVx + dmagVy*dmagVy)
       ELSE
          Lx = 0.0d0
          Ly = 0.0d0
       ENDIF

       !transform shape function derivatives to be w.r.t. Lx,Ly
       NaL = Nax*Lx + Nay*Ly
       NbL = Nbx*Lx + Nby*Ly
       NcL = Ncx*Lx + Ncy*Ly

       KLL = 0.0d0
       !some over all nodes
       KLL =  NaL*(uxyz(2,na)*Lx + uxyz(3,na)*Ly) + NbL*(uxyz(2,nb)*Lx + uxyz(3,nb)*Ly) + NcL*(uxyz(2,nc)*Lx + uxyz(3,nc)*Ly)
       KLL = c1*he2*KLL

       !calc. for each unknown
       VLap(:) = NaL*uxyz(:,na) + NbL*uxyz(:,nb) + NcL*uxyz(:,nc)
       
       !Multiply by ABS(KLL)
       VLap(:) = EA(i)*ABS(KLL)*VLap(:) 
       
       rhs_Lap(:,na) = rhs_Lap(:,na) - NaL*VLap(:)
       rhs_Lap(:,nb) = rhs_Lap(:,nb) - NbL*Vlap(:)
       rhs_Lap(:,nc) = rhs_Lap(:,nc) - NcL*Vlap(:)
!!$       IF(i .EQ. 0)THEN
!!$          WRITE(*,*)' for element: ',i
!!$          WRITE(*,*)' with nodes: ',na,nb,nc
!!$          WRITE(*,*)' unkno(na)= ',(uxyz(j,na),j=1,nunks)
!!$          WRITE(*,*)' unkno(nb)= ',(uxyz(j,nb),j=1,nunks)
!!$          WRITE(*,*)' unkno(nc)= ',(uxyz(j,nc),j=1,nunks)
!!$          WRITE(*,*)' Lx,Ly= ',Lx,Ly
!!$          WRITE(*,*)' Nax,Nay= ',Nax,Nay
!!$          WRITE(*,*)' Nbx,Nby= ',Nbx,Nby
!!$          WRITE(*,*)' Ncx,Ncy= ',Ncx,Ncy
!!$          WRITE(*,*)' NaL,NbL,NcL= ',NaL,NbL,NcL
!!$          WRITE(*,*)' KLL= ',KLL
!!$          WRITE(*,*)' EA,Vlap= ',EA(i),(Vlap(j),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_a= ',(rhs_Lap(j,na),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_b= ',(rhs_Lap(j,nb),j=1,nunks)
!!$          WRITE(*,*)' rhs_Lap_c= ',(rhs_Lap(j,nc),j=1,nunks)
!!$       ENDIF
    ENDDO

  END SUBROUTINE lapidusArtVis_eulerIncompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER BCs (COMPRESS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerBCs_farCompress(ndim,npoin,nelem,nnode,nboun,nunks,lnode,bnode,xyz,ulast,uxyz)
    !Farfield boundary conditions
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nboun,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    TYPE(NODE),INTENT(IN),DIMENSION(nboun) :: bnode           !boundary nodes
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !xyz coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: ulast !unknowns at nodes previous time
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION :: vx,vy,vz !velocity vector
    DOUBLE PRECISION :: rho,rholast !density
    DOUBLE PRECISION :: en              !energy
    DOUBLE PRECISION :: Pr              !pressure
    DOUBLE PRECISION :: vn,vt,v2 !normal and tangential component
    DOUBLE PRECISION :: vnp,Prp
    DOUBLE PRECISION :: cs,cslast    !speed of sound
    DOUBLE PRECISION :: rhocs_avg,cs_avg     !average values
    DOUBLE PRECISION :: vninf,vtinf     !infinity values velocity
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca,nhat,that

    nhat(:) = 0.0d0
    that(:) = 0.0d0
    vnp = 0.0d0
    Prp = 0.0d0

    DO i=1,nboun

       k = bnode(i)%bpt
       vninf = bnode(i)%vninf
       vtinf = bnode(i)%vtinf

       nhat(:) = bnode(i)%nhat(:)
       that(:) = bnode(i)%that(:)

       !calculate variables at previous timestep
       rholast = ulast(1,k)
       vx = ulast(2,k)/rholast
       vy = ulast(3,k)/rholast
       en = ulast(4,k)/rholast

       v2 = vx*vx + vy*vy
       
       !calculate pressure
       Pr = (gamma - 1.0d0)*rholast*(en - v2/2.0d0)
       
       !calculate speed of sound at previous timestep
       cslast = DSQRT(gamma*Pr/rholast)
       
       !calculate variables
       rho = uxyz(1,k)
       vx = uxyz(2,k)/rho
       vy = uxyz(3,k)/rho
       en = uxyz(4,k)/rho

       v2 = vx*vx + vy*vy
       
       !calculate pressure
       Pr = (gamma - 1.0d0)*rho*(en - v2/2.0d0)
       
       !calculate speed of sound
       cs = DSQRT(gamma*Pr/rho)
       
       !calculate average values
       cs_avg = (cslast + cs)/2.0d0
       rhocs_avg = (rholast*cslast + rho*cs)/2.0d0
       
       !define normal and tangent components
       vn = vx*nhat(1) + vy*nhat(2)
       vt = vx*that(1) + vy*that(2)

       !determine boundry condtions
       SELECT CASE(bnode(i)%btype)
          !farfield conditions
       CASE(4)
          !determine inflow or outflow
          SELECT CASE(bnode(i)%flow)
             !inflow
             CASE(0)

!!$          IF(vy .NE. vyinf)THEN
!!$             PRINT*,'1. ',i,vy,vyinf
!!$          ENDIF

                vnp = vn        !predicted
                Prp = Pr        !predicted
                Pr = (Prinf + Prp + rhocs_avg*(vninf - vnp))/2.0d0
                rho = rhoinf + (Pr - Prinf)/(cs_avg*cs_avg)
                vn = vninf + (Prinf - Pr)/rhocs_avg
                vt = vtinf

!!$          vx = vt*that(1) + vn*nhat(1)
!!$          vy = vt*that(2) + vn*nhat(2)
!!$
!!$          v2 = vx*vx + vy*vy
!!$
!!$          en = (Pr/((gamma - 1.0d0)*rho)) + v2/2.0d0
!!$          !correct unknown variables
!!$          uxyz(1,k) = rho
!!$          uxyz(2,k) = vx*rho
!!$          uxyz(3,k) = vy*rho
!!$          uxyz(4,k) = en*rho

!!$          IF(vy .NE. vyinf)THEN
!!$             PRINT*,'2. ',i,vy,vyinf
!!$          ENDIF
          v2 = vxinf*vxinf + vyinf*vyinf
          uxyz(1,k) = rhoinf
          uxyz(2,k) = vxinf*rhoinf
          uxyz(3,k) = vyinf*rhoinf
          uxyz(4,k) = ((Prinf/((gamma - 1.0d0)*rhoinf)) + v2/2.0d0)*rhoinf

             !outflow
             CASE(1)
                vnp = vn        !predicted
                Prp = Pr        !predicted
                vn = vninf ! = cs*Minf*DCOS(alpha)
                vt = vt
                Pr = Prp + rho*cs*(vn - vnp)
                rho = rho + (Pr - Prp)/(cs_avg*cs_avg)   

!!$          IF(vt .NE. vtinf)THEN
!!$             PRINT*,''
!!$             PRINT*,i,'Pr = ',Pr,'Prp =',Prp
!!$             PRINT*,i,'en = ',en,'rho=',rho
!!$             PRINT*,i,'vn = ',vn,'vninf=',vninf
!!$             PRINT*,i,'vt = ',vt,'vtinf=',vtinf
!!$          ENDIF

!!$          vx = vt*that(1) + vn*nhat(1)
!!$          vy = vt*that(2) + vn*nhat(2)
!!$
!!$          v2 = vx*vx + vy*vy
!!$
!!$          en = (Pr/((gamma - 1.0d0)*rho)) + v2/2.0d0
!!$          !correct unknown variables
!!$          uxyz(1,k) = rho
!!$          uxyz(2,k) = vx*rho
!!$          uxyz(3,k) = vy*rho
!!$          uxyz(4,k) = en*rho

          v2 = vxinf*vxinf + vyinf*vyinf
          uxyz(1,k) = rhoinf
          uxyz(2,k) = vxinf*rhoinf
          uxyz(3,k) = vyinf*rhoinf
          uxyz(4,k) = ((Prinf/((gamma - 1.0d0)*rhoinf)) + v2/2.0d0)*rhoinf


           ENDSELECT

        CASE(0)
           !wall boundary conditions

                !determine wing tip
                SELECT CASE(bnode(i)%flow)
                   !inflow
                   CASE(-1)
                      vn = 0.0d0
                      vt = vt

          !x,y components of velocity
          vx = vt*that(1) + vn*nhat(1)
          vy = vt*that(2) + vn*nhat(2)

          v2 = vx*vx + vy*vy

          en = (Pr/((gamma - 1.0d0)*rho)) + v2/2.0d0

          !correct unknown variables
          uxyz(1,k) = rho
          uxyz(2,k) = vx*rho
          uxyz(3,k) = vy*rho
          uxyz(4,k) = en*rho

                   !outflow
                   CASE(-2)
!!$                      vt = 0.0d0

                ENDSELECT
          ENDSELECT
          
     ENDDO

  END SUBROUTINE eulerBCs_farCompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER BCs (INCOMPRESS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerBCs_farIncompress(ndim,npoin,nelem,nnode,nboun,nunks,ca,lnode,bnode,xyz,uxyz)
    !Farfield boundary conditions incompressible Euler Eqs.
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nboun,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    DOUBLE PRECISION,INTENT(IN) :: ca     !artifical speed of sound
    DOUBLE PRECISION :: ca2_inv     !artifical speed of sound squared
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    TYPE(NODE),INTENT(IN),DIMENSION(nboun) :: bnode           !boundary nodes
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !xyz coordinates
    DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION :: vx,vy,vz !velocity vector
    DOUBLE PRECISION :: Pr              !pressure
    DOUBLE PRECISION :: vn,vt,v2 !normal and tangential component
    DOUBLE PRECISION :: vnp,Prp
    DOUBLE PRECISION :: vninf,vtinf     !infinity values velocity
    DOUBLE PRECISION,DIMENSION(ndim) :: a,b,c,Xba,Xca,nhat,that

    nhat(:) = 0.0d0
    that(:) = 0.0d0
    vnp = 0.0d0
    Prp = 0.0d0

    ca2_inv = 1/(ca*ca)

    DO i=1,nboun

       k = bnode(i)%bpt
       nhat(:) = bnode(i)%nhat(:)
       that(:) = bnode(i)%that(:)

       !define normal and tangent components
       vn = uxyz(2,k)*nhat(1) + uxyz(3,k)*nhat(2)
       vt = uxyz(2,k)*that(1) + uxyz(3,k)*that(2)

       !determine boundry condtions
       SELECT CASE(bnode(i)%btype)
          !farfield conditions
       CASE(4)

          uxyz(1,k) = Prinf
          uxyz(2,k) = vxinf
          uxyz(3,k) = vyinf

       CASE(0)
          !wall boundary conditions

          !determine wing tip
          SELECT CASE(bnode(i)%flow)
             !inflow
          CASE(-1)
             vn = 0.0d0
             vt = vt
             
             !x,y components of velocity
             vx = vt*that(1) + vn*nhat(1)
             vy = vt*that(2) + vn*nhat(2)
             
             v2 = vx*vx + vy*vy
             
             !correct unknown variables
!!$             uxyz(1,k) = Prinf
!!$             uxyz(2,k) = vx
!!$             uxyz(3,k) = vy
             
             !outflow
          CASE(-2)
!!$                      vt = 0.0d0
             
          ENDSELECT
       ENDSELECT
          
     ENDDO

  END SUBROUTINE eulerBCs_farIncompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Timestep FEM (EULER COMPRESS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE timestepFEM_eulerCompress(ndim,npoin,nelem,nnode,nunks,lnode,EA,xyz,Nxyz,uxyz,dt)
    !Calculate allowable timestep finite element method
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    DOUBLE PRECISION,PARAMETER :: beta = 0.90d0
    DOUBLE PRECISION :: he !height of element
    DOUBLE PRECISION :: cs    !speed of sound
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy,Nmag
    DOUBLE PRECISION :: rho !density
    DOUBLE PRECISION :: vx,vy
    DOUBLE PRECISION :: nhatx,nhaty !unit normals
    DOUBLE PRECISION :: vn 
    DOUBLE PRECISION :: lambdaMax    
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: EA !element area
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !xyz coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nelem) :: dt !timestep
    DOUBLE PRECISION,DIMENSION(npoin) :: Pr !pressure at nodes
    DOUBLE PRECISION,DIMENSION(nnode) :: lambda !characteristics

    Pr(:) = 0.0d0
    !calculate pressure
    CALL pressureCalc(ndim,npoin,nunks,uxyz,Pr)    

    dt(:) = 0.0d0
    DO i=1,nelem

       lambda(:) = 0.0d0

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       he = DSQRT(EA(i))

       !calculate lambda for first node
       rho = uxyz(1,na)
       vx = uxyz(2,na)/rho
       vy = uxyz(3,na)/rho

       cs = DSQRT(gamma*Pr(na)/rho)

       lambda(1) = DSQRT(vx*vx + vy*vy)
       lambda(1) = lambda(1) + cs

       !calculate lambda for second node
       rho = uxyz(1,nb)
       vx = uxyz(2,nb)/rho
       vy = uxyz(3,nb)/rho

       cs = DSQRT(gamma*Pr(nb)/rho)

       lambda(2) = DSQRT(vx*vx + vy*vy)
       lambda(2) = lambda(2) + cs

       rho = uxyz(1,nc)
       vx = uxyz(2,nc)/rho
       vy = uxyz(3,nc)/rho

       cs = DSQRT(gamma*Pr(nc)/rho)

       !calculate lambda for third node
       lambda(3) = DSQRT(vx*vx + vy*vy)
       lambda(3) = lambda(3) + cs

       lambdaMax = MAXVAL(lambda(:))

       dt(i) = beta*he/lambdaMax

    ENDDO

  END SUBROUTINE timestepFEM_eulerCompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Timestep FEM (EULER INCOMPRESS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE timestepFEM_eulerIncompress(ndim,npoin,nelem,nnode,nunks,ca,lnode,xyz,Nxyz,uxyz,dt)
    !Calculate allowable timestep finite element method
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nelem,nnode,nunks
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk
    DOUBLE PRECISION,INTENT(IN) :: ca     !artifical speed of sound
    DOUBLE PRECISION,PARAMETER :: beta = 0.10d0
    DOUBLE PRECISION :: alpha   !constant, vxyz dot Nxyz
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy,Nmag
    DOUBLE PRECISION :: vx,vy
    DOUBLE PRECISION :: vn 
    DOUBLE PRECISION :: lambdaMax    
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !xyz coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nelem) :: dt !timestep
    DOUBLE PRECISION,DIMENSION(nnode) :: lambda !characteristics

    dt(:) = 0.0d0
    DO i=1,nelem

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       !shape derivatives
       Nax = Nxyz(1,i)          !node a of element i, x-direction
       Nay = Nxyz(2,i)          !node a of element i, y-direction

       Nbx = Nxyz(3,i)          !node b of element i, x-direction
       Nby = Nxyz(4,i)          !node b of element i, y-direction

       Ncx = Nxyz(5,i)          !node c of element i, x-direction
       Ncy = Nxyz(6,i)          !node c of element i, y-direction 

       !calculate lambda for node a
       lambda(1) = DABS(uxyz(2,na)*Nax + uxyz(3,na)*Nay)

       !calculate lambda for node b
       lambda(2) = DABS(uxyz(2,nb)*Nbx + uxyz(3,nb)*Nby)

       !calculate lambda for node c
       lambda(3) = DABS(uxyz(2,nc)*Ncx + uxyz(3,nc)*Ncy)

       lambdaMax = MAXVAL(lambda(:)) + ca

       dt(i) = beta/lambdaMax

    ENDDO

  END SUBROUTINE timestepFEM_eulerIncompress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sum of Squares !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE sumSq_euler(k_in,nunks,npoin,duxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: k_in,npoin,nunks
    INTEGER :: i,j,k
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: duxyz
    DOUBLE PRECISION,DIMENSION(nunks) :: du2    

    du2(:) = 0.0d0

    DO i=1,nunks
       du2(i) = 0.0d0
       DO j=1,npoin
         du2(i) = du2(i) + duxyz(i,j)*duxyz(i,j)
       ENDDO
    ENDDO

    WRITE(*,*),k_in,(du2(j),j=1,nunks)

  END SUBROUTINE sumSq_euler


END MODULE solveOps

