!!$ Program: csi721_solveOps.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 04/01/2011
!!$ Description: Module includes operations related to CSI 721 project 2.

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
     MODULE PROCEDURE eulerHalfStep_linear
  END INTERFACE   

  INTERFACE solveEuler
     MODULE PROCEDURE solveEuler_compress
  END INTERFACE   

  INTERFACE lapidusArtVis
     MODULE PROCEDURE lapidusArtVis_euler
  END INTERFACE

  INTERFACE eulerBCs
     MODULE PROCEDURE eulerBCs_farfield
  END INTERFACE   

  INTERFACE timestepFEM
     MODULE PROCEDURE timestepFEM_euler
  END INTERFACE      

  INTERFACE sumSq
     MODULE PROCEDURE sumSq_euler
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! pressureCalc_compress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE pressureCalc_compress(ndim,nval,un,Pr)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nval
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
  SUBROUTINE calcGrad_linear(lnode,elvol,xyz,fxyz,gsol)
    IMPLICIT NONE

    INTEGER :: i,j,k
    INTEGER :: n1,n2,n3
    DOUBLE PRECISION :: u,v
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol ! Area of Element
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
       Nxyz(1,n1) = (-Xca(2) + Xba(2))/(2.0d0*elvol(i))
       Nxyz(1,n2) = Xca(2)/(2.0d0*elvol(i))
       Nxyz(1,n3) = -Xba(2)/(2.0d0*elvol(i))

       Nxyz(2,n1) = (Xca(1) - Xba(1))/(2.0d0*elvol(i))
       Nxyz(2,n2) = -Xca(1)/(2.0d0*elvol(i))
       Nxyz(2,n3) = Xba(1)/(2.0d0*elvol(i))     
     
       u=(elvol(i)/3.0d0)*(Nxyz(1,n1)*fxyz(n1) + Nxyz(1,n2)*fxyz(n2) + Nxyz(1,n3)*fxyz(n3))
       v=(elvol(i)/3.0d0)*(Nxyz(2,n1)*fxyz(n1) + Nxyz(2,n2)*fxyz(n2) + Nxyz(2,n3)*fxyz(n3))

       ! Update rhs for x and y coordinates
       rxyz(1,n1) = rxyz(1,n1) + u
       rxyz(1,n2) = rxyz(1,n2) + u
       rxyz(1,n3) = rxyz(1,n3) + u

       rxyz(2,n1) = rxyz(2,n1) + v
       rxyz(2,n2) = rxyz(2,n2) + v
       rxyz(2,n3) = rxyz(2,n3) + v

       ! Update Lumped Mass Matrix
       ML(n1) = ML(n1) + (elvol(i)/3.0d0)
       ML(n2) = ML(n2) + (elvol(i)/3.0d0)
       ML(n3) = ML(n3) + (elvol(i)/3.0d0)

    ENDDO

    ! Solve for gradient field
    DO i=1,npoin
       gsol(1,i) = rxyz(1,i)/ML(i)
       gsol(2,i) = rxyz(2,i)/ML(i)
    ENDDO
    
  END SUBROUTINE calcGrad_linear   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LAPLACE (LINeaR) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE solveLaplace_linear(nbpS,nbpL,lnode,bnodeS,bnodeL,elvol,xyz,fsol)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nbpS,nbpL
    INTEGER :: i,j,k,kstart
    INTEGER :: n1,n2,n3,nj,nk,INFO
    DOUBLE PRECISION :: u,v
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    INTEGER,INTENT(IN),DIMENSION(nbpS) :: bnodeS
    INTEGER,INTENT(IN),DIMENSION(nbpL) :: bnodeL
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol ! Area of Element
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
       Nxyz(1,n1) = (-Xca(2) + Xba(2))/(2.0d0*elvol(i))
       Nxyz(1,n2) = Xca(2)/(2.0d0*elvol(i))
       Nxyz(1,n3) = -Xba(2)/(2.0d0*elvol(i))

       Nxyz(2,n1) = (Xca(1) - Xba(1))/(2.0d0*elvol(i))
       Nxyz(2,n2) = -Xca(1)/(2.0d0*elvol(i))
       Nxyz(2,n3) = Xba(1)/(2.0d0*elvol(i))     

       ! Update Global Matrix
       kstart = 1
       DO j=1,nnode
          DO k=kstart,nnode
             nj = n_array(j)
             nk = n_array(k)
             KL(nj,nk) = KL(nj,nk) + elvol(i)*(Nxyz(1,nj)*Nxyz(1,nk) + Nxyz(2,nj)*Nxyz(2,nk))
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER HALF-STEP (LINeaR) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerHalfStep_linear(ndim,lnode,dt,Nxyz,uxyz,uele)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: nval
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
    nval = npoin
    CALL pressureCalc(ndim,nval,uxyz,Pr)    

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

  END SUBROUTINE eulerHalfStep_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! solveEuler_compress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE solveEuler_compress(ndim,nface,lnode,bface,dt,elvol,xyz,Nxyz,uxyz,uele,duxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nface
    INTEGER :: nval
    INTEGER :: i,j,k,kstart,niter
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION :: dx,dy
    DOUBLE PRECISION :: Fx,Fy
    DOUBLE PRECISION :: F1,F2,F3,F4
    DOUBLE PRECISION :: Nax,Nay,Nbx,Nby,Ncx,Ncy
    DOUBLE PRECISION :: KLL,Lx,Ly
    DOUBLE PRECISION :: vx,vy,vz !velocity vector
    DOUBLE PRECISION :: rho      !density
    DOUBLE PRECISION :: en       !energy
    DOUBLE PRECISION :: dtg      !global min. time step
    DOUBLE PRECISION :: ML_inv,Ct
    DOUBLE PRECISION,PARAMETER :: third = 1.0d0/3.0d0,twelfth = 1.0d0/12.0d0
    DOUBLE PRECISION,PARAMETER :: beta1 = 1.0d0
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    TYPE(FACE),INTENT(IN),DIMENSION(nface) :: bface    !boundary faces
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: dt !element timestep array
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol ! Area of Element
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
    DOUBLE PRECISION,DIMENSION(nunks) :: Fbnd !boundary fluxes
    DOUBLE PRECISION,DIMENSION(npoin) :: ML !lumped mass matrix
    DOUBLE PRECISION,DIMENSION(npoin,npoin) :: MC !consistent mass matrix
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: duh !high order element contribution
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: MCduh !consistent mass matrix mult. by duh
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs !Right Hand Side for rho,rho*vx,rho*vy,rho*e
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs_Lap  !Lapidus Artificial vis.
    DOUBLE PRECISION,DIMENSION(nunks,npoin) :: rhs_FCT !Right Hand Side for Flux Corrected Transport

    Pr(:) = 0.0d0
    !calculate pressure
    nval = nelem
    CALL pressureCalc(ndim,nval,uele,Pr)    

    Fx = 0.0d0
    Fy = 0.0d0

    dtg = 1E+16
    dtp(:) = 1E+16
    duxyz(:,:) = 0.0d0
!!$    duh(:,:) = 0.0d0

    rhs(:,:) = 0.0d0
    rhs_Lap(:,:) = 0.0d0
!!$    rhs_FCT(:,:) = 0.0d0

    !initialize mass matrix
    ML(:) = 0.0d0
!!$    MC(:,:) = 0.0d0

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
       Nax = elvol(i)*Nax
       Nay = elvol(i)*Nay

       Nbx = elvol(i)*Nbx
       Nby = elvol(i)*Nby

       Ncx = elvol(i)*Ncx
       Ncy = elvol(i)*Ncy

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
       ML(na) = ML(na) + (elvol(i)/3.0d0)
       ML(nb) = ML(nb) + (elvol(i)/3.0d0)
       ML(nc) = ML(nc) + (elvol(i)/3.0d0)

       !Update Consistent Mass Matrix
!!$       kstart = 1
!!$       DO j=1,nnode
!!$          DO k=kstart,nnode
!!$             nj = lnode(j,i)
!!$             nk = lnode(k,i)
!!$             MC(nj,nk) = MC(nj,nk) + elvol(i)*twelfth
!!$             MC(nk,nj) = MC(nj,nk) + elvol(i)*twelfth
!!$          ENDDO
!!$          kstart = kstart + 1
!!$       ENDDO

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

!!$    DO i=1,npoin
!!$       dtg = MIN(dtg,dtp(i))
!!$    ENDDO

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
       rhs(:,na) = rhs(:,na) - Fbnd(:)
       rhs(:,nb) = rhs(:,nb) - Fbnd(:)

    ENDDO

    !calculate HEC
!!$    MCduh(:,:) = 0.0d0
!!$    duh(:,:) = 0.0d0
!!$    niter = 3
!!$    DO k=1,niter
!!$       DO i=1,npoin
!!$          DO j=1,nunks
!!$             MCduh(j,:) = MATMUL(MC(:,:),duh(j,:))
!!$          ENDDO
!!$          ML_inv = 1.0d0/ML(i)
!!$          duh(1,i) = (dtp(i)*rhs(1,i) - MCduh(1,i))*ML_inv + duh(1,i)
!!$          duh(2,i) = (dtp(i)*rhs(2,i) - MCduh(2,i))*ML_inv + duh(2,i)
!!$          duh(3,i) = (dtp(i)*rhs(3,i) - MCduh(3,i))*ML_inv + duh(3,i)
!!$          duh(4,i) = (dtp(i)*rhs(4,i) - MCduh(4,i))*ML_inv + duh(4,i)
!!$       ENDDO
!!$    ENDDO   
!!$
!!$    !compute ML - MC
!!$    DO i=1,npoin
!!$       MC(i,i) = ML(i) - MC(i,i)
!!$    ENDDO
!!$
!!$    !compute AEC
!!$    DO i=1,npoin
!!$       Ct = dtg/dtp(i)
!!$       duh(:,i) = Ct*uxyz(:,i) + duh(:,i)
!!$    ENDDO
!!$    DO i=1,nunks
!!$       rhs_FCT(i,:) = MATMUL(MC,duh(i,:))
!!$    ENDDO

    !calulate Lapidus artificial viscosity
    CALL lapidusArtVis(ndim,lnode,elvol,xyz,Nxyz,uxyz,rhs_Lap)

    rhs(:,:) = rhs(:,:) + rhs_Lap(:,:)

    !Solve for unknowns, calculate LEC
!!$    duxyz(:,:) = 0.0d0
    DO i=1,npoin
       ML_inv = 1.0d0/ML(i)
       duxyz(1,i) = dtp(i)*rhs(1,i)*ML_inv !+ rhs_FCT(1,i)
       duxyz(2,i) = dtp(i)*rhs(2,i)*ML_inv !+ rhs_FCT(2,i)
       duxyz(3,i) = dtp(i)*rhs(3,i)*ML_inv !+ rhs_FCT(3,i)
       duxyz(4,i) = dtp(i)*rhs(4,i)*ML_inv !+ rhs_FCT(4,i)
    ENDDO

  END SUBROUTINE solveEuler_compress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! lapidusArtVis_euler !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE lapidusArtVis_euler(ndim,lnode,elvol,xyz,Nxyz,uxyz,rhs_Lap)
    !Artificial Viscosity: Lapidus
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: i,j,k,kstart
    INTEGER :: na,nb,nc,nj,nk,INFO
    DOUBLE PRECISION,PARAMETER :: c1 = 5.0d0
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
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol ! Area of Element
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

       he2 = elvol(i)

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
       !sum over all nodes
       KLL =  NaL*(vx(na)*Lx + vy(na)*Ly) + NbL*(vx(nb)*Lx + vy(nb)*Ly) + NcL*(vx(nc)*Lx + vy(nc)*Ly)
       KLL = c1*he2*KLL

       !calc. for each unknown
       VLap(:) = NaL*uxyz(:,na) + NbL*uxyz(:,nb) + NcL*uxyz(:,nc)
       
       !Multiply by ABS(KLL)
       VLap(:) = elvol(i)*ABS(KLL)*VLap(:) 
       
       rhs_Lap(:,na) = rhs_Lap(:,na) - NaL*VLap(:)
       rhs_Lap(:,nb) = rhs_Lap(:,nb) - NbL*Vlap(:)
       rhs_Lap(:,nc) = rhs_Lap(:,nc) - NcL*Vlap(:)
!!$       IF(i .EQ. 0)THEN
!!$          WRITE(*,*)' for element: ',i
!!$          WRITE(*,*)' with nodes: ',na,nb,nc
!!$          WRITE(*,*)' unkno(na)= ',(uxyz(j,na),j=1,4)
!!$          WRITE(*,*)' unkno(nb)= ',(uxyz(j,nb),j=1,4)
!!$          WRITE(*,*)' unkno(nc)= ',(uxyz(j,nc),j=1,4)
!!$          WRITE(*,*)' Lx,Ly= ',Lx,Ly
!!$          WRITE(*,*)' Nax,Nay= ',Nax,Nay
!!$          WRITE(*,*)' Nbx,Nby= ',Nbx,Nby
!!$          WRITE(*,*)' Ncx,Ncy= ',Ncx,Ncy
!!$          WRITE(*,*)' NaL,NbL,NcL= ',NaL,NbL,NcL
!!$          WRITE(*,*)' KLL= ',KLL
!!$          WRITE(*,*)' elvol,Vlap= ',elvol(i),(Vlap(j),j=1,4)
!!$          WRITE(*,*)' rhs_Lap_a= ',(rhs_Lap(j,na),j=1,4)
!!$          WRITE(*,*)' rhs_Lap_b= ',(rhs_Lap(j,nb),j=1,4)
!!$          WRITE(*,*)' rhs_Lap_c= ',(rhs_Lap(j,nc),j=1,4)
!!$       ENDIF
    ENDDO

  END SUBROUTINE lapidusArtVis_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EULER BCs (COMPRESS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE eulerBCs_farfield(ndim,lnode,bnode,xyz,ulast,uxyz)
    !Farfield boundary conditions
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
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

  END SUBROUTINE eulerBCs_farfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Timestep FEM (EULER) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE timestepFEM_euler(ndim,lnode,elvol,xyz,Nxyz,uxyz,dt)
    !Calculate allowable timestep finite element method
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: nval
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
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nelem) :: elvol !element area
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz !xyz coordinates
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim*nnode,nelem) :: Nxyz !shape function derivatives
    DOUBLE PRECISION,INTENT(IN),DIMENSION(nunks,npoin) :: uxyz !unknowns at nodes
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(nelem) :: dt !timestep
    DOUBLE PRECISION,DIMENSION(npoin) :: Pr !pressure at nodes
    DOUBLE PRECISION,DIMENSION(nnode) :: lambda !characteristics

    Pr(:) = 0.0d0
    !calculate pressure
    nval = npoin
    CALL pressureCalc(ndim,nval,uxyz,Pr)    

    dt(:) = 0.0d0
    DO i=1,nelem

       lambda(:) = 0.0d0

       na = lnode(1,i)
       nb = lnode(2,i)
       nc = lnode(3,i)

       he = DSQRT(elvol(i))

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

  END SUBROUTINE timestepFEM_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sum of Squares !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  SUBROUTINE sumSq_euler(k_in,duxyz)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: k_in
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

    PRINT*,k_in,du2(1),du2(2),du2(3),du2(4)

  END SUBROUTINE sumSq_euler


END MODULE solveOps

