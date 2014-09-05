!!$ Program: Euler2D_solver.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines to calculate 2D Euler fluxes.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).


MODULE fluxes

private                        

!!$public :: fluxRoe2D_v1
public :: fluxRoe2D
public :: fluxRHLL2D
public :: fluxHLLC2D

CONTAINS

  !********************************************************************************
  !* -- Roe's Flux Function with entropy fix (version 1) ---
  !*
  !* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
  !* Schemes, Journal of Computational Physics, 43, pp. 357-372.
  !*
  !* NOTE: 3D version of this subroutine is available for download at
  !*       http://cfdbooks.com/cfdcodes.html
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
  !*           primR(1:5) = right state (rhoR, uR, vR, pR)
  !*               njk(2) = Face normal (L -> R). Must be a unit vector.
  !*
  !* Output:    flux(1:5) = numerical flux
  !*                  wsn = half the max wave speed
  !*                        (to be used for time step calculations)
  !* ------------------------------------------------------------------------------
  !*
  !********************************************************************************
  subroutine fluxRoe2D_v1(nfv,wl,wr,Froe,nwsd)
!!$  subroutine fluxRoe2D(nfv,wl,wr,Froe,nwsd)
    
    use constants   , only : half,fifth,gamma
    use state_variables
    
    implicit none
    
    !Input:
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: nfv    !face normal, nf=[nx, ny]
    DOUBLE PRECISION,DIMENSION(4),INTENT(IN) :: wl,wr !left/right primative states 
    
    !Output:
    DOUBLE PRECISION,INTENT(OUT)              :: nwsd
    DOUBLE PRECISION,DIMENSION(4),INTENT(OUT) :: Froe !roe fluxes
    
    !Local variables
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: nx,ny                  !normal vector
    DOUBLE PRECISION :: tx,ty                  !tangent vector: tx*nx+ty*ny = 0
    DOUBLE PRECISION :: dl,dr          !density
    DOUBLE PRECISION :: vxl,vyl,vxr,vyr        !velocity components.
    DOUBLE PRECISION :: pgl,pgr          !gas pressures
    DOUBLE PRECISION :: enl,enr          !energy densities
    DOUBLE PRECISION :: vnl,vnr,vtl,vtr     !normal and tangent velocities
    DOUBLE PRECISION :: csl,csr,hl,hr    !speeds of sound.
    DOUBLE PRECISION :: kel,ker,keroe          !kinetic energies
    DOUBLE PRECISION :: RT,droe,vxroe,vyroe,hroe,csroe,vnroe,vtroe   ! Roe-averages
    DOUBLE PRECISION :: dd,dvn,dvt,dpg  ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: LdU    ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: fl,fr,fdiss !left/right fluxes and dissipation
    DOUBLE PRECISION,DIMENSION(4) :: dws ! User-specified width for entropy fix
    DOUBLE PRECISION,DIMENSION(4) :: ws !normal wave speeds
    DOUBLE PRECISION,DIMENSION(4,4) :: rem !right-eigenmatrix
    

    froe(:) = 0.0d0
    fl(:) = 0.0d0
    fr(:) = 0.0d0
    fdiss(:) = 0.0d0

    nx = nfv(1)
    ny = nfv(2)
    
    !Tangent vector (Do you like it? Actually, Roe flux can be implemented 
    ! without any tangent vector. See "I do like CFD, VOL.1" for details.)
    tx = -ny
    ty =  nx
    
    !left state
    dl = wl(1)
    vxl = wl(2)
    vyl = wl(3)
    vnl = vxl*nx+vyl*ny
    vtl = vxl*tx+vyl*ty
    pgl = wl(4)
    enl = pgl/(gamma - 1.0d0) + dl*kel
    csl = dsqrt(gamma*pgl/dl)
    kel = half*(vxl*vxl + vyl*vyl)
    enl = pgl/(gamma - 1.0d0) + dl*kel
    hl  = (enl + pgl)/dl
!!$    hl  = csl*csl/(gamma-1.0d0) + kel
    
    !right state
    dr  = wr(1)
    vxr = wr(2)
    vyr = wr(3)
    vnr = vxr*nx+vyr*ny
    vtr = vxr*tx+vyr*ty
    pgr = wr(4)
    csr = dsqrt(gamma*pgr/dr)
    ker = half*(vxr*vxr + vyr*vyr)
    enr = pgr/(gamma - 1.0d0) + dr*ker
    hr  = (enr + pgr)/dr
!!$    hr  = csr*csr/(gamma-1.0d0) + ker
   
    !compute roe averages
    RT    = dsqrt(dr/dl)
    droe  = RT*dl
    vxroe = (vxl + RT*vxr) / (1.0d0 + RT)
    vyroe = (vyl + RT*vyr) / (1.0d0 + RT)
    keroe = half*(vxroe*vxroe + vyroe*vyroe)
    hroe  = (hl + RT*hr) / (1.0d0 + RT)
    csroe = dsqrt((gamma-1.0d0)*(hroe - keroe))
    vnroe = vxroe*nx + vyroe*ny
    vtroe = vxroe*tx + vyroe*ty
    
    !Wave Strengths
    dd  = dr - dl 
    dvn = vnr - vnl
    dvt = vtr - vtl
    dpg = pgr - pgl
    
    LdU(1) = (dpg - droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(2) = droe*dvt
    LdU(3) = dd - dpg/(csroe*csroe)
    LdU(4) = (dpg + droe*csroe*dvn )/(2.0d0*csroe*csroe)
    
    !Wave Speed
    ws(1) = abs(vnroe - csroe)
    ws(2) = abs(vnroe)
    ws(3) = abs(vnroe)
    ws(4) = abs(vnroe + csroe)
    
    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    ! only for the nonlinear fields.
    dws(1) = fifth
    if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
    dws(4) = fifth
    if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )
    
    !Right Eigenvectors
    rem(1,1) = 1.0d0    
    rem(2,1) = vxroe - csroe*nx
    rem(3,1) = vyroe - csroe*ny
    rem(4,1) = hroe - vnroe*csroe
    
    rem(1,2) = 0.0d0
    rem(2,2) = tx
    rem(3,2) = ty
    rem(4,2) = vtroe
    
    rem(1,3) = 1.0d0
    rem(2,3) = vxroe
    rem(3,3) = vyroe
    rem(4,3) = keroe
    
    rem(1,4) = 1.0d0
    rem(2,4) = vxroe + csroe*nx
    rem(3,4) = vyroe + csroe*ny
    rem(4,4) = hroe + vnroe*csroe
    
    !dissipation term
    fdiss(:) = 0.0d0
    do i=1,4
       do j=1,4
          fdiss(i) = fdiss(i) + ws(j)*LdU(j)*rem(i,j)
       end do
    end do
    
    !Compute the flux.
    fl(1) = dl*vnl
    fl(2) = dl*vnl * vxl + pgl*nx
    fl(3) = dl*vnl * vyl + pgl*ny
    fl(4) = dl*vnl * hl
    
    fr(1) = dr*vnr
    fr(2) = dr*vnr * vxr + pgr*nx
    fr(3) = dr*vnr * vyr + pgr*ny
    fr(4) = dr*vnr * hr

    Froe = half*(fl + fr - fdiss)
    nwsd = half*(abs(vnroe) + csroe)  !half normal max wave speed
    
  end subroutine fluxRoe2D_v1
!!$  end subroutine fluxRoe2D
  !--------------------------------------------------------------------------------

  !********************************************************************************
  !* -- Roe's Flux Function with entropy fix (version 2) ---
  !*    * this version does not use tangent vector to compute flux
  !*
  !* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
  !* Schemes, Journal of Computational Physics, 43, pp. 357-372.
  !*
  !* NOTE: 3D version of this subroutine is available for download at
  !*       http://cfdbooks.com/cfdcodes.html
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   wl(1:4) =  left state (dl, vxl, vyl, pgl)
  !*           wr(1:4) = right state (dr, vxr, vyr, pgr)
  !*            nfv(2) = Face normal (l -> r). Must be a unit vector.
  !*
  !* Output: froe(1:4) = numerical flux
  !*              nwsd = half the max wave speed
  !*                     (to be used for time step calculations)
  !* ------------------------------------------------------------------------------
  !*
  !********************************************************************************
!!$  subroutine fluxRoe2D_v2(nfv,wl,wr,froe,nwsd)
  subroutine fluxRoe2D(nfv,wl,wr,froe,nwsd)
    
    use constants   , only : half,fifth,gamma
    use state_variables
    
    implicit none
    
    !Input:
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: nfv    !face normal, nf=[nx, ny]
    DOUBLE PRECISION,DIMENSION(4),INTENT(IN) :: wl,wr !left/right primative states 
    
    !Output:
    DOUBLE PRECISION,INTENT(OUT)              :: nwsd
    DOUBLE PRECISION,DIMENSION(4),INTENT(OUT) :: froe !roe fluxes
    
    !Local variables
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: nx,ny                  !normal vector
    DOUBLE PRECISION :: tx,ty                  !tangent vector: tx*nx+ty*ny = 0
    DOUBLE PRECISION :: dl,dr          !density
    DOUBLE PRECISION :: vxl,vyl,vxr,vyr        !velocity components.
    DOUBLE PRECISION :: pgl,pgr          !gas pressures
    DOUBLE PRECISION :: enl,enr          !energy densities
    DOUBLE PRECISION :: vnl,vnr,vtl,vtr     !normal and tangent velocities
    DOUBLE PRECISION :: csl,csr,hl,hr    !speeds of sound.
    DOUBLE PRECISION :: kel,ker,keroe          !kinetic energies
    DOUBLE PRECISION :: RT,droe,vxroe,vyroe,hroe,csroe,vnroe,vtroe   ! Roe-averages
    DOUBLE PRECISION :: dd,dvn,dvt,dpg  ! Wave strenghs
    DOUBLE PRECISION :: dvx,dvy  ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: LdU    ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: fl,fr,fdiss !left/right fluxes and dissipation
    DOUBLE PRECISION,DIMENSION(4) :: dws ! User-specified width for entropy fix
    DOUBLE PRECISION,DIMENSION(4) :: ws !normal wave speeds
    DOUBLE PRECISION,DIMENSION(4,4) :: rem !right-eigenmatrix
    

    froe(:) = 0.0d0
    fl(:) = 0.0d0
    fr(:) = 0.0d0
    fdiss(:) = 0.0d0

    nx = nfv(1)
    ny = nfv(2)
        
    !left state
    dl = wl(1)
    vxl = wl(2)
    vyl = wl(3)
    vnl = vxl*nx+vyl*ny
!!$    vtl = vxl*tx+vyl*ty
    pgl = wl(4)
    csl = dsqrt(gamma*pgl/dl)
    kel = half*(vxl*vxl + vyl*vyl)
    enl = pgl/(gamma - 1.0d0) + dl*kel
    hl  = (enl + pgl)/dl
!!$    hl  = csl*csl/(gamma-1.0d0) + kel
    
    !right state
    dr  = wr(1)
    vxr = wr(2)
    vyr = wr(3)
    vnr = vxr*nx+vyr*ny
!!$    vtr = vxr*tx+vyr*ty
    pgr = wr(4)
    csr = dsqrt(gamma*pgr/dr)
    ker = half*(vxr*vxr + vyr*vyr)
    enr = pgr/(gamma - 1.0d0) + dr*ker
    hr  = (enr + pgr)/dr
!!$    hr  = csr*csr/(gamma-1.0d0) + ker
   
    !compute roe averages
    RT    = dsqrt(dr/dl)
    droe  = RT*dl
    vxroe = (vxl + RT*vxr) / (1.0d0 + RT)
    vyroe = (vyl + RT*vyr) / (1.0d0 + RT)
    keroe = half*(vxroe*vxroe + vyroe*vyroe)
    hroe  = (hl + RT*hr) / (1.0d0 + RT)
    csroe = dsqrt((gamma-1.0d0)*(hroe - keroe))
    vnroe = vxroe*nx + vyroe*ny
!!$    vtroe = vxroe*tx + vyroe*ty
    
    !Wave Strengths
    dd  = dr - dl 
    dvn = vnr - vnl
!!$    dvt = vtr - vtl
    dpg = pgr - pgl
    
    LdU(1) = (dpg - droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(2) = dd - dpg/(csroe*csroe)
    LdU(3) = (dpg + droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(4) = droe    

    !Wave Speed
    ws(1) = abs(vnroe - csroe)  !left moving acoustic wave
    ws(2) = abs(vnroe)          !enthropy wave
    ws(3) = abs(vnroe + csroe)  !right moving acoustic wave
    ws(4) = abs(vnroe)          !shear waves
    
    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    ! only for the nonlinear fields.
    dws(1) = fifth
    if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
    dws(3) = fifth
    if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )
    
    !Right Eigenvectors

    !left moving acoustic waves
    rem(1,1) = 1.0d0    
    rem(2,1) = vxroe - csroe*nx
    rem(3,1) = vyroe - csroe*ny
    rem(4,1) = hroe - vnroe*csroe

    !enthorpy wave
    rem(1,2) = 1.0d0
    rem(2,2) = vxroe
    rem(3,2) = vyroe
    rem(4,2) = keroe

    !right moving acoustic waves    
    rem(1,3) = 1.0d0
    rem(2,3) = vxroe + csroe*nx
    rem(3,3) = vyroe + csroe*ny
    rem(4,3) = hroe + vnroe*csroe

    !shear waves
    dvx = vxr - vxl
    dvy = vyr - vyl
    rem(1,4) = 0.0d0
    rem(2,4) = dvx - dvn*nx
    rem(3,4) = dvy - dvn*ny
    rem(4,4) = vxroe*dvx + vyroe*dvy - vnroe*dvn
    
    !dissipation term
    fdiss(:) = 0.0d0
    do i=1,4
       do j=1,4
          fdiss(i) = fdiss(i) + ws(j)*LdU(j)*rem(i,j)
       end do
    end do
    
    !Compute the flux.
    fl(1) = dl*vnl
    fl(2) = dl*vnl * vxl + pgl*nx
    fl(3) = dl*vnl * vyl + pgl*ny
    fl(4) = dl*vnl * hl
    
    fr(1) = dr*vnr
    fr(2) = dr*vnr * vxr + pgr*nx
    fr(3) = dr*vnr * vyr + pgr*ny
    fr(4) = dr*vnr * hr

    froe(:) = half*(fl(:) + fr(:) - fdiss(:))
    nwsd = half*(abs(vnroe) + csroe)  !half normal max wave speed
    
!!$  end subroutine fluxRoe2D_v2
  end subroutine fluxRoe2D
  !--------------------------------------------------------------------------------

  !*****************************************************************************
  !* -- Rotated-RHLL Flux Function ---
  !*
  !* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
  !* Resolving, Rotated-Hybrid Riemann Solvers,
  !* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
  !*
  !* Robust Riemann solver for nonlinear instability (carbuncle).
  !*
  !* NOTE: 3D version of this subroutine is available for download at
  !*       http://cfdbooks.com/cfdcodes.html
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
  !*           primR(1:5) = right state (rhoR, uR, vR, pR)
  !*               njk(2) = Face normal (L -> R). Must be a unit vector.
  !*
  !* Output:    flux(1:5) = numerical flux
  !*                  wsn = half the max wave speed
  !*                        (to be used for time step calculations)
  !* ------------------------------------------------------------------------------
  !*
  !*****************************************************************************
  subroutine fluxRHLL2D_v1(nfv,wl,wr,frhll,nwsd)

    use constants, only : half,fifth,gamma,Minf
    use state_variables

    IMPLICIT NONE

    !Input:
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: nfv(2) !normal face vector, nfv=[nx, ny]
    DOUBLE PRECISION,DIMENSION(4),INTENT(IN) :: wl,wr  !left/right primative states
    
    !Output:
    DOUBLE PRECISION,DIMENSION(4),INTENT(OUT) :: frhll !RHLL flux
    DOUBLE PRECISION,INTENT(OUT) :: nwsd               !normal wave speed

    !Local variables
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: nx,ny                  !normal vector
    DOUBLE PRECISION :: tx,ty                  !tangent vector: tx*nx+ty*ny = 0
    DOUBLE PRECISION :: dl,dr          !density
    DOUBLE PRECISION :: vxl,vyl,vxr,vyr        !velocity components.
    DOUBLE PRECISION :: pgl,pgr          !gas pressures
    DOUBLE PRECISION :: enl,enr          !energy densities
    DOUBLE PRECISION :: vnl,vnr,vtl,vtr     !normal and tangent velocities
    DOUBLE PRECISION :: csl,csr,hl,hr    !speeds of sound.
    DOUBLE PRECISION :: kel,ker,keroe          !kinetic energies
    DOUBLE PRECISION :: RT,droe,vxroe,vyroe,hroe,csroe,vnroe,vtroe   ! Roe-averages
    DOUBLE PRECISION :: dd,dvn,dvt,dpg  ! Wave strenghs
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
    DOUBLE PRECISION :: alpha1, alpha2                 ! Projections of the new normals
    DOUBLE PRECISION :: abs_dq                         ! Magnitude of the velocity difference
    DOUBLE PRECISION :: srp,slm                        ! Wave speeds for the HLL part
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION,DIMENSION(4) :: LdU    ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: fl,fr,fdiss !left/right fluxes and dissipation
    DOUBLE PRECISION,DIMENSION(4) :: dws ! User-specified width for entropy fix
    DOUBLE PRECISION,DIMENSION(4) :: ws !normal wave speeds
    DOUBLE PRECISION,DIMENSION(4) :: abs_ws 
    DOUBLE PRECISION,DIMENSION(4,4) :: rem !right-eigenmatrix
    
    nx = nfv(1)
    ny = nfv(2)
    
    !Tangent vector (Do you like it? Actually, Roe flux can be implemented 
    ! without any tangent vector. See "I do like CFD, VOL.1" for details.)
    tx = -ny
    ty =  nx

    !left state
    dl = wl(1)
    vxl = wl(2)
    vyl = wl(3)
    vnl = vxl*nx+vyl*ny
    vtl = vxl*tx+vyl*ty
    pgl = wl(4)
    csl = dsqrt(gamma*pgl/dl)
    kel = half*(vxl*vxl + vyl*vyl)
    enl = pgl/(gamma - 1.0d0) + dl*kel
    hl  = (enl + pgl)/dl
!!$    hl  = csl*csl/(gamma-1.0d0) + kel

    !right state
    dr  = wr(1)
    vxr = wr(2)
    vyr = wr(3)
    vnr = vxr*nx+vyr*ny
    vtr = vxr*tx+vyr*ty
    pgr = wr(4)
    csr = dsqrt(gamma*pgr/dr)
    ker = half*(vxr*vxr + vyr*vyr)
    enr = pgr/(gamma - 1.0d0) + dr*ker
    hr  = (enr + pgr)/dr
!!$    hr  = csr*csr/(gamma-1.0d0) + ker
       
    !calculate left/right fluxes
    fl(1) = dl*vnl
    fl(2) = dl*vnl * vxl + pgl*nx
    fl(3) = dl*vnl * vyl + pgl*ny
    fl(4) = dl*vnl * hl
    
    fr(1) = dr*vnr
    fr(2) = dr*vnr * vxr + pgr*nx
    fr(3) = dr*vnr * vyr + pgr*ny
    fr(4) = dr*vnr * hr

    !Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
    !(NB: n1 and n2 may need to be frozen at some point during 
    !     a steady calculation to fully make it converge. For time-accurate 
    !     calculation, this is fine.)
    ! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
    
    eps = (1.0E-12)*Minf
    abs_dq = dsqrt( (vxr-vxl)**2.0d0 + (vyr-vyl)**2.0d0 )
    
    if ( abs_dq > eps) then
       nx1 = (vxr-vxl)/abs_dq
       ny1 = (vyr-vyl)/abs_dq
    else
       nx1 = -ny 
       ny1 =  nx
    endif

    !Rey = 1000.0d0
    !temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
    !      /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
    !nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
    !ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)
    
    alpha1 = nx * nx1 + ny * ny1 
    !To make alpha1 always positive.
    temp = sign(1.0d0,alpha1)
    nx1 = temp * nx1
    ny1 = temp * ny1
    alpha1 = temp * alpha1
    
    !Take n2 as perpendicular to n1.
    nx2 = -ny1
    ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
    !To make alpha2 always positive.
    temp = sign(1.0d0,alpha2)
    nx2 = temp * nx2
    ny2 = temp * ny2
    alpha2 = temp * alpha2
    
    !Now we are going to compute the Roe flux with n2 as the normal
    !and n1 as the tagent vector, with modified wave speeds (5.12)
    
    !compute roe averages
    RT    = dsqrt(dr/dl)
    droe  = RT*dl
    vxroe = (vxl + RT*vxr) / (1.0d0 + RT)
    vyroe = (vyl + RT*vyr) / (1.0d0 + RT)
    keroe = half*(vxroe*vxroe + vyroe*vyroe)
    hroe  = (hl + RT*hr) / (1.0d0 + RT)
    csroe = dsqrt((gamma-1.0d0)*(hroe - keroe))
    vnroe = vxroe*nx2 + vyroe*ny2
    vtroe = vxroe*nx1 + vyroe*ny1

    !Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnl = vxl*nx2 + vyl*ny2
    vnr = vxr*nx2 + vyr*ny2
    vtl = vxl*nx1 + vyl*ny1
    vtr = vxr*nx1 + vyr*ny1
    
!!$print*,nx1,ny1,nx2,ny2
!!$print*,vnl,vnr,vtl,vtr

    dd = dr - dl 
    dvn =  vnr - vnl
    dvt =  vtr - vtl
    dpg =  pgr - pgl

    LdU(1) = (dpg - droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(2) = droe*dvt
    LdU(3) = dd - dpg/(csroe*csroe)
    LdU(4) = (dpg + droe*csroe*dvn )/(2.0d0*csroe*csroe)
    
    !Wave Speeds
    ws(1) = vnroe - csroe
    ws(2) = vnroe
    ws(3) = vnroe
    ws(4) = vnroe + csroe
    abs_ws  = abs(ws)

    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    !only for the nonlinear fields.
    dws(1) = fifth
    if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
    dws(4) = fifth
    if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))    
    
    !HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    srp = max(0.0d0,vtr + csr,vtroe + csroe)
    slm = min(0.0d0,vtl - csl,vtroe - csroe)
    
    !Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
    ws = alpha2*abs_ws - ( alpha2*(srp+slm)*ws + 2.0d0*alpha1*srp*slm )/ (srp-slm)
    
    !Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1
    ty = ny1
    
    !Right Eigenvectors
    rem(1,1) = 1.0d0    
    rem(2,1) = vxroe - csroe*nx2
    rem(3,1) = vyroe - csroe*ny2
    rem(4,1) = hroe - vnroe*csroe

    rem(1,2) = 0.0d0
    rem(2,2) = tx
    rem(3,2) = ty
    rem(4,2) = vtroe

    rem(1,3) = 1.0d0
    rem(2,3) = vxroe
    rem(3,3) = vyroe
    rem(4,3) = keroe
        
    rem(1,4) = 1.0d0
    rem(2,4) = vxroe + csroe*nx2
    rem(3,4) = vyroe + csroe*ny2
    rem(4,4) = hroe + vnroe*csroe

    !dissipation term
    fdiss(:) = 0.0d0
    do i=1,4
       do j=1,4
          fdiss(i) = fdiss(i) + ws(j)*LdU(j)*rem(i,j)
       end do
    end do
    
    !Compute the Rotated-RHLL flux.
    frhll(:) = (srp*fl(:) - slm*fr(:))/(srp-slm) - half*fdiss(:)

    nwsd = half*( abs(vnroe) + abs(vtroe) + csroe)  !Normal max wave speed times half
    
  end subroutine fluxRHLL2D_v1
  !--------------------------------------------------------------------------------

  !*****************************************************************************
  !* -- Rotated-RHLL Flux Function ---
  !*
  !* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
  !* Resolving, Rotated-Hybrid Riemann Solvers,
  !* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
  !*
  !* Robust Riemann solver for nonlinear instability (carbuncle).
  !*
  !* NOTE: 3D version of this subroutine is available for download at
  !*       http://cfdbooks.com/cfdcodes.html
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
  !*           primR(1:5) = right state (rhoR, uR, vR, pR)
  !*               njk(2) = Face normal (L -> R). Must be a unit vector.
  !*
  !* Output:    flux(1:5) = numerical flux
  !*                  wsn = half the max wave speed
  !*                        (to be used for time step calculations)
  !* ------------------------------------------------------------------------------
  !*
  !*****************************************************************************
  subroutine fluxRHLL2D(nfv,wl,wr,frhll,nwsd)

    use constants, only : half,fifth,gamma,Minf
    use state_variables

    IMPLICIT NONE

    !Input:
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: nfv(2) !normal face vector, nfv=[nx, ny]
    DOUBLE PRECISION,DIMENSION(4),INTENT(IN) :: wl,wr  !left/right primative states
    
    !Output:
    DOUBLE PRECISION,DIMENSION(4),INTENT(OUT) :: frhll !RHLL flux
    DOUBLE PRECISION,INTENT(OUT) :: nwsd               !normal wave speed

    !Local variables
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: nx,ny                  !normal vector
    DOUBLE PRECISION :: tx,ty                  !tangent vector: tx*nx+ty*ny = 0
    DOUBLE PRECISION :: dl,dr          !density
    DOUBLE PRECISION :: vxl,vyl,vxr,vyr        !velocity components.
    DOUBLE PRECISION :: pgl,pgr          !gas pressures
    DOUBLE PRECISION :: enl,enr          !energy densities
    DOUBLE PRECISION :: vnl,vnr,vtl,vtr     !normal and tangent velocities
    DOUBLE PRECISION :: csl,csr,hl,hr    !speeds of sound.
    DOUBLE PRECISION :: kel,ker,keroe          !kinetic energies
    DOUBLE PRECISION :: RT,droe,vxroe,vyroe,hroe,csroe,vnroe,vtroe   ! Roe-averages
    DOUBLE PRECISION :: dd,dvn,dvt,dpg  ! Wave strenghs
    DOUBLE PRECISION :: dvx,dvy         ! Wave strenghs
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
    DOUBLE PRECISION :: alpha1, alpha2                 ! Projections of the new normals
    DOUBLE PRECISION :: abs_dq                         ! Magnitude of the velocity difference
    DOUBLE PRECISION :: srp,slm                        ! Wave speeds for the HLL part
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION,DIMENSION(4) :: LdU    ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: fl,fr,fdiss !left/right fluxes and dissipation
    DOUBLE PRECISION,DIMENSION(4) :: dws ! User-specified width for entropy fix
    DOUBLE PRECISION,DIMENSION(4) :: ws !normal wave speeds
    DOUBLE PRECISION,DIMENSION(4) :: abs_ws 
    DOUBLE PRECISION,DIMENSION(4,4) :: rem !right-eigenmatrix
    
    nx = nfv(1)
    ny = nfv(2)
    
    !Tangent vector (Do you like it? Actually, Roe flux can be implemented 
    ! without any tangent vector. See "I do like CFD, VOL.1" for details.)
    tx = -ny
    ty =  nx

    !left state
    dl = wl(1)
    vxl = wl(2)
    vyl = wl(3)
    vnl = vxl*nx+vyl*ny
    vtl = vxl*tx+vyl*ty
    pgl = wl(4)
    csl = dsqrt(gamma*pgl/dl)
    kel = half*(vxl*vxl + vyl*vyl)
    enl = pgl/(gamma - 1.0d0) + dl*kel
    hl  = (enl + pgl)/dl
!!$    hl  = csl*csl/(gamma-1.0d0) + kel

    !right state
    dr  = wr(1)
    vxr = wr(2)
    vyr = wr(3)
    vnr = vxr*nx+vyr*ny
    vtr = vxr*tx+vyr*ty
    pgr = wr(4)
    csr = dsqrt(gamma*pgr/dr)
    ker = half*(vxr*vxr + vyr*vyr)
    enr = pgr/(gamma - 1.0d0) + dr*ker
    hr  = (enr + pgr)/dr
!!$    hr  = csr*csr/(gamma-1.0d0) + ker
       
    !calculate left/right fluxes
    fl(1) = dl*vnl
    fl(2) = dl*vnl * vxl + pgl*nx
    fl(3) = dl*vnl * vyl + pgl*ny
    fl(4) = dl*vnl * hl
    
    fr(1) = dr*vnr
    fr(2) = dr*vnr * vxr + pgr*nx
    fr(3) = dr*vnr * vyr + pgr*ny
    fr(4) = dr*vnr * hr

    !Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
    !(NB: n1 and n2 may need to be frozen at some point during 
    !     a steady calculation to fully make it converge. For time-accurate 
    !     calculation, this is fine.)
    ! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
    
    eps = (1.0E-12)*Minf
    abs_dq = dsqrt( (vxr-vxl)**2.0d0 + (vyr-vyl)**2.0d0 )
    
    if ( abs_dq > eps) then
       nx1 = (vxr-vxl)/abs_dq
       ny1 = (vyr-vyl)/abs_dq
    else
       nx1 = -ny 
       ny1 =  nx
    endif

    !Rey = 1000.0d0
    !temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
    !      /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
    !nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
    !ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)
    
    alpha1 = nx * nx1 + ny * ny1 
    !To make alpha1 always positive.
    temp = sign(1.0d0,alpha1)
    nx1 = temp * nx1
    ny1 = temp * ny1
    alpha1 = temp * alpha1
    
    !Take n2 as perpendicular to n1.
    nx2 = -ny1
    ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
    !To make alpha2 always positive.
    temp = sign(1.0d0,alpha2)
    nx2 = temp * nx2
    ny2 = temp * ny2
    alpha2 = temp * alpha2
    
    !Now we are going to compute the Roe flux with n2 as the normal
    !and n1 as the tagent vector, with modified wave speeds (5.12)
    
    !compute roe averages
    RT    = dsqrt(dr/dl)
    droe  = RT*dl
    vxroe = (vxl + RT*vxr) / (1.0d0 + RT)
    vyroe = (vyl + RT*vyr) / (1.0d0 + RT)
    keroe = half*(vxroe*vxroe + vyroe*vyroe)
    hroe  = (hl + RT*hr) / (1.0d0 + RT)
    csroe = dsqrt((gamma-1.0d0)*(hroe - keroe))
    vnroe = vxroe*nx2 + vyroe*ny2
    vtroe = vxroe*nx1 + vyroe*ny1

    !Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnl = vxl*nx2 + vyl*ny2
    vnr = vxr*nx2 + vyr*ny2
    vtl = vxl*nx1 + vyl*ny1
    vtr = vxr*nx1 + vyr*ny1
    
!!$print*,nx1,ny1,nx2,ny2
!!$print*,vnl,vnr,vtl,vtr

    dd = dr - dl 
    dvn =  vnr - vnl
    dvt =  vtr - vtl
    dpg =  pgr - pgl

    LdU(1) = (dpg - droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(2) = dd - dpg/(csroe*csroe)
    LdU(3) = (dpg + droe*csroe*dvn )/(2.0d0*csroe*csroe)
    LdU(4) = droe    

    !Wave Speed
    ws(1) = vnroe - csroe !left moving acoustic wave
    ws(2) = vnroe          !enthropy wave
    ws(3) = vnroe + csroe  !right moving acoustic wave
    ws(4) = vnroe          !shear waves
    abs_ws  = abs(ws)

    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    !only for the nonlinear fields.
    dws(1) = fifth
    if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
    dws(3) = fifth
    if (abs_ws(3)<dws(3)) abs_ws(3) = half*(abs_ws(3)*abs_ws(3)/dws(3)+dws(3))    
    
    !HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    srp = max(0.0d0,vtr + csr,vtroe + csroe)
    slm = min(0.0d0,vtl - csl,vtroe - csroe)
    
    !Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
    ws = alpha2*abs_ws - ( alpha2*(srp+slm)*ws + 2.0d0*alpha1*srp*slm )/ (srp-slm)
    
    !Right Eigenvectors: with n2 as normal and n1 as tangent.
!!$    tx = nx1
!!$    ty = ny1
    
    !Right Eigenvectors

    !left moving acoustic waves
    rem(1,1) = 1.0d0    
    rem(2,1) = vxroe - csroe*nx2
    rem(3,1) = vyroe - csroe*ny2
    rem(4,1) = hroe - vnroe*csroe

    !enthorpy wave
    rem(1,2) = 1.0d0
    rem(2,2) = vxroe
    rem(3,2) = vyroe
    rem(4,2) = keroe

    !right moving acoustic waves    
    rem(1,3) = 1.0d0
    rem(2,3) = vxroe + csroe*nx2
    rem(3,3) = vyroe + csroe*ny2
    rem(4,3) = hroe + vnroe*csroe

    !shear waves
    dvx = vxr - vxl
    dvy = vyr - vyl
    rem(1,4) = 0.0d0
    rem(2,4) = dvx - dvn*nx2
    rem(3,4) = dvy - dvn*ny2
    rem(4,4) = vxroe*dvx + vyroe*dvy - vnroe*dvn


    !dissipation term
    fdiss(:) = 0.0d0
    do i=1,4
       do j=1,4
          fdiss(i) = fdiss(i) + ws(j)*LdU(j)*rem(i,j)
       end do
    end do
    
    !Compute the Rotated-RHLL flux.
    frhll(:) = (srp*fl(:) - slm*fr(:))/(srp-slm) - half*fdiss(:)

    nwsd = half*( abs(vnroe) + abs(vtroe) + csroe)  !Normal max wave speed times half
    
  end subroutine fluxRHLL2D
  !--------------------------------------------------------------------------------

  !*****************************************************************************
  !* -- HLLC (Harten, Lax, Van Leer, with contact wave) Flux Function ---
  !*
  !* T. Miyoshi & K. Kusano, "A multi-state HLL approximate Riemann solver
  !* for ideal MHD", JCP, 208, 315 (2005)
  !*
  !* A Riemann solver capable of resolving contact waves.
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
  !*           primR(1:5) = right state (rhoR, uR, vR, pR)
  !*               njk(2) = Face normal (L -> R). Must be a unit vector.
  !*
  !* Output:    flux(1:5) = numerical flux
  !*                  wsn = half the max wave speed
  !*                        (to be used for time step calculations)
  !* ------------------------------------------------------------------------------
  !*
  !*****************************************************************************
  subroutine fluxHLLC2D(nfv,wl,wr,fhllc,nwsd)

    use constants, only : half,fifth,gamma,Minf
    use state_variables

    IMPLICIT NONE

    !Input:
    DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: nfv(2) !normal face vector, nfv=[nx, ny]
    DOUBLE PRECISION,DIMENSION(4),INTENT(IN) :: wl,wr  !left/right primative states
    
    !Output:
    DOUBLE PRECISION,DIMENSION(4),INTENT(OUT) :: fhllc !HLLC flux
    DOUBLE PRECISION,INTENT(OUT) :: nwsd               !normal wave speed

    !Local variables
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: beta                  !parameter
    DOUBLE PRECISION :: nx,ny                  !normal vector
    DOUBLE PRECISION :: tx,ty                  !tangent vector: tx*nx+ty*ny = 0
    DOUBLE PRECISION :: dl,dr                  !density
    DOUBLE PRECISION :: vxl,vyl,vxr,vyr        !velocity components.
    DOUBLE PRECISION :: vxl0,vyl0        !velocity components.
    DOUBLE PRECISION :: pgl,pgr                !gas pressures
    DOUBLE PRECISION :: enl,enr                !energy densities
    DOUBLE PRECISION :: vnl,vnr,vtl,vtr !normal and tangent velocities
    DOUBLE PRECISION :: csl,csr,hl,hr   !speeds of sound.
    DOUBLE PRECISION :: kel,ker,keroe   !kinetic energies
    DOUBLE PRECISION :: RT,droe,vxroe,vyroe,hroe,csroe,vnroe,vtroe ! Roe-averages
    DOUBLE PRECISION :: dd,dvn,dvt,dpg  ! Wave strenghs
    DOUBLE PRECISION :: dvx,dvy         ! Wave strenghs
    DOUBLE PRECISION :: sl,sr,sm        ! signal strengths
    DOUBLE PRECISION :: dsl,dsr         !intermediate densities
    DOUBLE PRECISION :: vxsl,vxsr,vysl,vysr,vns !intermediate velocites
    DOUBLE PRECISION :: vnsl,vnsr,vtsl,vtsr !intermediate velocites
    DOUBLE PRECISION :: pgs       !intermediate gas pressure
    DOUBLE PRECISION :: ensl,ensr !intermediate energy densities
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: srp,slm        ! Wave speeds for the HLL part
    DOUBLE PRECISION :: al,ar,am
    DOUBLE PRECISION :: bp,bm
    DOUBLE PRECISION :: tl,tr
    DOUBLE PRECISION :: dmr,dml
    DOUBLE PRECISION :: pgc     !contact gas pressure
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION,DIMENSION(4) :: LdU ! Wave strenghs
    DOUBLE PRECISION,DIMENSION(4) :: fl,fr,fdiss !left/right fluxes and dissipation
    DOUBLE PRECISION,DIMENSION(4) :: dws ! User-specified width for entropy fix
    DOUBLE PRECISION,DIMENSION(4) :: ev  ! eigenvalues
    DOUBLE PRECISION,DIMENSION(4) :: ws !normal wave speeds
    DOUBLE PRECISION,DIMENSION(4) :: abs_ws 
    DOUBLE PRECISION,DIMENSION(4,4) :: rem !right-eigenmatrix
    
    !normal vetor
    nx = nfv(1)
    ny = nfv(2)
    
    !tangent vector 
    tx = -ny
    ty =  nx

    !left state
    dl = wl(1)
    vxl = wl(2)
    vyl = wl(3)
    vnl = vxl*nx+vyl*ny
    vtl = vxl*tx+vyl*ty
    pgl = wl(4)
    csl = dsqrt(gamma*pgl/dl)
    kel = half*(vxl*vxl + vyl*vyl)
    hl  = csl*csl/(gamma-1.0d0) + kel
    enl = pgl/(gamma - 1.0d0) + dl*kel

    !right state
    dr  = wr(1)
    vxr = wr(2)
    vyr = wr(3)
    vnr = vxr*nx+vyr*ny
    vtr = vxr*tx+vyr*ty
    pgr = wr(4)
    csr = dsqrt(gamma*pgr/dr)
    ker = half*(vxr*vxr + vyr*vyr)
    hr  = csr*csr/(gamma-1.0d0) + ker
    enr = pgr/(gamma - 1.0d0) + dr*ker

    !for subsonic cases
    beta = max(sqrt(2.0d0*kel)/csl,sqrt(2.0d0*ker)/csr)
    beta = min(1.0d0,beta)
    vxl0 = vxl
    vyl0 = vyl
    vxl = half*(vxl + vxr) + beta*half*(vxl - vxr)
    vyl = half*(vyl + vyr) + beta*half*(vyl - vyr)
    vxr = half*(vxl0 + vxr) - beta*half*(vxl0 - vxr)
    vyr = half*(vyl0 + vyr) - beta*half*(vyl0 - vyr)

    !compute roe averages, used to update normal wave speed (nwsd)
    RT    = dsqrt(dr/dl)
    droe  = RT*dl
    vxroe = (vxl + RT*vxr) / (1.0d0 + RT)
    vyroe = (vyl + RT*vyr) / (1.0d0 + RT)
    keroe = half*(vxroe*vxroe + vyroe*vyroe)
    hroe  = (hl + RT*hr) / (1.0d0 + RT)
    csroe = dsqrt((gamma-1.0d0)*(hroe - keroe))
    vnroe = vxroe*nx + vyroe*ny
    vtroe = vxroe*tx + vyroe*ty

    !eigenvalues
    ev(1) = vnroe - csroe
    ev(2) = vnroe
    ev(3) = vnroe + csroe
    ev(4) = vnroe        

    !Wave Speed
    ws(1) = abs(vnroe - csroe)  !left moving acoustic wave
    ws(2) = abs(vnroe)          !enthropy wave
    ws(3) = abs(vnroe + csroe)  !right moving acoustic wave
    ws(4) = abs(vnroe)          !shear waves
    
    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    ! only for the nonlinear fields.
    dws(1) = fifth
    if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
    dws(3) = fifth
    if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

    !calculate left/right signal speeds
    sl = min(ev(1),vnl - csl)
    sr = max(ev(3),vnr + csr)

    !middle signal speed
    !eq. (38) of Miyoshi and Kusano
    sm = ((sr - vnr)*dr*vnr - (sl - vnl)*dl*vnl - pgr + pgl)/((sr - vnr)*dr - (sl - vnl)*dl)

    !calculate left/right fluxes
!!$    fl(1) = dl*vnl
!!$    fl(2) = dl*vnl * vxl + pgl*nx
!!$    fl(3) = dl*vnl * vyl + pgl*ny
!!$    fl(4) = dl*vnl * hl
!!$    
!!$    fr(1) = dr*vnr
!!$    fr(2) = dr*vnr * vxr + pgr*nx
!!$    fr(3) = dr*vnr * vyr + pgr*ny
!!$    fr(4) = dr*vnr * hr

    fl(1) = dl*vnl
    fl(2) = dl*vnl * vnl + pgl
    fl(3) = dl*vnl * vtl + pgl
    fl(4) = dl*vnl * hl
    
    fr(1) = dr*vnr
    fr(2) = dr*vnr * vnr + pgr
    fr(3) = dr*vnr * vtr + pgr
    fr(4) = dr*vnr * hr

    !calculate left/right signal speeds
!!$    sl = min(vnl,vnr) - max(csl,csr)
!!$    sr = max(vnl,vnr) + max(csl,csr)

    !middle signal speed
    !eq. (38) of Miyoshi and Kusano
!!$    sm = ((sr - vnr)*dr*vnr - (sl - vnl)*dl*vnl - pgr + pgl)/((sr - vnr)*dr - (sl - vnl)*dl)

    !calulate frist intemediate state

    !pressure and normal velocity are constant over Riemann fan
    !eq. (23) and (39) of Miyoshi and Kusano
    pgs = pgl + dl*(sl - vnl)*(sm - vnl)
    vns = sm

    !eq. (43) of Miyoshi and Kusano
    dsl = dl*(sl - vnl)/(sl - sm)
    dsr = dr*(sr - vnr)/(sr - sm)

    !eqs. (44) and (46) of Miyoshi and Kusano revisted for general geometry
    !vxl^* = vxl + (sl - sm)*(pg^* - pgl)*nx/(dl*(sl - vnl)*(sl - sm))
    !vyl^* = vyl + (sl - sm)*(pg^* - pgl)*ny/(dl*(sl - vnl)*(sl - sm))
    !if (nx,ny) = (1,0), then vx^* = sm and vy^* = vy just as in
    !Miyoshi and Kusano.  Equations simplified below
!!$    vxsl = vxl + (sm - vnl)*nx
!!$    vxsr = vxr + (sm - vnr)*nx
!!$
!!$    vysl = vyl + (sm - vnl)*ny
!!$    vysr = vyr + (sm - vnr)*ny

!!$    vxsl = vxl + (sl - sm)*(pgs - pgl)*nx/(dl*(sl - vnl)*(sl - sm))
!!$    vxsr = vxr + (sr - sm)*(pgs - pgr)*nx/(dr*(sr - vnr)*(sr - sm))
!!$
!!$    vysl = vyl + (sl - sm)*(pgs - pgl)*ny/(dl*(sl - vnl)*(sl - sm))
!!$    vysr = vyr + (sr - sm)*(pgs - pgr)*ny/(dr*(sr - vnr)*(sr - sm))

    vnsl = vnl + (sl - sm)*(pgs - pgl)*nx/(dl*(sl - vnl)*(sl - sm))
    vnsr = vnr + (sr - sm)*(pgs - pgr)*nx/(dr*(sr - vnr)*(sr - sm))

    vtsl = vtl + (sl - sm)*(pgs - pgl)*ny/(dl*(sl - vnl)*(sl - sm))
    vtsr = vtr + (sr - sm)*(pgs - pgr)*ny/(dr*(sr - vnr)*(sr - sm))

    !eq. (48) of Miyoshi and Kusano    
    ensl = ((sl - vnl)*enl - pgl*vnl + pgs*sm)/(sl - sm)
    ensr = ((sr - vnr)*enr - pgr*vnl + pgs*sm)/(sr - sm)

    !calculate interface/intercell fluxes
    if(sl .le. 0)then
      fhllc(:) = fl(:)
    elseif((sl .lt. 0) .and. (0 .le. sm))then
      fhllc(1) = fl(1) + sl*(dsl - dl)
      fhllc(2) = fl(2) + sl*(dsl*vnsl - dl*vnl)
      fhllc(3) = fl(3) + sl*(dsl*vtsl - dl*vtl)
      fhllc(4) = fl(4) + sl*(ensl - enl)
    elseif((sm .lt. 0) .and. (0 .le. sr))then
      fhllc(1) = fr(1) + sr*(dsr - dl)
      fhllc(2) = fr(2) + sr*(dsr*vnsr - dr*vnr)
      fhllc(3) = fr(3) + sr*(dsr*vtsr - dr*vtr)
      fhllc(4) = fr(4) + sr*(ensr - enr)
    elseif(sr .lt. 0)then
      fhllc(:) = fr(:) 
    endif

    !rotate velocities
    fhllc(2) = fhllc(2)*nx + fhllc(3)*tx
    fhllc(3) = fhllc(2)*ny + fhllc(3)*ty

    nwsd = half*( abs(vnroe) + csroe)  !Normal max wave speed times half

!!$    al = max(ev(1),vnl + csl)
!!$    ar = max(ev(3),vnr + csr)
!!$    
!!$    bp = 0.0d0
!!$    bm = 0.0d0
!!$    if(al .lt. 0.0d0)then
!!$       bm = al
!!$    endif
!!$    if(ar .gt. 0.0d0)then
!!$       bp = ar
!!$    endif
!!$
!!$    tl = pgl +(vnl - al)*dl*vnl
!!$    tr = pgr +(vnr - ar)*dr*vnr
!!$
!!$    dml = dl*vnl - dl*al
!!$    dmr = -(dr*vnr - dr*ar)
!!$
!!$    am = (tl - tr)/(dml + dmr)
!!$
!!$    !contact gas pressure
!!$    pgc = (dml*tr + dmr*tl)/(dml + dmr)
!!$    pgc = max(pgc,0.0d0)
!!$
!!$    !calculate left/right fluxes
!!$    fl(1) = dl*vnl - bm*dl
!!$    fl(2) = dl*vxl * (vnl - bm) + pgl*nx
!!$    fl(3) = dl*vyl * (vnl - bm) + pgl*ny
!!$    fl(4) = enl*(vnl - bm) + pgl*vnl
!!$
!!$    fr(1) = dr*vnr - bp*dr
!!$    fr(2) = dr*vxr * (vnr - bp) + (pgr + sm*pgc)*nx
!!$    fr(3) = dr*vyr * (vnr - bp) + (pgr + sm*pgc)*ny
!!$    fr(4) = enr*(vnr - bp) + pgr*vnr + sm*am*pgc
!!$    
!!$    if(am .ge. 0.0d0)then
!!$       sl = am/(am - bm)
!!$       sr = 0.0d0
!!$       sm = -bm/(am - bm)
!!$    else
!!$       sl = 0.0d0
!!$       sr = -am/(bp - am)
!!$       sm = bp/(bp - am)
!!$    endif
!!$
!!$    fhllc(:) = sl*fl(:) + sr*fr(:)
!!$
!!$    nwsd = half*( abs(vnroe) + csroe)  !Normal max wave speed times half

    
  end subroutine fluxHLLC2D
  !--------------------------------------------------------------------------------

END MODULE fluxes
