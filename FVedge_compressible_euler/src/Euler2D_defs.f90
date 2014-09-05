!!$ Program:      Euler2D_defs.f90
!!$ Author:       Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description:  Contains global parameters and defines
!!$               structures.


MODULE constants

  !defined
  INTEGER,PARAMETER :: nedel = 3 !number of edges per element
  INTEGER,PARAMETER :: nnoed = 2 !number of nodes per edge
  INTEGER,PARAMETER :: nfael = 3 !number of faces per element
  INTEGER,PARAMETER :: nnofa = 2 !number of nodes per face

  !read-in
  INTEGER :: nnode              !number of nodes on mesh
  INTEGER :: nnoel              !number of nodes per element
  INTEGER :: nelem              !number of elements
  INTEGER :: nbcs               !number of boundary conditions
  INTEGER :: nboun              !number of boundary points

  !define from data
  INTEGER :: nedge              !number of edges
  INTEGER :: nface              !number of faces
  INTEGER :: nunks              !number of unknowns
  INTEGER :: nesup              !number of elements surrounding points
  INTEGER :: npsup              !number of points surrounding points

  !calculation specific values set in main file
  INTEGER          :: max_steps !maximum number of time-steps
  DOUBLE PRECISION :: tf        !final time
  DOUBLE PRECISION :: gamma     !ratio of specific heat, Cp/Cv 
  DOUBLE PRECISION :: alpha     !angle of attack in radians
  DOUBLE PRECISION :: pgmin     !minimum gas pressure
  DOUBLE PRECISION :: Cr        !Courant number
  DOUBLE PRECISION :: Clap      !Lapidus constant for scaling artificial viscosity
  CHARACTER(80)    :: flux_type
  CHARACTER(80)    :: limiter_type
  CHARACTER(80)    :: prob_type
  CHARACTER(80)    :: step_type

  !parameters
  DOUBLE PRECISION,PARAMETER :: tol   = 1.0E-14
  DOUBLE PRECISION,PARAMETER :: pi    = (ATAN(1.0d0))*4.0d0
  DOUBLE PRECISION,PARAMETER :: third = 1.0d0/3.0d0     
  DOUBLE PRECISION,PARAMETER :: half  = 1.0d0/2.0d0     
  DOUBLE PRECISION,PARAMETER :: fifth = 1.0d0/5.0d0
  DOUBLE PRECISION,PARAMETER :: sixth = 1.0d0/6.0d0

  !grid data
  DOUBLE PRECISION :: xmin,xmax,ymin,ymax !bounding box
  DOUBLE PRECISION :: Lx,Ly               !domain length in x and y

  !!!!!!!!!!!!!!! HYDRO FAR-FIELD VALUES !!!!!!!!!!!!!!!
  DOUBLE PRECISION :: Minf  !mach number
  DOUBLE PRECISION :: pginf !gas prssure
  DOUBLE PRECISION :: csinf !speed of sound
  DOUBLE PRECISION :: dinf  !density
  DOUBLE PRECISION :: vinf  !velocity
  DOUBLE PRECISION :: vxinf !x-component of velocity
  DOUBLE PRECISION :: vyinf !y-component of velocity
  DOUBLE PRECISION :: v2inf !velocity squared
  DOUBLE PRECISION :: eninf !energy density


END MODULE constants


MODULE data_types

  !Data type for conserved variables
  TYPE CONS
     DOUBLE PRECISION :: d      !density
     DOUBLE PRECISION :: mx     !x-component momentum
     DOUBLE PRECISION :: my     !y-component momentum
     DOUBLE PRECISION :: en     !energy density
  END TYPE CONS   

  !Data type for primative variables
  TYPE PRIM
     DOUBLE PRECISION :: d      !density
     DOUBLE PRECISION :: vx     !x-component velocity
     DOUBLE PRECISION :: vy     !y-component velocity
     DOUBLE PRECISION :: pg     !gas pressure
  END TYPE PRIM   

  !Data type for gradients
  TYPE GRAD
     DOUBLE PRECISION,DIMENSION(4) :: dx     
     DOUBLE PRECISION,DIMENSION(4) :: dy     
  END TYPE GRAD

  !Data type for nodes
  TYPE NODE
     INTEGER                         :: nsuel   !number of surrounding elements
     INTEGER                         :: nsuno   !number of surrounding nodess
     DOUBLE PRECISION                :: ml      !lumped mass entry (a.k.a. dual volume for FV)
     DOUBLE PRECISION                :: x,y     !point coordinates
     DOUBLE PRECISION                :: dt      !local time-step
     DOUBLE PRECISION                :: nwsd    !half normal wave speed
     DOUBLE PRECISION,DIMENSION(4)   :: res     !residual
     INTEGER,DIMENSION(:),POINTER    :: suno    !list of surrounding nodes
     INTEGER,DIMENSION(:),POINTER    :: suel    !list of surrounding elements
     DOUBLE PRECISION,DIMENSION(2,2) :: lsq_inv !2x2 least sqaures matrix inverse
  END TYPE NODE
  
  !Data type for elements
  TYPE ELEM
     INTEGER                      :: nnoel !number of nodes per element
     INTEGER                      :: nsuel !number of surrounding elements
     DOUBLE PRECISION             :: vol   !element volume (area in 2D)
     DOUBLE PRECISION             :: xc,yc !element centered point coordinates
     INTEGER,DIMENSION(:),POINTER :: lnode !list of nodes
     INTEGER,DIMENSION(:),POINTER :: suel  !list of surrounding elements
  END TYPE ELEM

  !Data type for edges
  TYPE EDGE
     INTEGER                       :: ni,nj   !nodes of edge
     INTEGER                       :: e1,e2   !elements connected through edge
     DOUBLE PRECISION              :: sa      !magnitude of unit scaled area vector
     DOUBLE PRECISION              :: se      !magnitude of unit scaled edge vector
     DOUBLE PRECISION,DIMENSION(2) :: sav     !unit scaled area vector
     DOUBLE PRECISION,DIMENSION(2) :: sev     !unit scaled edge vector
  END TYPE EDGE

  !Data type for boundaries
  TYPE BOUN
     CHARACTER(80)                         :: bname   !name of boundary condition
     INTEGER                               :: btype   !type of boundary (0:wall;
                                                      !                  1:periodic in x; 
                                                      !                  2:periodic in y;
                                                      !                  3:periodic in x and y;
                                                      !                  4:farfield)
     INTEGER                               :: nbnode  !number of boundary nodes
     INTEGER                               :: nbface  !number of boundary faces
     INTEGER,DIMENSION(:),POINTER          :: bnode   !list boundary nodes     
     INTEGER,DIMENSION(:),POINTER          :: belem   !list boundary elements          
     INTEGER,DIMENSION(:),POINTER          :: bface   !list boundary faces
     DOUBLE PRECISION,DIMENSION(:),POINTER :: bnx     !x-component of outward normal
     DOUBLE PRECISION,DIMENSION(:),POINTER :: bny     !y-component of outward normal
     DOUBLE PRECISION,DIMENSION(:),POINTER :: bn      !magnitude of outward normal
     DOUBLE PRECISION,DIMENSION(:),POINTER :: fnx     !x-component of face outward normal
     DOUBLE PRECISION,DIMENSION(:),POINTER :: fny     !y-component of face outward normal
     DOUBLE PRECISION,DIMENSION(:),POINTER :: fn      !magnitude of face normal
  END TYPE BOUN   

END MODULE data_types

MODULE mesh_data

  USE constants
  USE data_types

  IMPLICIT NONE
  TYPE(NODE),DIMENSION(:),POINTER           :: nodes   !nodes
  TYPE(ELEM),DIMENSION(:),POINTER           :: elems   !elements
  TYPE(EDGE),DIMENSION(:),POINTER           :: edges   !edges
  TYPE(BOUN),DIMENSION(:),POINTER           :: bouns   !boundaries
  TYPE(CONS),DIMENSION(:),POINTER           :: ucons   !conservative variables
  TYPE(PRIM),DIMENSION(:),POINTER           :: wprim   !primative variables
  TYPE(GRAD),DIMENSION(:),POINTER           :: gradw   !gradient of primative variables

END MODULE mesh_data
