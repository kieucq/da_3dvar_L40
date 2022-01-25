!
! NOTE:     
!      This is the Lorenz 95 model that is used a based
!      model for my various data assimilation systems.
!      The model has 40 vars and time integration is
!      RK-4th scheme. Also, the tangetial model for this 
!      L95 model will also be returned in the form of
!      a matrix. Also, the adjoint model (w.r.t. the 
!      Euler metric) is computed as transpose of the
!      tangential model.
!
! HISTORY:  
!    - Feb 2, 2009: Created  
!    - Feb 4, 2009: add namelist option
!
! REFERENCE:
!    - Lorenz 1996: 
!    - Talagran and Courier 1986:
!    - Krishnamurti 1998: Numerical technique
!
! AUTHOR: 
!      Chanh Q. Kieu, Research Associate
!      Dept. of Atmospheric and Oceanic science
!      Univ. of Maryland, College Park, MD
!      email: kieucq@atmos.umd.edu
!
! COPYRIGHT: (C) 2009
!
!=================================================================
!
    PROGRAM L95_model
    IMPLICIT NONE
    INTEGER, PARAMETER      :: N  = 40              ! number of vars
    INTEGER, PARAMETER      :: nt = 1000            ! number of integration steps
    REAL                    :: F                    ! forcing term
    REAL, DIMENSION(N,nt)   :: x                    ! main model var
    REAL, DIMENSION(N)      :: x1,x2,x3             ! buffer var for RK integration
    REAL, DIMENSION(N)      :: rhs1,rhs2,rhs3,rhs4  ! buffer total right hand side forcing
    REAL, DIMENSION(N,N)    :: M                    ! tangetial model for L95 model
    REAL, DIMENSION(N,N)    :: MT                   ! adjoint model
    REAL                    :: dt                   ! model time step (non-dim)
    INTEGER                 :: i,j,k,loop,irec      ! indexing
    INTEGER                 :: debug                ! debuging 
    INTEGER                 :: restart              ! restart inverval
    CHARACTER*100           :: ofile,temc           ! output restart
    irec        = 1
!
! open input file for the control run with no da cycles
!
    OPEN(10,file='fsc.dat',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=N*4)
    OPEN(11,file='namelist.L40')
    OPEN(12,file='ana.dat',status='old')
    OPEN(13,file='bgd.dat')
!
! reading namelist of model parameters
!
    READ(11,*)temc,temc,debug
    READ(11,*)temc,temc,restart
    READ(11,*)temc,temc,F
    READ(11,*)temc,temc,dt
    PRINT*,'Reading input namelist and returns'
    PRINT*,'debug      = ',debug
    PRINT*,'restart    = ',restart
    PRINT*,'F          = ',F
    PRINT*,'dt         = ',dt
!
! reading the analysis input data
!
    read(12,*)(x(i,1),i=1,N)
!
! write out first input
!
    WRITE(10,rec=irec)(x(i,1),i=1,N)
    irec        = irec + 1
!
! integrate the model now
!
    time_loop: DO loop = 2,nt
     IF (debug.eq.1) PRINT*,'Loop at time',loop
!
! compute forcing and advance modelfor RK step 1
!
     x1(:)      = x(:,loop-1)
     call rhs(F,N,x1,rhs1)
     DO i       = 1,N
      x1(i)     = x(i,loop-1) + rhs1(i)*dt/2
     ENDDO
!
! compute forcing and advance modelfor RK step 2
!
     call rhs(F,N,x1,rhs2)
     DO i       = 1,N
      x2(i)     = x(i,loop-1) + rhs2(i)*dt/2
     ENDDO
!
! compute forcing and advance model for RK step 3
!
     call rhs(F,N,x2,rhs3)
     DO i       = 1,N
      x3(i)     = x(i,loop-1) + rhs3(i)*dt
     ENDDO
!
! compute forcing and advance model for RK step 3
! 
     call rhs(F,N,x3,rhs4)
     DO i       = 1,N
      x(i,loop) = x(i,loop-1) + (rhs1(i)+2*rhs2(i)+2*rhs3(i)+rhs4(i))*dt/6
     ENDDO
     IF (debug.eq.1) THEN
      PRINT*,(x(i,loop),i=1,N)
      READ*
     ENDIF
!
! prinout the re-start invertval for creating obs data
!
     IF (loop.eq.restart) THEN
      WRITE(13,*)(x(i,loop),i=1,N)
     ENDIF
!
! printout for viewing now
!
     WRITE(10,rec=irec)(x(i,loop),i=1,N)
     irec       = irec + 1
    ENDDO time_loop
    PRINT*,'Lorenz 40-var model finishes perfectly'
    END


    SUBROUTINE rhs(F,N,x,rs)
    IMPLICIT NONE
    INTEGER N,i
    REAL x(N),rs(N),F
    rs(1)       = x(N)*(x(2)-x(N-1))     - x(1) + F
    rs(2)       = x(1)*(x(3)-x(N))       - x(2) + F
    rs(N)       = x(N-1)*(x(1)-x(N-2))   - x(N) + F
    DO i        = 3,N-1
     rs(i)      = x(i-1)*(x(i+1)-x(i-2)) - x(i) + F
    ENDDO
    RETURN
    END

    SUBROUTINE tangential_model(M,x,N)
    IMPLICIT NONE
    INTEGER N,i,j,k
    REAL M(N,N),x(N)
    M           = 0.
    DO i        = 3,N-1
     M(i,i-2)   = -x(i-1)
     M(i,i-1)   = x(i+1) - x(i-2)
     M(i,i)     = -1
     M(i,i+1)   = x(i-1)
    ENDDO
!
! deal with  i = 1
!
    M(1,N-1)   = -x(N)
    M(1,N)     = x(2) - x(N-1)
    M(1,1)     = -1
    M(1,2)     = x(N)
!
! deal with i = 2
!
    M(2,N)     = -x(1)
    M(2,1)     = x(3) - x(N)
    M(2,2)     = -1
    M(2,3)     = x(1)
!
! deal with i = N
!
    M(N,N-2)   = -x(N-1)
    M(N,N-1)   = x(1) - x(N-2)
    M(N,N)     = -1
    M(N,1)     = x(N-1)
    RETURN
    END


    SUBROUTINE adjoint_model(MT,x,N)
    IMPLICIT NONE
    INTEGER N,i,j,k
    REAL MT(N,N),M(N,N),x(N)
    call tangential_model(M,x,N)
    MT         = transpose(M)
    RETURN
    END

