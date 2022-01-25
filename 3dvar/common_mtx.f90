MODULE common_mtx
!=======================================================================
!
! [PURPOSE:] Matrix Functions
!
! [CREATED:] 07/20/2004 Takemasa Miyoshi
! [UPDATED:] 10/16/2004 Takemasa Miyoshi
!
! [PUBLIC:]
!   mtx_eigen  : eigenvalue decomposition
!   mtx_inv    : real symmetric matrix inverse
!   mtx_sqrt   : real symmetric matrix square root
!
! [REFERENCES:]
!    Core subroutines are adapted from netlib.org
!
! [HISTORY:]
!  07/20/2003 Takemasa Miyoshi  Created at University of Maryland, College Park
!
!=======================================================================
  USE common

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mtx_eigen, mtx_inv, mtx_sqrt

CONTAINS
!=======================================================================
!  Eigenvalue decomposition using subroutine rs
!    INPUT
!      INTEGER :: imode           : mode switch (0: only eiven values)
!      INTEGER :: n               : dimension of matrix
!      REAL(r_size) :: a(n,n)     : input matrix
!    OUTPUT
!      REAL(r_size) :: eival(n)   : eiven values in decending order
!                                   i.e. eival(1) is the largest
!      REAL(r_size) :: eivec(n,n) : eiven vectors
!      INTEGER :: nrank_eff       : number of positive eivenvalues
!=======================================================================
SUBROUTINE mtx_eigen(imode,n,a,eival,eivec,nrank_eff)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: imode ! 0: calculate only eigen values
  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: eival(1:n)
  REAL(r_size),INTENT(OUT) :: eivec(1:n,1:n)
  INTEGER,INTENT(OUT) :: nrank_eff

  REAL(r_dble),ALLOCATABLE :: a8(:,:)
  REAL(r_dble),ALLOCATABLE :: eival8(:)
  REAL(r_dble),ALLOCATABLE :: eivec8(:,:)
  REAL(r_dble),ALLOCATABLE :: wrk1(:)
  REAL(r_dble),ALLOCATABLE :: wrk2(:)
  INTEGER :: ierr,i,j

  ALLOCATE( a8(1:n,1:n) )
  a8 = a
  ALLOCATE( eival8(1:n) )
  ALLOCATE( eivec8(1:n,1:n) )
  ALLOCATE( wrk1(1:n),wrk2(1:n) )
  eivec8 = 0.0d0
  CALL rs(n,n,a8,eival8,imode,eivec8,wrk1,wrk2,ierr)
  IF( ierr/=0 ) THEN
    PRINT *,'!!! ERROR (mtx_eigen): rs error code is ',ierr
    STOP 2
  END IF
  DEALLOCATE( wrk1,wrk2,a8 )

  nrank_eff = n
  IF( eival8(n) > 0 ) THEN
    DO i=1,n
      IF( eival8(i) < ABS(eival8(n))*SQRT(EPSILON(eival8)) ) THEN
        nrank_eff = nrank_eff - 1
        eival8(i) = 0.0d0
        eivec8(:,i) = 0.0d0
      END IF
    END DO
  ELSE
    PRINT *,'!!! ERROR (mtx_eigen): All Eigenvalues are below 0'
    STOP 2
  END IF

  IF( nrank_eff<n .AND. eival8(1)/=0 ) THEN
    j = 0
    DO i=n,1,-1
      IF( eival8(i) == 0 ) THEN
        eival8(i) = eival8(n-nrank_eff-j)
        eivec(:,i) = eivec8(:,n-nrank_eff-j)
        eival8(n-nrank_eff-j) = 0.0d0
        eivec8(:,n-nrank_eff-j) = 0.0d0
        j = j+1
      END IF
    END DO
  END IF

  DO i=1,n
    eival(i) = eival8(n+1-i)
    eivec(:,i) = eivec8(:,n+1-i)
  END DO
  DEALLOCATE( eival8,eivec8 )

  RETURN
END SUBROUTINE mtx_eigen
!=======================================================================
!  Real symmetric matrix inversion using subroutine dspdi
!    INPUT
!      INTEGER :: n               : dimension of matrix
!      REAL(r_size) :: a(n,n)     : input matrix (real symmetric)
!    OUTPUT
!      REAL(r_size) :: ainv(n,n)  : inverse of a
!=======================================================================
SUBROUTINE mtx_inv(n,a,ainv)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: ainv(1:n,1:n)

  REAL(8),ALLOCATABLE :: acmp(:)
  REAL(8),ALLOCATABLE :: det(:)
  REAL(8),ALLOCATABLE :: work(:)
  INTEGER,ALLOCATABLE :: kpvt(:)
  INTEGER,ALLOCATABLE :: inert(:)
  INTEGER :: info
  INTEGER :: i,j,k

  IF(n==1) THEN
    ainv(1,1) = 1.0d0 / a(1,1)
  ELSE

  ALLOCATE( acmp(1:n*(n+1)/2) )
  ALLOCATE( det(1:2) )
  ALLOCATE( work(1:n) )
  ALLOCATE( kpvt(1:n) )
  ALLOCATE( inert(1:3) )
!-----------------------------------------------------------------------
!  Packed form of matrix
!-----------------------------------------------------------------------
  k=0
  DO j=1,n
    DO i=1,j
      k = k+1
      acmp(k) = a(i,j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  dspfa
!-----------------------------------------------------------------------
  CALL dspfa(acmp,n,kpvt,info)
  IF(info /= 0) THEN
    PRINT *,'!!! ERROR (mtx_inv): dspfa error code is ',info
    STOP 3
  END IF
!-----------------------------------------------------------------------
!  dspdi
!-----------------------------------------------------------------------
  CALL dspdi(acmp,n,kpvt,det,inert,work,001)
!-----------------------------------------------------------------------
!  unpack matrix
!-----------------------------------------------------------------------
  k=0
  DO j=1,n
    DO i=1,j
      k = k+1
      ainv(i,j) = acmp(k)
    END DO
  END DO

  DO j=1,n
    DO i=j+1,n
      ainv(i,j) = ainv(j,i)
    END DO
  END DO

  DEALLOCATE( acmp,det,work,kpvt,inert )

  END IF

  RETURN
END SUBROUTINE mtx_inv
!=======================================================================
!  Compute square root of real symmetric matrix
!    INPUT
!      INTEGER :: n                : dimension of matrix
!      REAL(r_size) :: a(n,n)      : input matrix (real symmetric)
!    OUTPUT
!      REAL(r_size) :: a_sqrt(n,n) : square root of a
!=======================================================================
SUBROUTINE mtx_sqrt(n,a,a_sqrt)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: n
  REAL(r_size),INTENT(IN) :: a(1:n,1:n)
  REAL(r_size),INTENT(OUT) :: a_sqrt(1:n,1:n)

  REAL(r_size),ALLOCATABLE :: eival(:)   ! holds eivenvalue of a
  REAL(r_size),ALLOCATABLE :: eivec(:,:) ! holds eivenvector of a
  INTEGER :: i,n_eff

  ALLOCATE( eival(1:n) )
  ALLOCATE( eivec(1:n,1:n) )
  CALL mtx_eigen(1,n,a,eival,eivec,n_eff)

  DO i=1,n
    a_sqrt(:,i) = eivec(:,i) * SQRT( eival(i) )
  END DO
  DEALLOCATE( eival )

  a_sqrt = matmul(a_sqrt,transpose(eivec))

  DEALLOCATE( eivec )

  RETURN
END SUBROUTINE mtx_sqrt

END MODULE common_mtx
