!
! This program is for creating an ensemble of input 
! for a cold start run of the breeding method. These
! input is supposed to be analysis at t = 0, so you
! dont need to run the assimilation step at t = 0.
! By default, there will have 11 members created.
!
  PROGRAM cold_start_ini 
  IMPLICIT NONE
  INTEGER,PARAMETER :: nx = 51, ny = 21 
  INTEGER,PARAMETER :: ne = 11 
  REAL,PARAMETER    :: obs_err_psi = 1.0d5
  REAL              :: psi(nx,ny),psip(nx,ny),psim(nx,ny)
  REAL(8)           :: rnd(nx*ny),psir(nx,ny)
  CHARACTER*50      :: ofile
  INTEGER           :: id,i,j,k
  OPEN(10,file='baroin.dat',status='old')
  READ(10,'(6e13.6)')((psi(i,j),i=1,nx),j=1,ny)
  ofile     = '../ana/ana00_00.dat'
!
! creating an ensemble of +/- enhanced fields.
! It turns out that this is the germ of singulariy
! for the matrxi inversion later used by etkf
!
  DO k      = 1,ne
   PRINT*,'Creating the ensemble ',k
   CALL com_randn(nx*ny,rnd)
   psir     = RESHAPE(rnd,(/nx,ny/))*obs_err_psi
   psip     = psi + psir
   id       = k 
   IF (id.le.9) THEN
    WRITE(ofile(12:12),'(1I1)')id
   ELSE
    WRITE(ofile(11:12),'(1I2)')id
   ENDIF
   OPEN(11,file=ofile,status='unknown')
   WRITE(11,'(6e13.6)')((psip(i,j),i=1,nx),j=1,ny)
   CLOSE(11)
  ENDDO
  END 


  SUBROUTINE com_randn(ndim,var)
  USE mt19937
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  REAL(8),INTENT(OUT) :: var(1:ndim)
  REAL(8) :: rnd(2),pi
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.
  pi = 4.*atan(1.)

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_genrand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF
  RETURN
END SUBROUTINE com_randn
